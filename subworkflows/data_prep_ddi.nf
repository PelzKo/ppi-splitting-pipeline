include { EMPTY_GO_ANNOTATIONS; FETCH_DOMAIN_META; GET_LENGTHS as GET_LENGTHS_DDI; GET_LENGTHS as GET_LENGTHS_SHARED_DDI; SUBSET_DOMAIN_DATA } from '../processes/data_prep'

// DATA_PREP's DDI-mode counterpart: turns the Pfam family accessions in the
// interaction file's protein1/protein2 columns into sampled domain instances
// (sequences, species, lengths) plus instances.tsv, the table mapping every
// instance back to its family, clan and parent protein.
//
// Same shape as DATA_PREP -- datasets needing a fetch are deduplicated into one
// union of family ids, fetched once under a synthetic "_shared" meta, then
// split back out per dataset -- and it emits the same channel names, which is
// what lets main.nf pick between the two with a branch instead of a second
// pipeline. `instances` is the one addition; DATA_PREP emits it as [].
workflow DATA_PREP_DDI {
    take:
    datasets_ch  // tuple(meta, files) -- reads files.ppis (Pfam family accessions here), .sequences, .species, .domain_instances

    main:
    // The precomputed hatch is all-or-nothing, so a row that supplies two of the
    // three files would fall into needs_fetch, ignore both, and stream 6.3 GB
    // from EBI -- a correct result by a route nobody asked for, hours later and
    // silently. Fail loudly instead, wording it like main.nf's --split_only check.
    checked_ch = datasets_ch.map { meta, f ->
        def supplied = [sequences: f.sequences, species: f.species, domain_instances: f.domain_instances]
        def missing  = supplied.findAll { k, v -> !v }.keySet()
        if (missing && missing.size() < supplied.size()) {
            error("--ddi_mode precomputed data prep needs sequences, species and domain_instances together (row '${meta.id}' supplies ${(supplied.keySet() - missing).join(', ')} but is missing ${missing.join(', ')}). Leave all three blank to fetch them from Pfam instead.")
        }
        tuple(meta, f)
    }

    // GO annotations are not part of this mode (they describe proteins, not
    // domain families), so unlike DATA_PREP's three-way gate the precomputed
    // hatch asks only for the files that cannot be derived.
    branched = checked_ch.branch { meta, f ->
        precomputed: f.sequences && f.species && f.domain_instances
            return tuple(meta, f.sequences, f.species, f.domain_instances)
        needs_fetch: true
            return tuple(meta, f.ppis)
    }

    precomputed_sequences = branched.precomputed.map { meta, sequences, species, domain_instances -> tuple(meta, sequences) }
    precomputed_species   = branched.precomputed.map { meta, sequences, species, domain_instances -> tuple(meta, species) }
    precomputed_instances = branched.precomputed.map { meta, sequences, species, domain_instances -> tuple(meta, domain_instances) }
    precomputed_lengths   = GET_LENGTHS_DDI(precomputed_sequences)
    precomputed_go        = EMPTY_GO_ANNOTATIONS(branched.precomputed.map { meta, sequences, species, domain_instances -> meta })

    // Same accession-pooling flatMap as DATA_PREP: DDI input reuses the
    // protein1/protein2 column names, the values are just Pfam accessions.
    families_list = branched.needs_fetch
        .flatMap { meta, ddis -> ddis.splitCsv(header: true).collectMany { row -> [row.protein1.trim(), row.protein2.trim()] } }
        .unique()
        .collectFile(name: 'families.txt', newLine: true, sort: true)

    // Pfam-A.clans.tsv and Pfam-A.fasta.gz are release-wide reference files
    // consumed by the one shared task, so they are run-global params rather than
    // samplesheet columns -- same idiom as gurobi_license in SPLIT_POSITIVES.
    // Both are staged as path inputs, which is what gets them checkIfExists and a
    // place in the task hash, and what makes them work on any executor.
    clans_ch = params.pfam_clans
        ? channel.value(file(params.pfam_clans, checkIfExists: true))
        : channel.value([])

    regions_ch = params.pfam_regions
        ? channel.value(file(params.pfam_regions, checkIfExists: true))
        : channel.value([])

    // The protein universe: Pfam-A.regions carries coordinates only, so sequence,
    // taxon and the review flag all come from UniProt flat files. Several may be
    // given -- Swiss-Prot first, then a TrEMBL file if an unreviewed stratum is in
    // --instance_tiers -- as a List or a comma-separated string. Staged like the
    // other two so they land in the task hash and get checkIfExists.
    def dat_raw = !params.uniprot_dat
        ? []
        : (params.uniprot_dat instanceof List ? params.uniprot_dat : params.uniprot_dat.toString().split(','))
    def dat_files = dat_raw.collect { it.toString().trim() }.findAll { it }.collect { file(it, checkIfExists: true) }
    dat_ch = channel.value(dat_files)

    // The cache is the exception: FETCH_DOMAIN_META *writes* to it, so it cannot
    // be staged, and a relative path would resolve inside the task's work dir and
    // be deleted with it -- making the documented "--interpro_cache makes re-runs
    // cheap" quietly untrue. Absolutised here. It must also be a filesystem every
    // compute node shares; on node-local scratch each task gets its own cold cache.
    cache_dir = params.interpro_cache ? file(params.interpro_cache).toAbsolutePath().toString() : ''
    if (cache_dir && !file(cache_dir).exists()) {
        log.warn "--interpro_cache ${cache_dir} does not exist yet; FETCH_DOMAIN_META will create it."
    }

    // FETCH_DOMAIN_META reads --instance_tiers straight from params (so an
    // embedding pipeline needs no wiring for it), which puts the value past
    // channel construction and into fetch_domains.py's argparse. Check it here
    // instead: a typo should not survive as far as a scheduled task. The name list
    // is duplicated from fetch_domains.py's TIERS -- Groovy cannot read it, and the
    // two must be changed together.
    def validTiers = ['human_reviewed', 'other_reviewed', 'human_unreviewed', 'other_unreviewed']
    def tiersRaw = params.instance_tiers instanceof List ? params.instance_tiers : params.instance_tiers.toString().split(',')
    def tiers = tiersRaw.collect { it.toString().trim() }.findAll { it }
    if (!tiers) {
        error("--instance_tiers is empty; name at least one of ${validTiers.join(', ')}.")
    }
    def unknownTiers = tiers - validTiers
    if (unknownTiers) {
        error("--instance_tiers has unknown stratum name(s) ${unknownTiers.join(', ')}; valid names are ${validTiers.join(', ')}.")
    }
    log.info "--instance_tiers ${tiers.join(',')}: families with nothing in a named stratum keep zero instances, and every DDI touching one drops out (see _shared/data/dropped_families.tsv)."
    // A warning, not an error, for the same reason the uniprot_dat guard below is
    // absent: this block runs for every DDI-mode workflow, including one that
    // supplies precomputed instances and never reaches FETCH_DOMAIN_META. On the
    // path that does need it, fetch_domains.py fails hard before the regions pass.
    if (tiers.any { it.endsWith('unreviewed') } && !dat_files) {
        log.warn "--instance_tiers names an unreviewed stratum but no --uniprot_dat was given; the Swiss-Prot download can only fill the reviewed strata and FETCH_DOMAIN_META will fail. Pass a TrEMBL flat file too."
    }
    // No guard here on --uniprot_dat / --interpro_cache. This block runs on *every*
    // DDI-mode workflow, including an embedding pipeline that supplies precomputed
    // instances and never reaches FETCH_DOMAIN_META at all -- for which the flat
    // file is irrelevant and demanding it would be a false failure. fetch_domains.py
    // raises the same error itself, on the only path that actually needs the file.

    shared_fetch   = FETCH_DOMAIN_META(
        families_list.map { families -> tuple([id: "_shared"], families) },
        clans_ch,
        regions_ch,
        dat_ch,
        channel.value(cache_dir),
    )
    shared_lengths = GET_LENGTHS_SHARED_DDI(shared_fetch.sequences)

    subset_out = SUBSET_DOMAIN_DATA(
        branched.needs_fetch,
        shared_fetch.sequences.map       { meta, f -> f }.first(),
        shared_fetch.go_annotations.map  { meta, f -> f }.first(),
        shared_fetch.species.map         { meta, f -> f }.first(),
        shared_lengths.map                { meta, f -> f }.first(),
        shared_fetch.instances.map       { meta, f -> f }.first(),
    )

    emit:
    sequences      = precomputed_sequences.mix(subset_out.sequences)
    go_annotations = precomputed_go.mix(subset_out.go_annotations)
    species        = precomputed_species.mix(subset_out.species)
    lengths        = precomputed_lengths.mix(subset_out.lengths)
    instances      = precomputed_instances.mix(subset_out.instances)
    // The one artefact that explains a shrunken DDI count, so an embedding
    // pipeline gets a channel rather than a hardcoded results/_shared/data/ path.
    // Two properties to respect: the fetch runs once per run under the synthetic
    // [id: "_shared"] meta, so this carries at most one item whose meta matches no
    // dataset -- join() it against a per-dataset channel and you get zero items --
    // and it is legitimately empty when every dataset supplied precomputed domain
    // data, so nothing may block on it.
    dropped_families = shared_fetch.dropped_families
}
