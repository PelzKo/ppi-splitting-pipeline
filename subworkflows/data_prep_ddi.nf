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
    datasets_ch  // tuple(meta, ddis, sequences, species, domain_instances)

    main:
    // The precomputed hatch is all-or-nothing, so a row that supplies two of the
    // three files would fall into needs_fetch, ignore both, and stream 6.3 GB
    // from EBI -- a correct result by a route nobody asked for, hours later and
    // silently. Fail loudly instead, wording it like main.nf's --split_only check.
    checked_ch = datasets_ch.map { meta, ddis, sequences, species, domain_instances ->
        def supplied = [sequences: sequences, species: species, domain_instances: domain_instances]
        def missing  = supplied.findAll { k, v -> !v }.keySet()
        if (missing && missing.size() < supplied.size()) {
            error("--ddi_mode precomputed data prep needs sequences, species and domain_instances together (row '${meta.id}' supplies ${(supplied.keySet() - missing).join(', ')} but is missing ${missing.join(', ')}). Leave all three blank to fetch them from Pfam instead.")
        }
        tuple(meta, ddis, sequences, species, domain_instances)
    }

    // GO annotations are not part of this mode (they describe proteins, not
    // domain families), so unlike DATA_PREP's three-way gate the precomputed
    // hatch asks only for the files that cannot be derived.
    branched = checked_ch.branch { meta, ddis, sequences, species, domain_instances ->
        precomputed: sequences && species && domain_instances
            return tuple(meta, sequences, species, domain_instances)
        needs_fetch: true
            return tuple(meta, ddis)
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

    fasta_ch = params.pfam_fasta
        ? channel.value(file(params.pfam_fasta, checkIfExists: true))
        : channel.value([])

    // The cache is the exception: FETCH_DOMAIN_META *writes* to it, so it cannot
    // be staged, and a relative path would resolve inside the task's work dir and
    // be deleted with it -- making the documented "--interpro_cache makes re-runs
    // cheap" quietly untrue. Absolutised here. It must also be a filesystem every
    // compute node shares; on node-local scratch each task gets its own cold cache.
    cache_dir = params.interpro_cache ? file(params.interpro_cache).toAbsolutePath().toString() : ''
    if (cache_dir && !file(cache_dir).exists()) {
        log.warn "--interpro_cache ${cache_dir} does not exist yet; FETCH_DOMAIN_META will create it."
    }

    shared_fetch   = FETCH_DOMAIN_META(
        families_list.map { families -> tuple([id: "_shared"], families) },
        clans_ch,
        fasta_ch,
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
}
