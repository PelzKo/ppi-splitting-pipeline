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
    // GO annotations are not part of this mode (they describe proteins, not
    // domain families), so unlike DATA_PREP's three-way gate the precomputed
    // hatch asks only for the files that cannot be derived.
    branched = datasets_ch.branch { meta, ddis, sequences, species, domain_instances ->
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

    // Pfam-A.clans.tsv is a release-wide reference file consumed by the one
    // shared task, so it is a run-global param rather than a samplesheet column
    // -- same idiom as gurobi_license in SPLIT_POSITIVES.
    clans_ch = params.pfam_clans
        ? channel.value(file(params.pfam_clans, checkIfExists: true))
        : channel.value([])

    shared_fetch   = FETCH_DOMAIN_META(families_list.map { families -> tuple([id: "_shared"], families) }, clans_ch)
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
