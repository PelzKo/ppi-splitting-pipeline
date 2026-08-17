#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { samplesheetToList } from 'plugin/nf-schema'

include { DATA_PREP }        from './subworkflows/data_prep'
include { DATA_PREP_DDI }    from './subworkflows/data_prep_ddi'
include { CLUSTERING }       from './subworkflows/clustering'
include { SPLIT_POSITIVES }  from './subworkflows/split_positives'
include { SAMPLE_NEGATIVES } from './subworkflows/sample_negatives'
include { TRAIN_BASELINE }   from './subworkflows/train_baseline'
include { QC }               from './subworkflows/qc'

// samplesheetToList() represents a blank cell as [] (not null), even for
// numeric fields where 0 is a legitimate override -- so check both.
def isGiven(v) {
    !(v == null || v == [])
}

// One row per PPI dataset. Anything left blank in the samplesheet falls
// back to the corresponding default in nextflow.config, so a single run
// can process several datasets in parallel, each with its own overrides.
def buildDatasetsChannel() {
    // samplesheetToList() returns each row as a positional list, not a map
    // -- order here must match assets/schema_input.json's `properties`.
    def fields = [
        "id", "ppis", "sequences", "go_annotations", "species", "domain_instances", "blast_results", "candidate_network",
        "partition", "node_mapping",
        "embedding_model", "cdhit_identity", "cdhit_wordsize", "split_method", "edge_weight",
        "kahip_k", "ilp_kahip_k", "train_split", "val_split", "test_split", "ilp_epsilon", "ilp_max_sec",
        "negative_sampling_method",
        "neg_ilp_time_limit", "neg_ilp_lambda_degree",
        "neg_ilp_lambda_taxon_pair", "neg_ilp_lambda_self_loop", "neg_ilp_lambda_jaccard",
    ]
    // moduleDir, not projectDir: when this pipeline is included as a subworkflow
    // projectDir is the *root* project's directory, so the schema would be looked
    // up in the embedding pipeline's tree. Only the standalone entry below reaches
    // this line today, but the trap is silent and free to close.
    def rows = samplesheetToList(params.samplesheet, "${moduleDir}/assets/schema_input.json")

    return channel.fromList(rows).map { rowList ->
        def row = [fields, rowList].transpose().collectEntries { k, v -> [(k): v] }

        // split_only skips FETCH_DATA/CLUSTERING/TRAIN_BASELINE/QC entirely,
        // so every one of those steps' precomputed-input escape hatches
        // becomes mandatory, and the split/negative-sampling method choice
        // is no longer per-dataset -- it's always the ILP path.
        if (params.split_only) {
            // DDI mode has no GO annotations at all (they describe proteins, not
            // domain families) and needs the domain instance table instead.
            def required = params.ddi_mode
                ? ["sequences", "species", "domain_instances", "partition", "node_mapping"]
                : ["sequences", "go_annotations", "species", "partition", "node_mapping"]
            def missing = required.findAll { !row[it] }
            if (missing) {
                error("--split_only requires every samplesheet row to supply ${required.join(', ')} (row '${row.id}' is missing ${missing.join(', ')}).")
            }
        }

        def meta = [
            id                       : row.id,
            embedding_model          : isGiven(row.embedding_model)          ? row.embedding_model          : params.embedding_model,
            cdhit_identity           : isGiven(row.cdhit_identity)           ? row.cdhit_identity           : params.cdhit_identity,
            cdhit_wordsize           : isGiven(row.cdhit_wordsize)           ? row.cdhit_wordsize           : params.cdhit_wordsize,
            split_method             : params.split_only ? "ilp" : (isGiven(row.split_method)             ? row.split_method             : params.split_method),
            edge_weight              : isGiven(row.edge_weight)              ? row.edge_weight              : params.edge_weight,
            kahip_k                  : isGiven(row.kahip_k)                  ? row.kahip_k                  : params.kahip_k,
            ilp_kahip_k              : isGiven(row.ilp_kahip_k)              ? row.ilp_kahip_k              : params.ilp_kahip_k,
            train_split              : isGiven(row.train_split)              ? row.train_split              : params.train_split,
            val_split                : isGiven(row.val_split)                ? row.val_split                : params.val_split,
            test_split               : isGiven(row.test_split)               ? row.test_split               : params.test_split,
            ilp_epsilon              : isGiven(row.ilp_epsilon)              ? row.ilp_epsilon              : params.ilp_epsilon,
            ilp_max_sec              : isGiven(row.ilp_max_sec)              ? row.ilp_max_sec              : params.ilp_max_sec,
            negative_sampling_method : params.split_only ? "ilp" : (isGiven(row.negative_sampling_method) ? row.negative_sampling_method : params.negative_sampling_method),
            neg_ilp_time_limit       : isGiven(row.neg_ilp_time_limit)       ? row.neg_ilp_time_limit       : params.neg_ilp_time_limit,
            neg_ilp_lambda_degree    : isGiven(row.neg_ilp_lambda_degree)     ? row.neg_ilp_lambda_degree     : params.neg_ilp_lambda_degree,
            neg_ilp_lambda_taxon_pair: isGiven(row.neg_ilp_lambda_taxon_pair) ? row.neg_ilp_lambda_taxon_pair : params.neg_ilp_lambda_taxon_pair,
            neg_ilp_lambda_self_loop : isGiven(row.neg_ilp_lambda_self_loop)  ? row.neg_ilp_lambda_self_loop  : params.neg_ilp_lambda_self_loop,
            // The Jaccard bias term matches GO-term overlap, and DDI mode has no
            // GO annotations at all (they describe proteins, not domain families
            // -- DATA_PREP_DDI emits a header-only table). 0 makes
            // assemble_active_biases() skip the term outright rather than
            // computing an all-zero one and reporting it as perfectly matched,
            // so it overrides the samplesheet the same way --split_only
            // overrides split_method above.
            neg_ilp_lambda_jaccard   : params.ddi_mode ? 0 : (isGiven(row.neg_ilp_lambda_jaccard) ? row.neg_ilp_lambda_jaccard : params.neg_ilp_lambda_jaccard),
        ]
        // A named map, not a positional tail: consumers say `f.blast_results`
        // rather than counting commas, so adding a samplesheet column touches
        // the `fields` list above and this literal, and nothing else. An absent
        // optional file is [] and never null -- a `path` input accepts [] as
        // "no file", and null would blow up staging.
        def files = [
            ppis             : file(row.ppis, checkIfExists: true),
            sequences        : row.sequences         ? file(row.sequences,         checkIfExists: true) : [],
            go_annotations   : row.go_annotations    ? file(row.go_annotations,    checkIfExists: true) : [],
            species          : row.species           ? file(row.species,           checkIfExists: true) : [],
            domain_instances : row.domain_instances  ? file(row.domain_instances,  checkIfExists: true) : [],
            blast_results    : row.blast_results     ? file(row.blast_results,     checkIfExists: true) : [],
            candidate_network: row.candidate_network ? file(row.candidate_network, checkIfExists: true) : [],
            partition        : row.partition         ? file(row.partition,         checkIfExists: true) : [],
            node_mapping     : row.node_mapping      ? file(row.node_mapping,      checkIfExists: true) : [],
        ]
        tuple(meta, files)
    }
}

// The whole pipeline, as a named workflow so another project can `include` it.
// Nextflow can only import named workflows, and only from a local path -- the
// anonymous `workflow { }` at the bottom is the standalone entry point and the
// only thing that reads params.samplesheet.
workflow PPI_SPLITTING {
    take:
    datasets_ch   // tuple(meta, filesMap); filesMap keys: ppis, sequences, go_annotations,
                  // species, domain_instances, blast_results, candidate_network, partition,
                  // node_mapping. An absent optional file is [], never null. meta carries
                  // `id` plus every per-dataset parameter override -- see
                  // buildDatasetsChannel() for the full list, and note that mutating meta
                  // downstream breaks the joins, which key on the whole map.

    main:
    ppis_ch = datasets_ch.map { meta, f -> tuple(meta, f.ppis) }
    // Read by both SPLIT_POSITIVES (DDI mode's SELECT_EXAMPLES, where the pairs'
    // parent proteins join the same one-split-per-parent accounting as the
    // positives) and SAMPLE_NEGATIVES.
    candidate_network_ch = datasets_ch.map { meta, f -> tuple(meta, f.candidate_network) }

    // DDI mode swaps the whole data-prep front end: Pfam family accessions in
    // place of UniProt ones, domain instances in place of full chains. Both
    // subworkflows emit the same channel names, so nothing downstream branches.
    // Both take the whole files map and pick the keys they need -- which keys
    // those are is stated in each one's `take:` comment.
    if (params.ddi_mode) {
        data = DATA_PREP_DDI(datasets_ch)
    } else {
        data = DATA_PREP(datasets_ch)
    }

    if (params.split_only) {
        // --split_only: partition/node_mapping are precomputed and required
        // (validated in buildDatasetsChannel), so CLUSTERING (FETCH_DATA/
        // RUN_BLAST/MAKE_METIS/RUN_KAHIP) never needs to run at all.
        partition_ch    = datasets_ch.map { meta, f -> tuple(meta, f.partition) }
        node_mapping_ch = datasets_ch.map { meta, f -> tuple(meta, f.node_mapping) }
    } else {
        clustered = CLUSTERING(
            data.sequences, data.lengths,
            datasets_ch.map { meta, f -> tuple(meta, f.blast_results) },
            data.instances
        )
        partition_ch    = clustered.partition
        node_mapping_ch = clustered.node_mapping
    }

    split = SPLIT_POSITIVES(ppis_ch, data.sequences, partition_ch, node_mapping_ch, data.instances, candidate_network_ch)

    neg = SAMPLE_NEGATIVES(
        split.train_ppis, split.val_ppis, split.test_ppis,
        data.species, data.go_annotations,
        candidate_network_ch,
        // DDI mode only (both empty channels in PPI mode): the per-split example
        // tables, protein universes and reserves of never-in-play proteins
        // EXPAND_NEGATIVES turns family pairs into instance pairs with.
        split.ddi_files, data.instances
    )

    // --split_only stops here: SOLVE_ILP (via SPLIT_POSITIVES) + CDHIT2D +
    // REMOVE_REDUNDANT + SAMPLE_NEGATIVES_ILP have already produced and
    // published the four split files; TRAIN_BASELINE/QC add nothing this
    // mode asks for. There is then no report to emit either.
    multiqc_report_ch = channel.empty()
    if (!params.split_only) {
        baseline = TRAIN_BASELINE(
            split.train_fasta, split.val_fasta, split.test_fasta,
            neg.train, neg.val, neg.test_balanced, neg.test_realistic
        )

        qc = QC(
            neg.train, neg.val, neg.test_balanced, neg.test_realistic,
            clustered.blast_out, baseline.embeddings, data.go_annotations, data.species,
            split.train_fasta, split.val_fasta, split.test_fasta,
            split.sorted_mqc, split.nr_mqc, neg.mqc, baseline.mqc
        )
        multiqc_report_ch = qc.multiqc_report
    }

    emit:
    // Everything an embedding pipeline needs to ingest the result. The two
    // labelled channels carry the negative-set name as its own tuple field
    // rather than inside meta, because every join()/combine(by: 0) above keys on
    // the unmodified meta map -- see PLAN_domainsplit_integration.md §5. With one
    // negative-sampling method per row, which is all the samplesheet accepts
    // today, that field is just that method's name.
    instances      = data.instances     // tuple(meta, instances.tsv); (meta, []) in PPI mode
    sequences      = data.sequences     // tuple(meta, sequences.fasta)
    labelled       = neg.labelled       // tuple(meta, negset, label, csv) at family level
    labelled_inst  = neg.labelled_inst  // same, at domain-instance level; empty in PPI mode
    multiqc_report = multiqc_report_ch   // path multiqc_report.html; empty under --split_only
}

workflow {
    PPI_SPLITTING( buildDatasetsChannel() )
}
