#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * TYPED PARAMETER DECLARATIONS
 *
 * Types only -- every default stays in conf/params.config, which remains "every
 * parameter of the pipeline, and nothing else".
 *
 * Nextflow's v2 (strict) parser, the default since 26.x, stopped inferring types
 * for command-line parameters. `--seed 7` arrives as the String "7" and
 * `--split_only false` as the String "false". This pipeline never calls
 * validateParameters(), so nothing rejects them -- they flow straight into the
 * code, where `if (params.split_only)` on the String "false" is *true* and a
 * numeric override silently becomes a string in a meta map and a task hash.
 * Declaring the type makes Nextflow coerce the value before any of that.
 *
 * Only the *entry* script's declarations count. A pipeline that includes
 * PPI_SPLITTING as a subworkflow does not get this block -- same rule as the
 * params in conf/params.config, which is why that file is separated out in the
 * first place -- so an embedding pipeline has to repeat these declarations in
 * its own main.nf. daisybio/domainsplit does.
 *
 * `Float`, not `Double`/`BigDecimal`/`Number`: those reject the String outright
 * instead of coercing it. `mqc_order` is deliberately absent -- resolveMqcOrder()
 * accepts both a list and a comma-separated string on purpose.
 */
params {
    split_only: Boolean
    seed: Integer
    heatmap_max_per_split: Integer

    // clustering
    cdhit_identity: Float
    cdhit_wordsize: Integer
    kahip_seed: Integer
    kahip_k: Integer
    ilp_kahip_k: Integer

    // positive split
    train_split: Float
    val_split: Float
    test_split: Float
    ilp_epsilon: Float
    ilp_max_sec: Integer

    // negative sampling
    neg_ilp_lambda_degree: Float
    neg_ilp_lambda_taxon_pair: Float
    neg_ilp_lambda_self_loop: Float
    neg_ilp_lambda_jaccard: Float
    neg_ilp_time_limit: Integer
    neg_ilp_mip_gap: Float

    // DDI mode
    ddi_mode: Boolean
    ddi_examples_target: Integer
    ddi_examples_pool_factor: Integer
    ddi_select_max_sec: Integer
    ddi_max_ilp_candidates: Integer
    ddi_lambda_diversity: Float
    ddi_shortlist_factor: Integer
    ddi_select_verbose: Boolean
    ddi_candidate_factor: Integer
}

include { samplesheetToList } from 'plugin/nf-schema'

include { DATA_PREP }        from './subworkflows/data_prep'
include { DATA_PREP_DDI }    from './subworkflows/data_prep_ddi'
include { CLUSTERING }       from './subworkflows/clustering'
include { SPLIT_POSITIVES }  from './subworkflows/split_positives'
include { SAMPLE_NEGATIVES } from './subworkflows/sample_negatives'
include { TRAIN_BASELINE }   from './subworkflows/train_baseline'
include { QC }               from './subworkflows/qc'

include { mqcLabel }        from './helpers/mqc_labels'

// samplesheetToList() represents a blank cell as [] (not null), even for
// numeric fields where 0 is a legitimate override -- so check both.
def isGiven(v) {
    !(v == null || v == [])
}

// Parse and validate one row's negative-set list. A row's
// negative_sampling_method is a comma-separated list of sampler names -- one
// entry per negative set it wants out of the single positive split, so the
// positive halves of the resulting datasets are byte-identical.
// "ilp_candidates" is the ILP sampler restricted to the row's candidate_network
// pool; "ilp" is the same sampler unrestricted.
//
// Called from PPI_SPLITTING's body rather than from buildDatasetsChannel(), so a
// pipeline that includes this workflow and builds datasets_ch in Groovy is held
// to the same rules -- the samplesheet schema only checks the string's shape, not
// the cross-column ones. (The known-name list lives inside the function on
// purpose: a top-level `def x = ...` in a Nextflow script is a local, invisible
// here, which the linter rejects outright.)
def parseNegsets(meta, candidate_network) {
    def known = ["default", "ilp", "ilp_candidates", "uniform"]
    def names = meta.negative_sampling_method.toString().tokenize(',')*.trim().findAll { n -> n }
    if (!names) {
        error("Dataset '${meta.id}': negative_sampling_method is empty -- name at least one of ${known.join(', ')}.")
    }
    def unknown = names.findAll { n -> !(n in known) }
    if (unknown) {
        error("Dataset '${meta.id}': unknown negative_sampling_method ${unknown.join(', ')} (known: ${known.join(', ')}).")
    }
    if (names.size() != new HashSet(names).size()) {
        error("Dataset '${meta.id}': negative_sampling_method '${meta.negative_sampling_method}' names a method more than once -- each negative set is produced once.")
    }
    if ("ilp_candidates" in names && !candidate_network) {
        error("Dataset '${meta.id}': negative_sampling_method lists 'ilp_candidates', which draws negatives from the candidate_network pool, but the row supplies no candidate_network.")
    }
    return names
}

// Warn when a caller's meta.mqc_labels does not cover exactly the row's negative
// sets. Not fatal: mqcLabel() falls back to "<id><negsuffix>" for the sets the
// map misses, which is a cosmetically wrong report rather than a wrong result.
//
// Lives here rather than in buildDatasetsChannel() for the same reason
// parseNegsets() does: an embedding pipeline builds meta in Groovy and is the
// only caller that ever sets mqc_labels at all.
def checkMqcLabels(meta, negsets) {
    if (!(meta.mqc_labels instanceof Map)) {
        return
    }
    def given    = meta.mqc_labels.keySet().collect { k -> k.toString() }
    def expected = negsets.collect { n -> n.toString() }
    if (given.toSorted() != expected.toSorted()) {
        log.warn "Dataset '${meta.id}': mqc_labels covers ${given.toSorted().join(', ')} but the row asks for negative sets ${expected.toSorted().join(', ')} -- the sets with no entry fall back to '<id>_<negset>' display labels."
    }
}

// Resolve the run's report order once, from every display label it will produce.
//
// The numbering is global, so no task can compute it: this is the only place that
// sees the whole run. params.mqc_order accepts a Groovy list from a config file
// and a comma-separated string from the command line. Every warning is advisory
// -- a mis-set order is a cosmetic problem, never a failed run.
def resolveMqcOrder(labels) {
    def present   = labels.collect { l -> l.toString() }.unique()
    def requested = (params.mqc_order instanceof CharSequence)
        ? params.mqc_order.toString().tokenize(',')*.trim().findAll { r -> r }
        : (params.mqc_order ?: []).collect { r -> r.toString().trim() }.findAll { r -> r }

    def ordered = []
    requested.each { r ->
        if (!(r in present)) {
            log.warn "mqc_order names '${r}', which is not a display label of this run -- ignoring it. Labels present: ${present.toSorted().join(', ')}"
        }
        else if (r in ordered) {
            log.warn "mqc_order names '${r}' more than once -- keeping the first occurrence."
        }
        else {
            ordered << r
        }
    }
    // Everything the list did not mention, appended alphabetically -- so an
    // incomplete order is still a total one, and an empty order is alphabetical.
    def leftover = (present - ordered).toSorted()
    if (requested && !ordered) {
        log.warn "mqc_order was given but matched none of this run's display labels, so the report order is alphabetical. Labels present: ${present.toSorted().join(', ')}"
    }
    return (ordered + leftover).join(',')
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
        //
        // PPI_SPLITTING's body re-checks this for every caller (an embedding
        // pipeline never runs this function). Kept here as well because it fires
        // before any channel exists and names the offending samplesheet *row*,
        // which is the better message for a standalone run -- so the two must be
        // changed together.
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
                  //
                  // meta may also carry the optional `mqc_labels` map (negative-set name ->
                  // MultiQC display label), which only affects the report -- see
                  // helpers/mqc_labels.nf. buildDatasetsChannel() never sets it, so a
                  // samplesheet run's meta, and therefore every task hash, is unchanged.

    main:
    // --split_only skips CLUSTERING/TRAIN_BASELINE/QC outright, which makes every
    // one of their precomputed-input escape hatches mandatory and pins both method
    // choices to the ILP path. buildDatasetsChannel() forces the two methods and
    // checks the five files for the samplesheet entry, but an embedding pipeline
    // builds meta itself and never reaches it -- same reason parseNegsets() lives
    // here rather than there.
    //
    // Validated, never rewritten. Forcing a method here would mean silently
    // editing a meta the caller built deliberately, and meta is the key for every
    // join()/combine(by: 0) below -- so the closure hands back the two objects it
    // was given, untouched, and a caller that asked for kahip under --split_only
    // is told rather than quietly given something else.
    checked_ch = datasets_ch.map { meta, f ->
        if (params.split_only) {
            // DDI mode has no GO annotations at all (they describe proteins, not
            // domain families) and needs the domain instance table instead.
            def required = params.ddi_mode
                ? ["sequences", "species", "domain_instances", "partition", "node_mapping"]
                : ["sequences", "go_annotations", "species", "partition", "node_mapping"]
            // The files map, not a samplesheet row: an absent optional file is []
            // here, never null, so plain falsiness is the right test.
            def missing = required.findAll { k -> !f[k] }
            if (missing) {
                error("--split_only requires every dataset to supply ${required.join(', ')} (dataset '${meta.id}' is missing ${missing.join(', ')}).")
            }
            if (meta.split_method != "ilp" || meta.negative_sampling_method != "ilp") {
                error("--split_only runs the ILP path only, but dataset '${meta.id}' asks for split_method='${meta.split_method}', negative_sampling_method='${meta.negative_sampling_method}'. A samplesheet run gets both forced in buildDatasetsChannel(); a caller that builds meta itself must set them.")
            }
        }
        tuple(meta, f)
    }

    ppis_ch = checked_ch.map { meta, f -> tuple(meta, f.ppis) }

    // Validate each row's negative-set list once, and gate its candidate network
    // on it in the same pass: the network is read only when the row actually asked
    // for the negative set that uses it, so "ilp" behaves identically whether or
    // not a stray network is supplied. Gating here rather than in each consumer is
    // what keeps SELECT_EXAMPLES and SAMPLE_NEGATIVES_ILP from disagreeing about
    // whether the network is in play.
    prepared_ch = checked_ch.map { meta, f ->
        def negsets = parseNegsets(meta, f.candidate_network)
        checkMqcLabels(meta, negsets)
        if (f.candidate_network && !("ilp_candidates" in negsets)) {
            log.warn "${meta.id}: candidate_network supplied but negative_sampling_method does not list 'ilp_candidates' -- ignoring the network everywhere, including SELECT_EXAMPLES"
        }
        tuple(meta, negsets, ("ilp_candidates" in negsets) ? f.candidate_network : [])
    }
    // tuple(meta, [negset, ...]) -- the negative sets this row wants, in list order.
    negsets_ch = prepared_ch.map { meta, negsets, cand -> tuple(meta, negsets) }
    // Read by both SPLIT_POSITIVES (DDI mode's SELECT_EXAMPLES, where the pairs'
    // parent proteins join the same one-split-per-parent accounting as the
    // positives) and SAMPLE_NEGATIVES. [] unless the row asked for ilp_candidates.
    candidate_network_ch = prepared_ch.map { meta, negsets, cand -> tuple(meta, cand) }

    // Every MultiQC display label this run will produce -- one per (dataset row,
    // negative set) -- collected into the single comma-joined order string that
    // MULTIQC's relabel step numbers the report by. One item, so it pairs 1:1 with
    // MULTIQC's own collect()ed file list.
    //
    // Only MULTIQC reads it. Handing it to the tasks that emit the *_mqc.tsv files
    // instead would put the run's dataset set into the hash of SELECT_EXAMPLES and
    // SAMPLE_NEGATIVES_ILP, so adding one dataset would re-solve every other one on
    // -resume.
    mqc_order_ch = prepared_ch
        .flatMap { meta, negsets, cand -> negsets.collect { ns -> mqcLabel(meta, ns) } }
        .collect()
        .map { labels -> resolveMqcOrder(labels) }

    // DDI mode swaps the whole data-prep front end: Pfam family accessions in
    // place of UniProt ones, domain instances in place of full chains. Both
    // subworkflows emit the same channel names, so nothing downstream branches.
    // Both take the whole files map and pick the keys they need -- which keys
    // those are is stated in each one's `take:` comment.
    if (params.ddi_mode) {
        data = DATA_PREP_DDI(checked_ch)
    } else {
        data = DATA_PREP(checked_ch)
    }

    if (params.split_only) {
        // --split_only: partition/node_mapping are precomputed and required
        // (validated above, for every caller), so CLUSTERING (FETCH_DATA/
        // RUN_BLAST/MAKE_METIS/RUN_KAHIP) never needs to run at all.
        partition_ch    = checked_ch.map { meta, f -> tuple(meta, f.partition) }
        node_mapping_ch = checked_ch.map { meta, f -> tuple(meta, f.node_mapping) }
    } else {
        clustered = CLUSTERING(
            data.sequences, data.lengths,
            checked_ch.map { meta, f -> tuple(meta, f.blast_results) },
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
        split.ddi_files, data.instances,
        negsets_ch
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
            split.sorted_mqc, split.nr_mqc, neg.mqc, baseline.mqc,
            mqc_order_ch
        )
        multiqc_report_ch = qc.multiqc_report
    }

    emit:
    // Everything an embedding pipeline needs to ingest the result. The two
    // labelled channels carry the negative-set name as its own tuple field
    // rather than inside meta, because every join()/combine(by: 0) above keys on
    // the unmodified meta map -- see PLAN_domainsplit_integration.md §5. That
    // field is one entry of the row's comma-separated negative_sampling_method:
    // a row listing N of them emits N items per split, all over the same positive
    // rows, so a consumer that cares which set a CSV belongs to must join on
    // `by: [0, 1]` rather than on meta alone.
    instances        = data.instances     // tuple(meta, instances.tsv); (meta, []) in PPI mode
    sequences        = data.sequences     // tuple(meta, sequences.fasta)
    labelled         = neg.labelled       // tuple(meta, negset, label, csv) at family level
    labelled_inst    = neg.labelled_inst  // same, at domain-instance level; empty in PPI mode
    multiqc_report   = multiqc_report_ch  // path multiqc_report.html; empty under --split_only
    // tuple([id: "_shared"], dropped_families.tsv) -- the families Pfam had nothing
    // usable for, with a reason column. At most one item, and its meta is the
    // synthetic shared one, so do not join() it against any per-dataset channel.
    // Empty in PPI mode and whenever no dataset needed a fetch.
    dropped_families = data.dropped_families
}

workflow {
    PPI_SPLITTING( buildDatasetsChannel() )
}
