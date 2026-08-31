include { mqcLabel } from '../helpers/mqc_labels'

// nsuffix is "" for a row asking for a single negative set and "_<negset>"
// otherwise, computed once in subworkflows/sample_negatives.nf. It goes on the
// output filename, but never on --split-name: that string is looked up in
// utils.SPLIT_ORDER, which is what forces train -> val -> test order in the
// report, and an unknown key would sort the negative sets' rows into each other.
//
// --id takes the MultiQC *display* label instead of "${meta.id}${nsuffix}", so an
// embedding pipeline can name the report's rows in its own vocabulary. The two
// coincide unless the caller supplied meta.mqc_labels -- see helpers/mqc_labels.nf.
process SAMPLE_NEGATIVES_DEGREE {
    publishDir(path: { "${params.outdir}/${meta.id}" }, mode: 'copy', saveAs: { f -> f.endsWith('_mqc.tsv') ? null : f })
    tag "${meta.id}${nsuffix}_${label}"
    label 'error_retry'

    input:
    // Staged under a fixed name distinct from "${label}${nsuffix}.csv" -- otherwise
    // for label in {train, val} the output would overwrite the staged input
    // symlink, corrupting the upstream task's cached output on -resume.
    tuple val(meta), val(negset), val(label), val(nsuffix), path(positives, stageAs: 'positives_in.csv'), val(ratio), val(uniform)  // label: "train" | "val" | "test_balanced" | "test_realistic"

    output:
    tuple val(meta), val(negset), val(label), val(nsuffix), path("${label}${nsuffix}.csv"), emit: labelled
    tuple val(meta), path("${label}*_mqc.tsv"),                                             emit: mqc

    script:
    def uniform_flag = uniform ? '--uniform' : ''
    """
    sample_negatives.py \\
        --positives  ${positives} \\
        --output     ${label}${nsuffix}.csv \\
        --split-name ${label} \\
        --ratio      ${ratio} \\
        ${uniform_flag} \\
        --seed       ${params.seed} \\
        --id         ${mqcLabel(meta, negset)}
    """
}

process SAMPLE_NEGATIVES_ILP {
    publishDir(path: { "${params.outdir}/${meta.id}" }, mode: 'copy', saveAs: { f -> f.endsWith('_mqc.tsv') ? null : f })
    tag "${meta.id}${nsuffix}_${label}"
    label 'error_retry'
    label 'gurobi'

    input:
    // positives is staged under a fixed name distinct from any "${label}${nsuffix}.csv"
    // output -- see SAMPLE_NEGATIVES_DEGREE above for why that matters. candidate_network
    // is [] for every negative set but "ilp_candidates".
    tuple val(meta), val(negset), val(label), val(nsuffix), path(positives, stageAs: 'positives_in.csv'), val(neg_ratio), path(species), path(go_annotations), path(candidate_network)  // label: "train" | "val" | "test_balanced" | "test_realistic"
    path gurobi_license  // optional; [] if params.gurobi_license is unset

    output:
    tuple val(meta), val(negset), val(label), val(nsuffix), path("${label}${nsuffix}.csv"), emit: labelled
    tuple val(meta), path("${label}*_mqc.tsv"),                                             emit: mqc

    script:
    def cand_arg = candidate_network ? "--candidate-network ${candidate_network}" : ''
    def lic_arg  = gurobi_license    ? "--gurobi-license ${gurobi_license}"        : ''
    """
    n_positives=\$(( \$(wc -l < ${positives}) - 1 ))  # -1 for the header row
    max_candidates=\$(( 4 * n_positives * ${task.attempt} ))

    sample_negatives_ilp.py \\
        --positives          ${positives} \\
        --output             ${label}${nsuffix}.csv \\
        --split-name         ${label} \\
        --neg-ratio          ${neg_ratio} \\
        --species            ${species} \\
        --go-annotations     ${go_annotations} \\
        ${cand_arg} \\
        ${lic_arg} \\
        --lambda-degree      ${meta.neg_ilp_lambda_degree} \\
        --lambda-taxon-pair  ${meta.neg_ilp_lambda_taxon_pair} \\
        --lambda-self-loop   ${meta.neg_ilp_lambda_self_loop} \\
        --lambda-jaccard     ${meta.neg_ilp_lambda_jaccard} \\
        --solver             ${params.neg_ilp_solver} \\
        --time-limit         ${meta.neg_ilp_time_limit} \\
        --mip-gap            ${params.neg_ilp_mip_gap} \\
        --threads            ${task.cpus} \\
        --seed               ${params.seed} \\
        --diagnostics-out    ${label}_mqc.tsv \\
        --residuals-out      ${label}_residuals_mqc.tsv \\
        --max-candidates     \$max_candidates \\
        --verbose \\
        --id                 ${mqcLabel(meta, negset)}
    """
}

// DDI mode only. Both samplers above work at Pfam family level; this turns the
// family pairs they labelled into the domain-instance pairs the classifier needs
// -- positives by looking up the pairs SELECT_EXAMPLES already vetted, negatives
// by drawing carriers from that split's protein universe.
process EXPAND_NEGATIVES {
    publishDir(path: { "${params.outdir}/${meta.id}" }, mode: 'copy', saveAs: { f -> f.endsWith('_mqc.tsv') ? null : f })
    tag "${meta.id}${nsuffix}_${label}"
    label 'error_retry'

    input:
    // split is the SELECT_EXAMPLES split whose files are staged here (test for
    // both test_balanced and test_realistic); label is the sampling label. The
    // labelled CSV is staged under a fixed name because it arrives as
    // "${label}${nsuffix}.csv" -- see SAMPLE_NEGATIVES_DEGREE for why that would
    // corrupt the upstream task's cached output on -resume. Every negative set of
    // one split is handed the same SELECT_EXAMPLES bundle: it runs once per row.
    tuple val(meta), val(split), val(negset), val(label), val(nsuffix),
          path(labelled, stageAs: 'labelled_in.csv'),
          path(examples), path(candidate_examples), path(universe), path(reserve), path(fasta),
          path(instances)

    output:
    tuple val(meta), val(negset), val(label), val(nsuffix), path("${label}${nsuffix}_instances.csv"), emit: labelled
    tuple val(meta), path("${label}*_mqc.tsv"),                                                       emit: mqc

    script:
    // No --allow-shared-parents here: the only thing this task did with it was
    // decide how to slice the unclaimed pool, and SELECT_EXAMPLES now hands over
    // a per-split reserve file with that decision already made (under
    // split_method=random that file is the whole pool, in every split).
    """
    expand_negatives.py \\
        --labelled           ${labelled} \\
        --examples           ${examples} \\
        --candidate-examples ${candidate_examples} \\
        --universe           ${universe} \\
        --reserve            ${reserve} \\
        --instances          ${instances} \\
        --fasta              ${fasta} \\
        --output             ${label}${nsuffix}_instances.csv \\
        --split-name         ${label} \\
        --examples-target    ${params.ddi_examples_target} \\
        --seed               ${params.seed} \\
        --id                 ${mqcLabel(meta, negset)}
    """
}
