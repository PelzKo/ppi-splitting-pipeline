// One task per (dataset, attribute, negative set). The negative set is carried
// only so the caller can key the scatter plots by it -- nothing in the script
// reads it, and the tag names it unconditionally so a row's several sets are
// distinguishable in the log.
process BIAS_ANALYSIS {
    tag "${meta.id}_${negset}_${attribute}"
    label 'error_retry'

    input:
    tuple val(meta), val(attribute), val(negset),
          path(train_csv), path(val_csv), path(test_balanced_csv), path(test_realistic_csv),
          path(blast_tsv), path(embeddings), path(go_annotations), path(species)

    output:
    tuple val(meta), val(negset), path("*_bias_mqc.tsv"), emit: mqc, optional: true

    script:
    """
    bias_analysis.py \\
        --attribute       ${attribute} \\
        --train           ${train_csv} \\
        --val             ${val_csv} \\
        --test_balanced   ${test_balanced_csv} \\
        --test_realistic  ${test_realistic_csv} \\
        --blast           ${blast_tsv} \\
        --embeddings      ${embeddings} \\
        --go_annotations  ${go_annotations} \\
        --species         ${species} \\
        --seed            ${params.seed}
    """
}

// id is the MultiQC display label for one (dataset, negative set), computed in
// subworkflows/qc.nf -- not meta.id. See helpers/mqc_labels.nf.
process COLLECT_BIAS {
    tag "${id}_bias_scatter"

    input:
    tuple val(id), path(tsvs)

    output:
    tuple val(id), path("bias_scatter_mqc.html"), emit: mqc

    script:
    """
    collect_bias.py ${tsvs} --id ${id}
    """
}

// DDI mode only. Reduces one dataset's splitting-stage MultiQC bars into a
// single attrition waterfall; it picks the two bars it needs out of whatever it
// is handed by their column headers, so it neither depends on a filename nor
// needs the upstream scripts to emit anything extra.
// id is the whole row's MultiQC display label, computed in subworkflows/qc.nf.
process DDI_ATTRITION {
    tag "${id}_ddi_attrition"

    input:
    tuple val(id), path(tsvs)

    output:
    tuple val(id), path("ddi_attrition_mqc.tsv"), emit: mqc

    script:
    """
    ddi_attrition.py ${tsvs} --id ${id}
    """
}

// Two ids, and they must not be conflated: `id` is meta.id and is this process's
// publishDir component, `label` is the row's MultiQC display label and goes only
// into the report. They differ whenever a row asks for several negative sets (the
// label then carries the first set's suffix), so folding them into one would move
// a published directory -- which this whole feature must not do.
process SIMILARITY_HEATMAP {
    publishDir(path: { "${params.outdir}/${id}/multiqc" }, mode: 'copy')
    tag "${label}_heatmap"

    input:
    tuple val(id), val(label), path(train_fasta), path(val_fasta), path(test_fasta), path(blast_tsv)

    output:
    tuple val(id), path("similarity_heatmap_mqc.html")

    script:
    """
    plot_similarity_heatmap.py \\
        --train_fasta ${train_fasta} \\
        --val_fasta   ${val_fasta} \\
        --test_fasta  ${test_fasta} \\
        --blast       ${blast_tsv} \\
        --max_per_split ${params.heatmap_max_per_split} \\
        --seed        ${params.seed} \\
        --id          ${label}
    """
}

// One combined report for the whole run. stageAs: "?/*" stages each
// dataset's contributions into its own numbered subdirectory, since
// same-named files (e.g. classifier_metrics_mqc.tsv) would otherwise
// collide -- MultiQC scans subdirectories recursively, so this is transparent.
//
// relabel_mqc.py runs first and is the only place that knows the run's whole
// dataset order: it rewrites each file's Sample column and section id into the
// globally numbered form and writes the report_section_order config. Doing it
// here rather than in the emitting tasks keeps the run's dataset set out of the
// hash of SELECT_EXAMPLES and SAMPLE_NEGATIVES_ILP, which would otherwise
// re-solve for every dataset whenever one was added.
//
// The numbered stageAs directories are exactly why the config is needed: MultiQC
// discovers custom-content sections in directory order, which is channel arrival
// order, so without an explicit order the sections shuffle from run to run.
process MULTIQC {
    publishDir(path: { "${params.outdir}/multiqc" }, mode: 'copy')
    tag "multiqc"

    input:
    path multiqc_files, stageAs: "?/*"
    val  mqc_order  // comma-joined display labels, in report order

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data

    script:
    """
    relabel_mqc.py . \\
        --out       relabelled \\
        --mqc-order '${mqc_order}' \\
        --config    multiqc_config.yml

    multiqc relabelled \\
        -c multiqc_config.yml \\
        --title "PPI Splitting Pipeline" \\
        --filename multiqc_report.html
    """
}
