// Datasets requesting the same embedding_model share one embedding call over
// the union of their train/val/test sequences, avoiding recomputation for
// proteins in more than one dataset. stageAs auto-numbers fasta_files since
// every dataset's split fasta shares the same name (train_nr.fasta, etc.)
// and would otherwise collide; embed_sequences.py merges/dedupes them.
process EMBED_SEQUENCES {
    publishDir(path: { "${params.outdir}/_shared/embeddings" }, mode: 'copy', saveAs: { f -> "embeddings_${embedding_model}.npz" })
    tag "embed_${embedding_model}"
    label "process_gpu"

    input:
    tuple val(embedding_model), path(fasta_files, stageAs: "input_*")

    output:
    tuple val(embedding_model), path("embeddings.npz")

    script:
    def gpu_arg = workflow.profile.tokenize(',')*.trim().contains('gpu') ? '--require-gpu' : ''
    """
    embed_sequences.py \\
        --fasta ${fasta_files} \\
        --model ${embedding_model} \\
        ${gpu_arg}
    """
}

// Same rule as subworkflows/sample_negatives.nf's negsuffix(), re-derived from
// meta because this process sees one negative set at a time rather than the row's
// whole list: "" for a single-negative-set row, so its MultiQC sample names are
// unchanged, "_<negset>" otherwise so the sets do not collide in one report.
// Keep the two in step.
def negsuffix(meta, negset) {
    meta.negative_sampling_method.toString().tokenize(',').size() > 1 ? "_${negset}" : ""
}

process TRAIN_CLASSIFIER {
    tag "${meta.id}${negsuffix(meta, negset)}"
    label 'error_retry'

    input:
    tuple val(meta), val(negset), path(train_csv), path(val_csv), path(test_balanced_csv), path(test_realistic_csv), path(embeddings)

    output:
    tuple val(meta), path("classifier_metrics_*_mqc.tsv"), emit: mqc

    script:
    // test_realistic is one shared file across a row's negative sets, so its
    // metrics row repeats -- truthfully -- once per set.
    def mqc_id = "${meta.id}${negsuffix(meta, negset)}"
    """
    train_classifier.py \\
        --train          ${train_csv} \\
        --val            ${val_csv} \\
        --test_balanced  ${test_balanced_csv} \\
        --test_realistic ${test_realistic_csv} \\
        --embeddings     ${embeddings} \\
        --seed           ${params.seed} \\
        --id             ${mqc_id}
    """
}
