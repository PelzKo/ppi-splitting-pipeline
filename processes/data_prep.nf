process FETCH_DATA {
    publishDir(path: { "${params.outdir}/${meta.id}/data" }, mode: 'copy')
    tag "${meta.id}"

    input:
    tuple val(meta), path(proteins)  // plain text file, one protein ID per line

    output:
    tuple val(meta), path("sequences.fasta"),    emit: sequences
    tuple val(meta), path("go_annotations.tsv"), emit: go_annotations
    tuple val(meta), path("species.tsv"),        emit: species

    script:
    """
    fetch_data.py ${proteins} sequences.fasta go_annotations.tsv species.tsv
    """
}

// DDI mode's FETCH_DATA: turns a list of Pfam family accessions into sampled
// domain instances plus the same sequences/species/GO artefacts DATA_PREP emits
// in PPI mode, from bulk Pfam/UniProt downloads and one streaming pass.
process FETCH_DOMAIN_META {
    publishDir(path: { "${params.outdir}/${meta.id}/data" }, mode: 'copy')
    tag "${meta.id}"
    // One ~6.3 GB stream and no checkpoint. _open_url's own retry loop covers
    // the open only: a mid-stream truncation surfaces as an EOFError out of the
    // gzip reader and kills the task, so without this label one EBI hiccup takes
    // out a six-hour job and everything downstream of it.
    label 'error_retry'

    input:
    tuple val(meta), path(families)  // plain text, one Pfam family accession per line
    path clans                       // Pfam-A.clans.tsv(.gz), or [] to download it
    path pfam_fasta                  // local Pfam-A.fasta.gz, or [] to stream it from EBI
    val cache_dir                    // absolute cache directory, or '' for no cache

    output:
    tuple val(meta), path("sequences.fasta"),        emit: sequences
    tuple val(meta), path("go_annotations.tsv"),     emit: go_annotations
    tuple val(meta), path("species.tsv"),            emit: species
    tuple val(meta), path("instances.tsv"),          emit: instances
    // Always written, header-only when nothing dropped, so this stays a required
    // output. Published rather than turned into a MultiQC section: the fetch runs
    // once per run under the synthetic "_shared" meta, so a section could only be
    // run-wide, and QC does not run at all under --split_only.
    tuple val(meta), path("dropped_families.tsv"),   emit: dropped_families

    script:
    def pool_size   = (params.ddi_examples_target as int) * (params.ddi_examples_pool_factor as int)
    def clans_arg   = clans      ? "--clans ${clans}"           : ''
    def fasta_arg   = pfam_fasta ? "--pfam-fasta ${pfam_fasta}" : ''
    // The one input that cannot be staged, because the task writes to it: a
    // relative --interpro-cache would resolve inside the work directory and be
    // thrown away with it. DATA_PREP_DDI absolutises the param for that reason.
    def cache_arg   = cache_dir  ? "--interpro-cache ${cache_dir}" : ''
    def release_arg = params.pfam_release ? "--pfam-release ${params.pfam_release}" : ''
    """
    fetch_domains.py \\
        --families  ${families} \\
        --pool-size ${pool_size} \\
        --seed      ${params.seed} \\
        --tier      ${params.instance_tier} \\
        ${clans_arg} \\
        ${cache_arg} \\
        ${fasta_arg} \\
        ${release_arg}
    """
}

process GET_LENGTHS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("lengths.tsv")

    script:
    """
    { printf 'protein_id\\tlength\\n'; \
      awk '/^>/{if(acc) print acc"\\t"len; acc=substr(\$1,2); len=0; next} {len+=length(\$0)} END{if(acc) print acc"\\t"len}' ${fasta} \
          | sort; \
    } > lengths.tsv
    """
}

// Splits the shared FETCH_DATA/GET_LENGTHS batch back out per dataset, so
// downstream steps (e.g. BLAST) see only this dataset's own proteins.
process SUBSET_FETCHED_DATA {
    publishDir(path: { "${params.outdir}/${meta.id}/data" }, mode: 'copy', saveAs: { f -> f == 'lengths.tsv' ? null : f })
    tag "${meta.id}"

    input:
    tuple val(meta), path(ppis)
    path shared_sequences,      stageAs: 'shared_sequences.fasta'
    path shared_go_annotations, stageAs: 'shared_go_annotations.tsv'
    path shared_species,        stageAs: 'shared_species.tsv'
    path shared_lengths,        stageAs: 'shared_lengths.tsv'

    output:
    tuple val(meta), path("sequences.fasta"),    emit: sequences
    tuple val(meta), path("go_annotations.tsv"), emit: go_annotations
    tuple val(meta), path("species.tsv"),        emit: species
    tuple val(meta), path("lengths.tsv"),        emit: lengths

    script:
    """
    subset_fetched_data.py \\
        --ppis ${ppis} \\
        --sequences ${shared_sequences} \\
        --go_annotations ${shared_go_annotations} \\
        --species ${shared_species} \\
        --lengths ${shared_lengths} \\
        --out_sequences sequences.fasta \\
        --out_go_annotations go_annotations.tsv \\
        --out_species species.tsv \\
        --out_lengths lengths.tsv
    """
}

// DDI mode's SUBSET_FETCHED_DATA. Same script, plus --instances: the DDI CSV's
// id columns hold Pfam families while sequences/lengths are instance-keyed, so
// the keep-set is resolved through the shared instances.tsv rather than taken
// from the CSV directly (which would silently subset to nothing).
process SUBSET_DOMAIN_DATA {
    publishDir(path: { "${params.outdir}/${meta.id}/data" }, mode: 'copy', saveAs: { f -> f == 'lengths.tsv' ? null : f })
    tag "${meta.id}"

    input:
    tuple val(meta), path(ddis)
    path shared_sequences,      stageAs: 'shared_sequences.fasta'
    path shared_go_annotations, stageAs: 'shared_go_annotations.tsv'
    path shared_species,        stageAs: 'shared_species.tsv'
    path shared_lengths,        stageAs: 'shared_lengths.tsv'
    path shared_instances,      stageAs: 'shared_instances.tsv'

    output:
    tuple val(meta), path("sequences.fasta"),    emit: sequences
    tuple val(meta), path("go_annotations.tsv"), emit: go_annotations
    tuple val(meta), path("species.tsv"),        emit: species
    tuple val(meta), path("lengths.tsv"),        emit: lengths
    tuple val(meta), path("instances.tsv"),      emit: instances

    script:
    """
    subset_fetched_data.py \\
        --ppis ${ddis} \\
        --sequences ${shared_sequences} \\
        --go_annotations ${shared_go_annotations} \\
        --species ${shared_species} \\
        --lengths ${shared_lengths} \\
        --instances ${shared_instances} \\
        --out_sequences sequences.fasta \\
        --out_go_annotations go_annotations.tsv \\
        --out_species species.tsv \\
        --out_lengths lengths.tsv \\
        --out_instances instances.tsv
    """
}

// DDI mode drops the GO bias attributes (they describe proteins, not domain
// families), but --go_annotations stays a required flag downstream. A dataset
// that supplies its own domain data therefore still needs the file to exist;
// FETCH_DOMAIN_META writes the same header-only table for the fetched path.
process EMPTY_GO_ANNOTATIONS {
    publishDir(path: { "${params.outdir}/${meta.id}/data" }, mode: 'copy')
    tag "${meta.id}"

    input:
    val meta

    output:
    tuple val(meta), path("go_annotations.tsv")

    script:
    """
    printf 'protein_id\\tgo_bp\\tgo_mf\\tgo_cc\\n' > go_annotations.tsv
    """
}
