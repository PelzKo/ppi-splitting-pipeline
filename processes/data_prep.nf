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

    input:
    tuple val(meta), path(families), path(clans)  // families: plain text, one Pfam accession per line; clans: [] to download

    output:
    tuple val(meta), path("sequences.fasta"),    emit: sequences
    tuple val(meta), path("go_annotations.tsv"), emit: go_annotations
    tuple val(meta), path("species.tsv"),        emit: species
    tuple val(meta), path("instances.tsv"),      emit: instances

    script:
    def pool_size  = (params.ddi_examples_target as int) * (params.ddi_examples_pool_factor as int)
    def clans_arg  = clans              ? "--clans ${clans}"                        : ''
    def cache_arg  = params.interpro_cache ? "--interpro-cache ${params.interpro_cache}" : ''
    def fasta_arg  = params.pfam_fasta  ? "--pfam-fasta ${params.pfam_fasta}"       : ''
    """
    fetch_domains.py \\
        --families  ${families} \\
        --pool-size ${pool_size} \\
        --seed      ${params.seed} \\
        ${clans_arg} \\
        ${cache_arg} \\
        ${fasta_arg}
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
