process SORT_PPIS {
    tag "${meta.id}"

    input:
    // instances is DDI mode's family/clan/instance table, [] in PPI mode. It is
    // declared LAST because the channel appends it with .join() -- Nextflow
    // matches tuple elements positionally and never checks, so a slot out of
    // order silently stages the wrong file.
    tuple val(meta), path(ppis), path(fasta), path(partition), path(node_mapping), path(instances)

    output:
    tuple val(meta), path("train.csv"),   emit: train_ppis
    tuple val(meta), path("val.csv"),     emit: val_ppis
    tuple val(meta), path("test.csv"),    emit: test_ppis
    tuple val(meta), path("train.fasta"), emit: train_fasta
    tuple val(meta), path("val.fasta"),   emit: val_fasta
    tuple val(meta), path("test.fasta"),  emit: test_fasta
    tuple val(meta), path("*_mqc.tsv"),   emit: mqc, optional: true

    script:
    def inst_arg = instances ? "--instances ${instances}" : ''
    """
    sort_ppis.py \\
        --ppis         ${ppis} \\
        --partition    ${partition} \\
        --fasta        ${fasta} \\
        --node_mapping ${node_mapping} \\
        ${inst_arg}
    """
}

// Deliberately naive baseline: shuffles PPIs randomly instead of using a
// KaHIP partition, so the same protein can land in more than one split --
// see bin/bias_analysis.py's "topology_shortcut" attribute. No CD-HIT runs
// downstream, since it would just strip the shared proteins back out.
process SPLIT_RANDOM {
    tag "${meta.id}"

    input:
    tuple val(meta), path(ppis), path(fasta), path(instances)  // instances last; [] in PPI mode

    output:
    tuple val(meta), path("train.csv"),   emit: train_ppis
    tuple val(meta), path("val.csv"),     emit: val_ppis
    tuple val(meta), path("test.csv"),    emit: test_ppis
    tuple val(meta), path("train.fasta"), emit: train_fasta
    tuple val(meta), path("val.fasta"),   emit: val_fasta
    tuple val(meta), path("test.fasta"),  emit: test_fasta
    tuple val(meta), path("*_mqc.tsv"),   emit: mqc

    script:
    def inst_arg = instances ? "--instances ${instances}" : ''
    """
    sort_ppis_random.py \\
        --ppis        ${ppis} \\
        --fasta       ${fasta} \\
        --train-split ${meta.train_split} \\
        --val-split   ${meta.val_split} \\
        --test-split  ${meta.test_split} \\
        --seed        ${params.seed} \\
        --id          ${meta.id} \\
        ${inst_arg}
    """
}

process CDHIT2D {
    tag "${meta.id}_${label}"

    input:
    tuple val(meta), val(label), path(db1_fasta), path(db2_fasta)  // label: "train_val" | "train_test"

    output:
    tuple val(meta), val(label), path("cdhit.out"), emit: sim

    script:
    // cd-hit-2d errors out on an empty -i2, which happens whenever db2's split has a
    // 0 fraction. An empty cdhit.out is the correct answer anyway -- REMOVE_REDUNDANT
    // reads it as "nothing to keep" and the split stays empty -- so guard here rather
    // than filtering the channel, which would drop the dataset from the 9-way join.
    """
    if [ -s ${db2_fasta} ]; then
        cd-hit-2d \\
            -i  ${db1_fasta} \\
            -i2 ${db2_fasta} \\
            -o  cdhit.out \\
            -c  ${meta.cdhit_identity} \\
            -n  ${meta.cdhit_wordsize} \\
            -T  ${task.cpus} \\
            -M  4000
    else
        touch cdhit.out
    fi
    """
}

process SOLVE_ILP {
    tag "${meta.id}"
    label 'error_retry'
    label 'gurobi'

    input:
    tuple val(meta), path(ppis), path(fasta), path(partition), path(node_mapping), path(instances)  // instances last
    path gurobi_license

    output:
    tuple val(meta), path("train.csv"),   emit: train_ppis
    tuple val(meta), path("val.csv"),     emit: val_ppis
    tuple val(meta), path("test.csv"),    emit: test_ppis
    tuple val(meta), path("train.fasta"), emit: train_fasta
    tuple val(meta), path("val.fasta"),   emit: val_fasta
    tuple val(meta), path("test.fasta"),  emit: test_fasta
    tuple val(meta), path("*_mqc.tsv"),   emit: mqc, optional: true

    script:
    def license_export = gurobi_license ? "export GRB_LICENSE_FILE=\$PWD/${gurobi_license}" : ""
    def inst_arg       = instances ? "--instances ${instances}" : ''
    """
    ${license_export}
    solve_ilp.py \\
        --ppis          ${ppis} \\
        --fasta         ${fasta} \\
        --partition     ${partition} \\
        --node_mapping  ${node_mapping} \\
        ${inst_arg} \\
        --train-split ${meta.train_split} \\
        --val-split   ${meta.val_split} \\
        --test-split  ${meta.test_split} \\
        --epsilon     ${meta.ilp_epsilon} \\
        --max-sec     ${meta.ilp_max_sec} \\
        --seed        ${params.seed} \\
        ${params.ilp_solver ? "--solver ${params.ilp_solver}" : ""}
    """
}

process REMOVE_REDUNDANT {
    tag "${meta.id}"

    input:
    tuple val(meta),
          path(orig_ppis),
          path(train_ppis), path(val_ppis), path(test_ppis),
          path(train_fasta), path(val_fasta), path(test_fasta),
          path(sim_train_val,  stageAs: 'sim_train_val.out'),
          path(sim_train_test, stageAs: 'sim_train_test.out'),
          path(instances)  // DDI mode's instances.tsv, [] in PPI mode -- LAST, matching the .join() order

    output:
    tuple val(meta), path("train_nr.csv"),   emit: train_ppis
    tuple val(meta), path("val_nr.csv"),     emit: val_ppis
    tuple val(meta), path("test_nr.csv"),    emit: test_ppis
    tuple val(meta), path("train_nr.fasta"), emit: train_fasta
    tuple val(meta), path("val_nr.fasta"),   emit: val_fasta
    tuple val(meta), path("test_nr.fasta"),  emit: test_fasta
    tuple val(meta), path("*_mqc.tsv"),      emit: mqc

    script:
    def inst_arg = instances ? "--instances ${instances}" : ''
    """
    remove_redundant.py \\
        --ppis           ${orig_ppis} \\
        --train_ppis     ${train_ppis} \\
        --val_ppis       ${val_ppis} \\
        --test_ppis      ${test_ppis} \\
        --train_fasta    ${train_fasta} \\
        --val_fasta      ${val_fasta} \\
        --test_fasta     ${test_fasta} \\
        --sim_train_val  ${sim_train_val} \\
        --sim_train_test ${sim_train_test} \\
        --id             ${meta.id} \\
        ${inst_arg}
    """
}
