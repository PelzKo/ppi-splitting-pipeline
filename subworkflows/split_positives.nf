include { SORT_PPIS; SOLVE_ILP; SPLIT_RANDOM; CDHIT2D; REMOVE_REDUNDANT; SELECT_EXAMPLES } from '../processes/splitting'

// Assigns PPIs to train/val/test (KaHIP, ILP, or random shuffle, per
// meta.split_method), then runs CD-HIT-2D redundancy removal for the two
// leakage-aware methods. "random" skips CD-HIT -- see SPLIT_RANDOM.
workflow SPLIT_POSITIVES {
    take:
    ppis_ch          // tuple(meta, ppis)
    sequences_ch     // tuple(meta, fasta)
    partition_ch     // tuple(meta, partition)
    node_mapping_ch  // tuple(meta, node_mapping)
    instances_ch     // tuple(meta, instances_or_[]) -- [] in PPI mode
    candidate_network_ch  // tuple(meta, candidate_network_or_[]) -- DDI mode's SELECT_EXAMPLES only

    main:
    // instances.tsv reconciles the three id vocabularies DDI mode splits across:
    // the interaction file's Pfam families, the partition's clans and the FASTA's
    // domain instances. Every dataset carries an item (DATA_PREP emits []), so
    // these joins can never drop one.
    joined = ppis_ch.join(sequences_ch).join(partition_ch).join(node_mapping_ch).join(instances_ch)
    // tuple(meta, ppis, fasta, partition, node_mapping, instances)

    branched = joined.branch { meta, ppis, fasta, partition, node_mapping, instances ->
        ilp:    meta.split_method == "ilp"
        random: meta.split_method == "random"
        kahip:  true
    }

    gurobi_license_ch = params.gurobi_license
        ? channel.value(file(params.gurobi_license, checkIfExists: true))
        : channel.value([])

    ilp_out    = SOLVE_ILP(branched.ilp, gurobi_license_ch)
    kahip_out  = SORT_PPIS(branched.kahip)
    random_out = SPLIT_RANDOM(branched.random.map { meta, ppis, fasta, partition, node_mapping, instances -> tuple(meta, ppis, fasta, instances) })

    // ilp/kahip go through CD-HIT redundancy removal below; random doesn't
    // -- it feeds straight into the final emit channels further down.
    homology_train_ppis  = ilp_out.train_ppis.mix(kahip_out.train_ppis)
    homology_val_ppis    = ilp_out.val_ppis.mix(kahip_out.val_ppis)
    homology_test_ppis   = ilp_out.test_ppis.mix(kahip_out.test_ppis)
    homology_train_fasta = ilp_out.train_fasta.mix(kahip_out.train_fasta)
    homology_val_fasta   = ilp_out.val_fasta.mix(kahip_out.val_fasta)
    homology_test_fasta  = ilp_out.test_fasta.mix(kahip_out.test_fasta)
    splitter_mqc = ilp_out.mqc.mix(kahip_out.mqc).mix(random_out.mqc)

    // One channel, one process: each (meta, label, fasta1, fasta2) item
    // becomes its own task, so Nextflow runs both CD-HIT-2D comparisons
    // for every dataset in parallel.
    cdhit_inputs_ch = homology_train_fasta.join(homology_val_fasta).map { meta, t, v -> tuple(meta, "train_val", t, v) }
        .mix(homology_train_fasta.join(homology_test_fasta).map { meta, t, te -> tuple(meta, "train_test", t, te) })

    cdhit_out = CDHIT2D(cdhit_inputs_ch)

    cdhit_branched = cdhit_out.sim.branch {
        meta, label, f ->
            train_val:  label == "train_val"
            train_test: label == "train_test"
    }
    sim_tv = cdhit_branched.train_val.map { meta, label, f -> tuple(meta, f) }
    sim_tt = cdhit_branched.train_test.map { meta, label, f -> tuple(meta, f) }

    // ppis_ch (pre-split) is threaded in so REMOVE_REDUNDANT can compute the
    // KaHIP/ILP discard count itself for the PPI Partitioning chart.
    // instances is appended LAST here, which is the order REMOVE_REDUNDANT's
    // input block has to declare -- the tuple is positional and unchecked, so
    // an element out of place stages sim_train_test.out as the instance table.
    nr_inputs = ppis_ch.join(homology_train_ppis).join(homology_val_ppis).join(homology_test_ppis)
        .join(homology_train_fasta).join(homology_val_fasta).join(homology_test_fasta)
        .join(sim_tv).join(sim_tt).join(instances_ch)

    nr = REMOVE_REDUNDANT(nr_inputs)

    // The final per-split family CSVs and instance FASTAs, whichever path produced them.
    fam_train   = nr.train_ppis.mix(random_out.train_ppis)
    fam_val     = nr.val_ppis.mix(random_out.val_ppis)
    fam_test    = nr.test_ppis.mix(random_out.test_ppis)
    fasta_train = nr.train_fasta.mix(random_out.train_fasta)
    fasta_val   = nr.val_fasta.mix(random_out.val_fasta)
    fasta_test  = nr.test_fasta.mix(random_out.test_fasta)

    // DDI mode's last splitting step: family pairs become domain-instance pairs
    // and each parent protein is claimed by at most one split. It runs
    // over all three splits jointly -- the claim table spans them -- so it is one
    // task per dataset rather than the per-split fan-out used elsewhere.
    if (params.ddi_mode) {
        sel_inputs = fam_train.join(fam_val).join(fam_test)
            .join(fasta_train).join(fasta_val).join(fasta_test)
            .join(instances_ch).join(candidate_network_ch)

        sel = SELECT_EXAMPLES(sel_inputs, gurobi_license_ch)

        // SELECT_EXAMPLES drops DDIs that reached zero examples, so the family
        // CSVs the negative sampler sees are its filtered ones, not REMOVE_REDUNDANT's.
        out_train_ppis = sel.train_ppis
        out_val_ppis   = sel.val_ppis
        out_test_ppis  = sel.test_ppis
        // Everything EXPAND_NEGATIVES needs about one split, in a single item
        // keyed (meta, split), so SAMPLE_NEGATIVES can pick it up with one
        // combine(by: [0, 1]) -- and so the test entry is broadcast to both
        // test_balanced and test_realistic rather than joined 1:1 with one of them.
        ddi_files = sel.train_examples.join(sel.train_candidates).join(sel.train_universe).join(fasta_train)
                .map { meta, ex, cand, uni, fa -> tuple(meta, "train", ex, cand, uni, fa) }
            .mix(sel.val_examples.join(sel.val_candidates).join(sel.val_universe).join(fasta_val)
                .map { meta, ex, cand, uni, fa -> tuple(meta, "val", ex, cand, uni, fa) })
            .mix(sel.test_examples.join(sel.test_candidates).join(sel.test_universe).join(fasta_test)
                .map { meta, ex, cand, uni, fa -> tuple(meta, "test", ex, cand, uni, fa) })
        unclaimed      = sel.unclaimed
        // Mixed into sorted_mqc rather than threaded through QC as a 17th take:
        // QC only collects that channel for MULTIQC, and both are splitting-stage
        // diagnostics.
        sorted_mqc     = splitter_mqc.mix(sel.mqc)
    } else {
        out_train_ppis = fam_train
        out_val_ppis   = fam_val
        out_test_ppis  = fam_test
        ddi_files      = channel.empty()
        unclaimed      = channel.empty()
        sorted_mqc     = splitter_mqc
    }

    emit:
    train_ppis  = out_train_ppis
    val_ppis    = out_val_ppis
    test_ppis   = out_test_ppis
    train_fasta = fasta_train
    val_fasta   = fasta_val
    test_fasta  = fasta_test
    // DDI mode only; both empty in PPI mode.
    // ddi_files: tuple(meta, split, examples, candidate_examples, universe, fasta)
    // unclaimed: tuple(meta, path) -- run-wide for the dataset, not per split
    ddi_files   = ddi_files
    unclaimed   = unclaimed
    sorted_mqc  = sorted_mqc
    nr_mqc      = nr.mqc
}
