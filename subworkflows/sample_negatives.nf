include { SAMPLE_NEGATIVES_DEGREE; SAMPLE_NEGATIVES_ILP; EXPAND_NEGATIVES } from '../processes/negative_sampling'

// "" when the row asked for a single negative set -- so every output filename is
// byte-identical to what the pipeline produced before the fan-out existed -- and
// "_<negset>" when it asked for several.
//
// Output filenames only. The same rule, keyed off meta rather than the list, is
// the fallback inside helpers/mqc_labels.nf's mqcLabel(), which is what names the
// MultiQC rows -- so a caller can rename those without touching any path.
def negsuffix(negsets, negset) {
    negsets.size() > 1 ? "_${negset}" : ""
}

// Re-emit the one test_realistic item under every negative-set name the row asked
// for. That split cannot depend on the negative set: the ILP branch excludes it by
// design and "uniform"/"default" both hand it to the uniform sampler, so every
// negative set routes it to the same uniform 10x sampling over the same shared
// positives with the same seed. It is therefore sampled -- and in DDI mode
// expanded -- once, and published without a suffix.
def broadcastRealistic(ch, negsets_ch) {
    def branched = ch.branch { meta, negset, label, suffix, f ->
        shared: label == "test_realistic"
        other:  true
    }
    return branched.other.map { meta, negset, label, suffix, f -> tuple(meta, negset, label, f) }
        .mix(branched.shared.combine(negsets_ch, by: 0).flatMap { meta, negset, label, suffix, f, negsets ->
            negsets.collect { ns -> tuple(meta, ns, label, f) }
        })
}

// Samples negative PPIs for each split, using the degree-preserving
// sampler, a fully-uniform variant, or the bias-aware ILP sampler --
// selected per negative set, which is what meta.negative_sampling_method's
// comma-separated list names.
//
// In DDI mode the samplers still run at Pfam family level -- families are
// split-exclusive and species.tsv carries family rows, so every bias term works
// unchanged -- and EXPAND_NEGATIVES then turns each labelled family pair into
// domain-instance pairs, which is the vocabulary the embeddings are keyed by.
workflow SAMPLE_NEGATIVES {
    take:
    train_ppis            // tuple(meta, path)
    val_ppis               // tuple(meta, path)
    test_ppis               // tuple(meta, path)
    species_ch              // tuple(meta, path)
    go_annotations_ch        // tuple(meta, path)
    candidate_network_ch      // tuple(meta, candidate_network_or_[]); [] unless the row asked for ilp_candidates
    ddi_files_ch               // tuple(meta, split, examples, cand_examples, universe, reserve, fasta); empty in PPI mode
    instances_ch                // tuple(meta, instances_or_[])
    negsets_ch                   // tuple(meta, [negset, ...]) -- validated in main.nf's parseNegsets()

    main:
    // One channel, one process: each (meta, negset, label) item becomes its own
    // task, so Nextflow runs all four splits of every negative set of every
    // dataset in parallel. The negative set is its own tuple field and never a
    // meta key -- adding it to meta would rekey every join()/combine(by: 0)
    // against the fasta/species/instances channels and silently empty them.
    fan_splits = train_ppis.map { meta, f -> tuple(meta, "train", f, 1.0) }
        .mix(val_ppis.map  { meta, f -> tuple(meta, "val", f, 1.0) })
        .mix(test_ppis.map { meta, f -> tuple(meta, "test_balanced", f, 1.0) })
        .combine(negsets_ch, by: 0)
        .flatMap { meta, label, f, ratio, negsets ->
            negsets.collect { ns -> tuple(meta, ns, label, negsuffix(negsets, ns), f, ratio) }
        }

    // ... except test_realistic, which is sampled once for every negative set at
    // one shared, unsuffixed filename -- see broadcastRealistic() above for why
    // its content cannot differ between them. It still needs *a* negative set to
    // pick a sampler branch below; the first one gives the same call as any other.
    realistic_split = test_ppis.combine(negsets_ch, by: 0)
        .map { meta, f, negsets ->
            if (negsets.size() > 1) {
                log.info "${meta.id}: test_realistic is identical for every negative set (${negsets.join(', ')}) -- sampling it once and publishing it unsuffixed as test_realistic.csv"
            }
            tuple(meta, negsets[0], "test_realistic", "", f, 10.0)
        }

    splits_ch = fan_splits.mix(realistic_split)

    branched = splits_ch.branch { meta, negset, label, suffix, f, ratio ->
        // Both ILP negative sets use the same sampler; only the candidate-network
        // restriction differs, applied below. test_realistic always falls through
        // to "default" (uniform), even under "ilp" -- it simulates an uncontrolled
        // random screen, so bias-matching it would defeat the point, which is also
        // why that split is negative-set invariant.
        ilp:     negset in ["ilp", "ilp_candidates"] && label != "test_realistic"
        // Deliberately naive baseline: uniform negatives for every split,
        // paired with split_method=random, to showcase the degree shortcut
        // instead of avoiding it -- see bin/bias_analysis.py's
        // topology_shortcut attribute.
        uniform: negset == "uniform"
        default: true
    }

    neg_gurobi_license_ch = params.gurobi_license
        ? channel.value(file(params.gurobi_license, checkIfExists: true))
        : channel.value([])

    // species/go_annotations/candidate_network are one-per-dataset;
    // combine(by: 0) broadcasts each dataset's single file to every one
    // of that dataset's (up to 4 per negative set) splits, rather than a
    // full cross-join.
    ilp_inputs = branched.ilp
        .combine(species_ch, by: 0)
        .combine(go_annotations_ch, by: 0)
        .combine(candidate_network_ch, by: 0)
        // The candidate network belongs to the "ilp_candidates" negative set alone.
        // A row asking for both gets one unrestricted and one restricted sampler
        // over the same positives, which is the whole point of the fan-out.
        .map { meta, negset, label, suffix, f, ratio, species, go, cand ->
            tuple(meta, negset, label, suffix, f, ratio, species, go, negset == "ilp_candidates" ? cand : [])
        }
    ilp_out = SAMPLE_NEGATIVES_ILP(ilp_inputs, neg_gurobi_license_ch)

    uniform_inputs = branched.uniform.map { meta, negset, label, suffix, f, ratio ->
        tuple(meta, negset, label, suffix, f, ratio, true)
    }
    // label == "test_realistic" here also catches the ILP negative sets'
    // test_realistic item, since the branch above routes it here instead of
    // to SAMPLE_NEGATIVES_ILP.
    default_inputs = branched.default.map { meta, negset, label, suffix, f, ratio ->
        tuple(meta, negset, label, suffix, f, ratio, label == "test_realistic")
    }
    degree_out = SAMPLE_NEGATIVES_DEGREE(uniform_inputs.mix(default_inputs))

    fam_labelled = ilp_out.labelled.mix(degree_out.labelled)  // tuple(meta, negset, label, suffix, csv)
    fam_mqc      = ilp_out.mqc.mix(degree_out.mqc)

    if (params.ddi_mode) {
        // combine(by: [0, 1]) rather than join: test_balanced and test_realistic
        // both key onto the "test" bundle, and join is 1:1 so it would silently
        // drop one of them -- as would the several negative sets of one split.
        exp_inputs = fam_labelled
            .map { meta, negset, label, suffix, f ->
                tuple(meta, label.startsWith("test") ? "test" : label, negset, label, suffix, f)
            }
            .combine(ddi_files_ch, by: [0, 1])
            .combine(instances_ch, by: 0)
        // tuple(meta, split, negset, label, suffix, labelled, examples, cand_examples, universe, reserve, fasta, instances)

        exp_out       = EXPAND_NEGATIVES(exp_inputs)
        neg_labelled  = broadcastRealistic(exp_out.labelled, negsets_ch)
        neg_mqc       = fam_mqc.mix(exp_out.mqc)
        inst_labelled = neg_labelled
        out_fam       = broadcastRealistic(fam_labelled, negsets_ch)
    } else {
        neg_labelled  = broadcastRealistic(fam_labelled, negsets_ch)
        neg_mqc       = fam_mqc
        inst_labelled = channel.empty()
        out_fam       = neg_labelled
    }

    neg_branched = neg_labelled.branch {
        meta, negset, label, f ->
            train:          label == "train"
            val:            label == "val"
            test_balanced:  label == "test_balanced"
            test_realistic: label == "test_realistic"
    }

    // The four per-split emits below are what TRAIN_BASELINE and QC consume, keyed
    // by (meta, negset) so their joins stay 1:1 across the fan-out; the two
    // `labelled*` ones exist for a pipeline that includes PPI_SPLITTING and wants
    // every split in one channel, keyed by negative set.

    emit:
    train          = neg_branched.train.map          { meta, negset, label, f -> tuple(meta, negset, f) }
    val            = neg_branched.val.map            { meta, negset, label, f -> tuple(meta, negset, f) }
    test_balanced  = neg_branched.test_balanced.map  { meta, negset, label, f -> tuple(meta, negset, f) }
    test_realistic = neg_branched.test_realistic.map { meta, negset, label, f -> tuple(meta, negset, f) }
    mqc            = neg_mqc
    labelled       = out_fam        // tuple(meta, negset, label, csv) at family level
    labelled_inst  = inst_labelled  // same, at domain-instance level; empty in PPI mode
}
