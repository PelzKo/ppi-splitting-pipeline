include { BIAS_ANALYSIS; COLLECT_BIAS } from '../processes/qc'
include { mqcLabel } from '../helpers/mqc_labels'

// Runs BIAS_ANALYSIS for every attribute applicable to each dataset and
// negative set, then collects each (dataset, negative set)'s results into one
// scatter plot. Shared by the full pipeline (PPI_SPLITTING feeds it
// SAMPLE_NEGATIVES' labelled splits and hands the scatter on to QC) and the
// --bias_only shortcut (fed straight from the samplesheet), so the attribute
// list and the same_species rule live in exactly one place.
workflow BIAS_DIAGNOSTICS {
    take:
    train_ppis            // tuple(meta, negset, path)
    val_ppis
    test_balanced_ppis
    test_realistic_ppis
    blast_out             // tuple(meta, path)
    embeddings
    go_annotations_ch     // tuple(meta, path_or_[]) -- [] only under --bias_only --ddi_mode
    species_ch

    main:
    // Whether to include "same_species" depends on each dataset's own
    // species.tsv, so it's computed per-dataset here rather than with a
    // single run-wide collect().
    // DDI mode drops the three GO-based attributes -- domain families carry no
    // GO annotations at all, so DATA_PREP_DDI emits a header-only table -- and
    // adds parent_degree. The other four need no change: sequence_similarity,
    // embedding_similarity and same_species act on the domain instances the rows
    // hold, while self_interactions and topology_shortcut act on the node pair,
    // which bias_analysis.py reads from the rows' own family1/family2 columns.
    // These names must match bias_analysis.py's ATTRIBUTES dict exactly -- it is
    // also the argparse `choices`, so a mismatch is a hard task failure.
    attrs_ch = species_ch.map { meta, sp ->
        def taxa = sp.splitCsv(header: true, sep: '\t').collect { it.taxon_id }.unique()
        def attrs = params.ddi_mode
            ? ["sequence_similarity", "embedding_similarity", "self_interactions",
               "topology_shortcut", "parent_degree"]
            : ["sequence_similarity", "embedding_similarity",
               "functional_relatedness_BP", "functional_relatedness_MF",
               "functional_relatedness_CC", "self_interactions",
               "topology_shortcut"]
        if (taxa.size() > 1) attrs << "same_species"
        tuple(meta, attrs)
    }.flatMap { meta, attrs -> attrs.collect { a -> tuple(meta, a) } }

    // One negative set's four labelled CSVs in one item, keyed (meta, negset) --
    // join(by: [0, 1]), because a plain 1:1 join on meta would pair one negative
    // set's train CSV with another's val CSV once a row asks for several.
    negset_splits = train_ppis.join(val_ppis, by: [0, 1])
        .join(test_balanced_ppis,  by: [0, 1])
        .join(test_realistic_ppis, by: [0, 1])
    // tuple(meta, negset, train, val, test_balanced, test_realistic)

    // blast/embeddings/go/species are one-per-dataset; combine(by: 0) broadcasts
    // each dataset's single set of files to every one of that dataset's
    // (attribute, negative set) pairs, rather than a full cross-join. The bias
    // analysis runs per negative set because that is the half of each dataset the
    // sets differ in -- the positive half is shared by construction.
    bias_inputs = negset_splits
        .combine(attrs_ch, by: 0)
        .map { meta, negset, train, val, tb, tr, attr -> tuple(meta, attr, negset, train, val, tb, tr) }
        .combine(blast_out,         by: 0)
        .combine(embeddings,        by: 0)
        .combine(go_annotations_ch, by: 0)
        .combine(species_ch,        by: 0)

    bias = BIAS_ANALYSIS(bias_inputs)

    // One scatter per (dataset, negative set): the negset-qualified id is what
    // keeps a row's two sets from being averaged into a single plot.
    bias_by_negset = bias.mqc.flatMap { meta, negset, f ->
        def files = (f instanceof List) ? f : [f]
        files.collect { ff -> tuple(mqcLabel(meta, negset), ff) }
    }
    scatter = COLLECT_BIAS(bias_by_negset.groupTuple())

    emit:
    mqc     = bias.mqc      // tuple(meta, negset, *_bias_mqc.tsv) -- published, not sent to MultiQC
    scatter = scatter.mqc   // tuple(display label, bias_scatter_mqc.html)
}
