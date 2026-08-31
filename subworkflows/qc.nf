include { BIAS_ANALYSIS; COLLECT_BIAS; DDI_ATTRITION; SIMILARITY_HEATMAP; MULTIQC } from '../processes/qc'
include { mqcLabel; mqcRowLabel } from '../helpers/mqc_labels'

// Some mqc-emitting processes glob-match more than one file per task, which
// Nextflow packs into a List -- flatten to one (id, file) pair per file so
// groupTuple() below doesn't nest a List inside the grouped list.
//
// Keyed by the row's MultiQC display label rather than meta.id: for everything
// that reaches MULTIQC the key is only a grouping handle, but DDI_ATTRITION also
// passes it on as its --id, and that label is what params.mqc_order names.
def flattenMqc(ch) {
    ch.flatMap { meta, f ->
        def files = (f instanceof List) ? f : [f]
        files.collect { ff -> tuple(mqcRowLabel(meta), ff) }
    }
}

// Runs the per-attribute bias analyses, collects them into a scatter plot,
// builds the train/val/test similarity heatmap, and assembles one combined
// MultiQC report for the whole run from every dataset's diagnostics.
workflow QC {
    take:
    train_ppis            // tuple(meta, negset, path)
    val_ppis
    test_balanced_ppis
    test_realistic_ppis
    blast_out
    embeddings
    go_annotations_ch
    species_ch
    train_fasta
    val_fasta
    test_fasta
    sorted_mqc
    nr_mqc
    neg_mqc
    clf_mqc
    mqc_order             // comma-joined display labels, in report order (resolved in main.nf)

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

    // meta.id AND the display label: the first is SIMILARITY_HEATMAP's publishDir
    // component, the second only names the report section. They differ for a row
    // with several negative sets, so passing one for both would move a published
    // directory.
    heatmap_inputs = train_fasta.join(val_fasta).join(test_fasta).join(blast_out)
        .map { meta, t, v, te, b -> tuple(meta.id, mqcRowLabel(meta), t, v, te, b) }
    heatmap = SIMILARITY_HEATMAP(heatmap_inputs)

    splitting_mqc = flattenMqc(sorted_mqc).mix(flattenMqc(nr_mqc))

    // sort_ppis.py's write_mqc() is only ever reached through sort_ppis_random.py,
    // and the random path skips REMOVE_REDUNDANT entirely -- so the two writers of
    // the shared "split_bar_<label>" section id are on mutually exclusive paths and
    // this never fires today. It is here because if that ever changes MultiQC
    // silently merges the two files into one chart, keeping whichever it parses
    // last, and nothing else in the run would say so.
    splitting_mqc.groupTuple().subscribe { id, files ->
        def names = files.collect { ff -> ff.name }
        if (names.any { n -> n.startsWith('sort_ppis_bar') } && names.any { n -> n.startsWith('remove_redundant_bar') }) {
            log.warn "${id}: both sort_ppis_bar_mqc.tsv and remove_redundant_bar_mqc.tsv were produced. They declare the same MultiQC section id (split_bar_<label>), so MultiQC will merge them into one chart and keep whichever file it parses last."
        }
    }

    // One stacked bar per dataset -- not per negative set: every bar it reads back
    // is a splitting-stage bar, and the splitting stage runs once per row whatever
    // the negative sets are. Accounting for every input DDI: discarded by
    // the partitioner, removed by CD-HIT-2D, dropped because no domain-instance
    // example was left for it, or kept. It reads the counts back out of the
    // splitting stage's own MultiQC bars rather than re-deriving them, so the
    // waterfall and the per-stage charts cannot disagree -- and neither
    // splitter nor SELECT_EXAMPLES needs new instrumentation.
    ddi_attrition = params.ddi_mode ? DDI_ATTRITION(splitting_mqc.groupTuple()).mqc : channel.empty()

    // Bias tables are deliberately excluded here -- they don't add value
    // over the bias_scatter plot above, which is what's kept. bias.mqc
    // still feeds COLLECT_BIAS unconditionally, just not this final mix.
    mqc_files = splitting_mqc
        .mix(flattenMqc(neg_mqc))
        .mix(flattenMqc(clf_mqc))
        .mix(scatter.mqc)
        .mix(heatmap)
        .mix(ddi_attrition)
        .map { id, f -> f }
        .collect()

    multiqc = MULTIQC(mqc_files, mqc_order)

    emit:
    // Emitted rather than only published, so a pipeline including PPI_SPLITTING can
    // re-publish it under its own name and location without knowing this one's.
    multiqc_report = multiqc.report
}
