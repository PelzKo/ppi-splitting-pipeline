include { BIAS_ANALYSIS; COLLECT_BIAS; DDI_ATTRITION; SIMILARITY_HEATMAP; MULTIQC } from '../processes/qc'

// Some mqc-emitting processes glob-match more than one file per task, which
// Nextflow packs into a List -- flatten to one (id, file) pair per file so
// groupTuple() below doesn't nest a List inside the grouped list.
def flattenMqc(ch) {
    ch.flatMap { meta, f ->
        def files = (f instanceof List) ? f : [f]
        files.collect { ff -> tuple(meta.id, ff) }
    }
}

// Runs the per-attribute bias analyses, collects them into a scatter plot,
// builds the train/val/test similarity heatmap, and assembles one combined
// MultiQC report for the whole run from every dataset's diagnostics.
workflow QC {
    take:
    train_ppis            // tuple(meta, path)
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

    // train/val/test/blast/embeddings/go/species are one-per-dataset;
    // combine(by: 0) broadcasts each dataset's single set of files to
    // every one of that dataset's attributes, rather than a full cross-join.
    bias_inputs = attrs_ch
        .combine(train_ppis,          by: 0)
        .combine(val_ppis,            by: 0)
        .combine(test_balanced_ppis,  by: 0)
        .combine(test_realistic_ppis, by: 0)
        .combine(blast_out,           by: 0)
        .combine(embeddings,          by: 0)
        .combine(go_annotations_ch,   by: 0)
        .combine(species_ch,          by: 0)

    bias = BIAS_ANALYSIS(bias_inputs)

    scatter = COLLECT_BIAS(flattenMqc(bias.mqc).groupTuple())

    heatmap_inputs = train_fasta.join(val_fasta).join(test_fasta).join(blast_out)
        .map { meta, t, v, te, b -> tuple(meta.id, t, v, te, b) }
    heatmap = SIMILARITY_HEATMAP(heatmap_inputs)

    splitting_mqc = flattenMqc(sorted_mqc).mix(flattenMqc(nr_mqc))

    // One stacked bar per dataset accounting for every input DDI: discarded by
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

    MULTIQC(mqc_files)
}
