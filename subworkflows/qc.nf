include { DDI_ATTRITION; SIMILARITY_HEATMAP; MULTIQC } from '../processes/qc'
include { mqcRowLabel } from '../helpers/mqc_labels'

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

// Builds the train/val/test similarity heatmap and assembles one combined
// MultiQC report for the whole run from every dataset's diagnostics. The bias
// analyses and their scatter plots are computed upstream by BIAS_DIAGNOSTICS
// (subworkflows/bias_diagnostics.nf), which --bias_only also runs on its own.
workflow QC {
    take:
    bias_scatter          // tuple(display label, bias_scatter_mqc.html) -- BIAS_DIAGNOSTICS.scatter
    blast_out
    train_fasta
    val_fasta
    test_fasta
    sorted_mqc
    nr_mqc
    neg_mqc
    clf_mqc
    mqc_order             // comma-joined display labels, in report order (resolved in main.nf)

    main:
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
    // over the bias_scatter plot, which is what's kept; BIAS_ANALYSIS publishes
    // them under <id>/bias/<negset>/ instead.
    mqc_files = splitting_mqc
        .mix(flattenMqc(neg_mqc))
        .mix(flattenMqc(clf_mqc))
        .mix(bias_scatter)
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
