// Shared MultiQC display-label helpers.
//
// A functions-only module rather than a copy per consumer: main.nf, two process
// files and subworkflows/qc.nf all need the same rule, and neither existing
// layer can host it -- processes/ cannot include subworkflows/sample_negatives.nf
// (which is where negsuffix() lives) without making the include graph circular,
// since that file includes processes/negative_sampling.nf.
//
// Nothing here reads params or channels, so it is safe for every layer. Each
// function is deliberately self-contained rather than calling its neighbours:
// Nextflow binds included functions one at a time, so an intra-module call is a
// dependency on how that binding works, for no gain in three short lines.

// The MultiQC display label for one (dataset row, negative set).
//
// meta.mqc_labels is an optional Map from negative-set name to display label, so
// an embedding pipeline can make the report speak its own vocabulary
// ('ilp' -> 'minimal_leakage', 'ilp_candidates' -> 'minimal_leakage_hcni').
// An absent key, a non-Map value, or a negative set the map does not cover all
// fall back to this pipeline's own "<id><negsuffix>" -- which is what keeps a
// samplesheet-driven standalone run's labels exactly as they were.
//
// mqc_labels must be set when meta is built and never mutated afterwards: every
// join()/combine(by: 0) in the pipeline keys on the whole meta map, so a late
// mutation silently drops branches.
//
// This affects the report only. Published filenames and directories key off
// meta.id and the negative-set suffix computed in subworkflows/sample_negatives.nf,
// and --split-name keeps receiving the bare split label, which utils.SPLIT_ORDER
// is looked up by.
def mqcLabel(meta, negset) {
    def mapped = (meta.mqc_labels instanceof Map) ? meta.mqc_labels[negset] : null
    if (mapped) {
        return mapped.toString()
    }
    def multi = meta.negative_sampling_method.toString().tokenize(',').size() > 1
    // .toString() on every branch: a GString never equals an equal String, and
    // these labels are compared against params.mqc_order entries and used as
    // groupTuple() keys.
    return (multi ? "${meta.id}_${negset}" : "${meta.id}").toString()
}

// The label for an artefact belonging to the whole row rather than to one
// negative set -- the splitting stage's bars, the similarity heatmap, the DDI
// attrition waterfall. Those tasks run once per row however many negative sets
// it asked for, so they take the first set's label. A bare meta.id would not be
// one of params.mqc_order's names at all and would sort after every real label.
def mqcRowLabel(meta) {
    def names = meta.negative_sampling_method.toString().tokenize(',')*.trim().findAll { n -> n }
    def first = names ? names[0] : ""
    def mapped = (meta.mqc_labels instanceof Map) ? meta.mqc_labels[first] : null
    if (mapped) {
        return mapped.toString()
    }
    return (names.size() > 1 ? "${meta.id}_${first}" : "${meta.id}").toString()
}
