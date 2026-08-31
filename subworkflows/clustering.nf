include { RUN_BLAST; MAKE_METIS; RUN_KAHIP } from '../processes/clustering'

// Builds the protein similarity graph (BLAST all-vs-all -> METIS graph) and
// partitions it with KaHIP, ready for SPLIT_POSITIVES. Each dataset can
// skip BLAST via the samplesheet's blast_results column.
workflow CLUSTERING {
    take:
    sequences_ch      // tuple(meta, fasta)
    lengths_ch        // tuple(meta, lengths)
    blast_results_ch  // tuple(meta, blast_results_or_[])
    instances_ch      // tuple(meta, instances_or_[]) -- [] in PPI mode

    main:
    branched = sequences_ch.join(blast_results_ch).branch { meta, fasta, blast_results ->
        precomputed: blast_results
            return tuple(meta, blast_results)
        needs_blast: true
            return tuple(meta, fasta)
    }

    blast_out = RUN_BLAST(branched.needs_blast).mix(branched.precomputed)

    // In DDI mode instances.tsv contracts the graph to one node per Pfam clan.
    // Both DATA_PREP and DATA_PREP_DDI emit `instances` for every dataset ([] in
    // PPI mode), so this join can never drop one.
    metis_out = MAKE_METIS(blast_out.join(lengths_ch).join(instances_ch))

    // The ILP splitter clusters proteins into many small KaHIP partitions
    // first, whereas the default splitter partitions straight into train/val/test.
    //
    // k is clamped to the graph's node count. kaffpa exits 0 when k exceeds it,
    // but the surplus blocks come back empty, so the ILP then sees fewer
    // non-empty clusters than it has splits to fill and fails as
    // infeasible_or_unbounded. That bites in DDI mode, where MAKE_METIS
    // contracts the graph to one node per Pfam clan and the node count can drop
    // far below ilp_kahip_k's default of 100. Computed here rather than in the
    // process body so the task's command line is untouched whenever n >= k --
    // no PPI-mode input hash moves and -resume still finds everything cached.
    kahip_inputs = metis_out.graph.map { meta, graph ->
        def k      = (meta.split_method == "ilp") ? meta.ilp_kahip_k : meta.kahip_k
        def header = graph.withReader { reader -> reader.readLine() }
        def nodes  = header ? (header.trim().split(/\s+/)[0] as int) : 0
        if (nodes > 0 && k > nodes) {
            log.warn "${meta.id}: similarity.graph has only ${nodes} nodes -- clamping k from ${k} to ${nodes}"
            k = nodes
        }
        tuple(meta, graph, k)
    }
    partition = RUN_KAHIP(kahip_inputs)

    emit:
    blast_out    = blast_out
    node_mapping = metis_out.node_mapping
    partition    = partition
}
