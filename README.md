![Logo](logo.png)

Automated leakage-aware splitting of a protein–protein interaction (PPI) dataset into train, validation, and test sets, with redundancy removal, negative sampling, embedding-based classification, and bias analysis.

With `--ddi_mode` the same pipeline splits **domain–domain interactions** (Pfam family pairs) instead of PPIs — see [DDI mode](#ddi-mode) below.

Have a look at the [Wiki](https://github.com/bionetslab/ppi-splitting-pipeline/wiki) for more information.


![Pipeline overview](metro_map.svg)

## Quick Start

### Input PPI File

To custom-split your PPI dataset, you need to provide it to the pipeline as a CSV file with at least two columns (`protein1`, `protein2`) containing UniProt accession IDs. Additional columns (e.g., STRING evidence scores) are preserved throughout the pipeline.

```
protein1,protein2
P45985,Q14315
Q86TC9,P35609
O14836-2,P12345
...
```

In DDI mode the file has exactly the same shape, but the two columns hold **Pfam
family accessions** and one row is one domain–domain interaction:

```
protein1,protein2
PF00069,PF00017
PF00069,PF00028
...
```

The column names don't change, so every downstream step reads the file the same way.

### Samplesheet preparation

You provide all parameters for the pipeline via a samplesheet CSV where one row corresponds to one run. E.g., :

| id        | ppis             | split_method | negative_sampling_method | neg_ilp_solver | gurobi_license     |
|-----------|------------------|--------------|--------------------------|----------------|--------------------|
| fast-run  | data/my_ppis.csv | kahip        | default                  |                |                    |
| split-ilp | data/my_ppis.csv | ilp          | default                  |                |                    |
| all-ilp   | data/my_ppis.csv | ilp          | ilp                      | gurobi         | path/to/gurobi.lic |

### Run the pipeline

If you have a GPU, `-profile gpu` will submit the embedding step to a GPU, as specified by your nextflow config.
It also makes `EMBED_SEQUENCES` **fail** rather than fall back to CPU when torch finds no usable CUDA
device — a CPU run is roughly 60× slower and would only hit the scheduler's walltime hours later. The
usual cause is a torch wheel built for a newer CUDA than the node's driver (`environment.yml` pins a
cu12x build for this reason). Without `-profile gpu`, CPU is the expected device and the step just warns.

```
nextflow run main.nf --samplesheet samplesheet.csv --outdir results -profile gpu -c my_config.config
```

In `my_config.config`, you have to specify where the GPU is located. SLURM example:

```
profiles {
    ...
    gpu {
        process {
            withLabel:process_gpu {
                queue = 'shared-gpu'
                clusterOptions = '--qos=limitgpus --gpus=a40:1 --exclude small-gpu'
            }
        }
    }
}
```

### View the report

The MultiQC report can be found at `results/multiqc/multiqc_report.html`, which you can view in a browser.

---

## Workflow

![Pipeline overview](metro_map.svg)


### Step descriptions

**FETCH_DATA** — Queries UniProt for the union of unique proteins across every samplesheet dataset that needs a fetch (extracted directly from each dataset's `protein1`/`protein2` columns and deduplicated, no PPI CSVs concatenated). Retrieves sequences (canonical + isoform-specific via the FASTA endpoint), GO annotations (biological process, molecular function, cellular component), and NCBI taxon IDs. Outputs `sequences.fasta`, `go_annotations.tsv`, and `species.tsv`, published to `results/_shared/data/`, then split back out per dataset (`SUBSET_FETCHED_DATA`) — see [Multiple datasets (samplesheet)](#multiple-datasets-samplesheet) below.

**FETCH_DOMAIN_META** — DDI mode only, replacing `FETCH_DATA`. Pools the Pfam family accessions across every dataset that needs a fetch and resolves them in one streaming pass over Pfam's bulk files, with no per-family API requests: `Pfam-A.fasta.gz` for the domain instances (already cut, so sequence and coordinates come from one release), `Pfam-A.clans.tsv.gz` for family → clan, and `speclist.txt` plus UniProt's reviewed-accession list for the sampling tiers. Samples up to `M` instances per family — restricted to human instances by `--instance_tier human_only` — and writes `sequences.fasta` keyed by instance, `species.tsv`, a header-only `go_annotations.tsv`, `instances.tsv` and `dropped_families.tsv` (every requested family that kept no instance, with the reason) — published to `results/_shared/data/`, then split back out per dataset (`SUBSET_DOMAIN_DATA`).

**GET_LENGTHS** — Computes per-protein sequence lengths for length-normalized BLAST scores. Runs once on the shared fetch batch (see above) for datasets needing a fetch, and once per dataset for datasets supplying a precomputed `sequences.fasta`.

**RUN_BLAST** — Runs all-against-all BLASTp with `makeblastdb` + `blastp` to quantify pairwise sequence similarity. In DDI mode the sequences are domain instances rather than full chains; the invocation is identical.

**MAKE_METIS** — Converts the BLAST results into a weighted similarity graph in METIS format. Edge weights are either raw bitscore or bitscore normalised by the geometric mean of protein lengths. In DDI mode it contracts instances onto their Pfam **clan**: instance–instance hits collapse to one weight per clan pair, taking the best hit, and hits within a clan are skipped. Clan-mates are then the same node, so no split can separate them.

**RUN_KAHIP** — Partitions the similarity graph into `k` parts using KaHIP's `kaffpa`. 

- k=3 is used togehter with **SORT_PPIS** to produce the final splits. The largest partition is used as the training set, the second-largest as validation, and the smallest as test.
- k=100 is used together with **SOLVE_ILP**. The clusters are moved to training, validation, and test set while minimizing the data loss and satisfying the constraints.

**SORT_PPIS** — Assigns each PPI to a split based on the KaHIP partition. PPIs are only retained if both partners occur in the same block of the partition. Writes per-split CSV and FASTA files. Note this path **ignores `train_split`/`val_split`/`test_split`** — the split sizes are whatever the three ranked KaHIP blocks happen to be.

**SOLVE_ILP** — Assigns each PPI to a split by solving a mixed-integer linear program (CVXPY) that maximizes the number of retained PPIs while satisfying constraints (each cluster is assigned to 1 split, and each split holds its `train_split`/`val_split`/`test_split` share of the retained interactions — 0.8/0.1/0.1 by default — within `ilp_epsilon`). A fraction of `0` is legal: that split is left out of the model and comes out empty rather than missing, so every later step still sees three splits. The ILP solver is Gurobi, by default, with the gurobi license file specified via `--gurobi_license`. If no license is available, the open-source solvers SCIP or HiGHS can be used instead.

**SPLIT_RANDOM** — `split_method=random`: a deliberately naive baseline that shuffles PPIs and slices them into train/val/test by proportion only, ignoring sequence similarity or topology entirely. Doesn't discard any PPIs, and its output skips `CDHIT2D`/`REMOVE_REDUNDANT` entirely — see [Naive baseline: the topology shortcut](#naive-baseline-the-topology-shortcut-optional) below.

**CDHIT2D** – Calls CD-HIT 2D between train/val and train/test to identify proteins in val/test that are too similar to any training protein (above the CD-HIT identity threshold). Writes a TSV of redundant proteins for each split. Only runs for `kahip`/`ilp` splits. In DDI mode it compares domain instances.

**REMOVE_REDUNDANT** — Removes proteins from val and test that are too similar to any training protein using the CD-HIT 2D TSVs. Only runs for `kahip`/`ilp` splits. Its kept-vs-removed counts feed into the same "PPI Partitioning" chart `SORT_PPIS`/`SOLVE_ILP` started (stacked `Kept` vs `Removed (CD-HIT)` for the `train`/`val`/`test` bars). In DDI mode the verdict is taken one level up and strictly: a family survives only if *every* one of its instances in that split survived, and dropping a family drops every DDI touching it.

**SELECT_EXAMPLES** — DDI mode only. Picks up to `N` domain-instance pairs ("examples") per surviving DDI under one rule: **no parent protein may be used by more than one split**. That rule couples the splits — claiming a protein for train takes it away from test — so it is solved as an ILP (CVXPY) over all three splits at once, and only over the proteins two splits actually compete for; the rest is a local pick. What is left is then cut into connected components over the contested proteins they share and each component solved on its own, which is exact — units sharing no contested protein share no constraint and no objective term — and is what keeps the model tractable at real DDI counts. A component too large for `ddi_max_ilp_candidates`, or reached after `ddi_select_max_sec` is spent, falls back to a deterministic greedy that still applies the one-split-per-protein rule exactly and loses only optimality; the **DDI Example Selection ILP** MultiQC table reports the component count, the largest component, and any fallback. Candidate pairs are any instance of family A × any instance of family B, preferring distinct parents over reusing one, and `candidate_network` pairs claim their parents in the same ILP. Emits the per-split example tables, each split's **protein universe** (the parents it claimed), the proteins no candidate ever reached (`unclaimed.txt`, plus a `{split}_reserve.txt` share of them per split, weighted by the DDIs that split kept), the DDI lists with zero-example DDIs removed, and a drop report.

**SAMPLE_NEGATIVES_DEGREE** — Samples random negative pairs for each split. By default, negatives are drawn such that each protein's degree distribution is approximately preserved, producing a balanced test set (1:1 positive:negative) and a realistic test set (1:10 ratio). With `negative_sampling_method=uniform`, endpoints are instead drawn fully uniformly at random for *every* split (not just the realistic test set) — see [Naive baseline: the topology shortcut](#naive-baseline-the-topology-shortcut-optional) below.

**SAMPLE_NEGATIVES_ILP** – An ILP-based alternative satisfying the size constraints while minimizing biases between the positive and negative sets; see [Bias-aware ILP negative sampling](#bias-aware-ilp-negative-sampling-optional) below. Under `negative_sampling_method=ilp_candidates` the same sampler draws only from the row's `candidate_network` pool, and a row may ask for several negative sets at once — see [Several negative sets from one positive split](#several-negative-sets-from-one-positive-split-optional).

**EXPAND_NEGATIVES** — DDI mode only. Both samplers draw negatives as family pairs; this step turns them into the instance pairs the classifier trains on. Positives are copied from `SELECT_EXAMPLES` rather than redrawn, and each negative family pair gets up to `N` instance pairs from that split's own protein universe, drawn by the same rule the positives were — topped up from this split's reserve of never-claimed proteins when the universe is too thin, which `SELECT_EXAMPLES` has already partitioned so no two splits can be handed the same one. Nothing is resampled to repair the ratio: a negative pair can end up with fewer than `N` examples, so the example-level ratio can drift from the family-level one. Both are reported.

**EMBED_SEQUENCES** — Computes per-protein embeddings using the selected model:
- `none` — 21-dimensional mean-pooled one-hot amino acid composition
- `esm2` — ESM-2 650M (dimension 1280), mean-pooled over residues
- `prot_t5` — ProtT5-XL (dimension 1024), mean-pooled over residues
- A path to a pre-computed `.npz` file skips this step entirely.

Every samplesheet dataset requesting the same model is embedded together in
one call over the union of their train/val/test sequences, published once to
`results/_shared/embeddings/embeddings_<model>.npz` (not duplicated per
dataset). The two neural models run in batches, grouped longest-first under a
padded-token budget (`--batch-tokens`, `--max-batch` in `bin/embed_sequences.py`);
sequences over 1024 residues are truncated, as they were one at a time.

**TRAIN_CLASSIFIER** — Trains a Random Forest classifier on concatenated pair embeddings. Hyperparameters are tuned on the validation AUROC over 3 configurations (max_depth 5/10/30, max_samples 0.2), then the best model is retrained on train+val and evaluated on the balanced and realistic test sets. When the labelled CSVs carry the DDI mode `family1`/`family2` columns it additionally averages each DDI's example predictions and scores at family level, in a second table per test split ("Classifier Performance, DDI level"). Hyperparameter selection stays on example-level validation AUROC in both modes.

**BIAS_ANALYSIS** — Runs in parallel for each attribute, computing:
- *Utility* — NMI(A; Y) = MI / √(H(A)·H(Y)): how much the attribute is correlated with the PPI label
- *Detectability* — Spearman ρ of a Ridge regressor predicting the attribute from pair embeddings

Attributes analyzed:

| Attribute                         | Description                                                                                                                                                                                                                                                                                                                                          |
|-----------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `sequence_similarity`             | BLASTp pident between the two proteins, normalized to [0, 1]                                                                                                                                                                                                                                                                                         |
| `embedding_similarity`            | Cosine similarity of the two individual protein embeddings                                                                                                                                                                                                                                                                                           |
| `functional_relatedness_BP/MF/CC` | Jaccard similarity of GO term sets (biological process / molecular function / cellular component)                                                                                                                                                                                                                                                    |
| `self_interactions`               | 1 if both proteins are identical, 0 otherwise                                                                                                                                                                                                                                                                                                        |
| `same_species`                    | 1 if both proteins share the same NCBI taxon ID, 0 otherwise (only included if the dataset contains proteins from more than one species)                                                                                                                                                                                                             |
| `topology_shortcut`               | Each endpoint's training positive-rate (pos / (pos + neg) training degree), using whichever endpoint(s) occurred in training; only included if at least one val/test/test_balanced/test_realistic pair has an endpoint that occurred in training — see [Naive baseline: the topology shortcut](#naive-baseline-the-topology-shortcut-optional) below |
| `parent_degree`                   | DDI mode only. How often each of the pair's two parent proteins is reused within the same split, averaged over the two. Positives and negatives are drawn the same way from the same protein universe, so this should carry no label information — a nonzero NMI means one class reuses its parents more than the other                                 |

DDI mode drops the three `functional_relatedness_*` attributes (Pfam families carry no GO
annotations) and adds `parent_degree`. `sequence_similarity`, `embedding_similarity` and
`same_species` are unchanged, acting on the two domain instances the row holds.
`self_interactions` and `topology_shortcut` instead act on the row's Pfam family pair (its
`family1`/`family2` columns): `self_interactions` becomes "is this a self-DDI", which also
catches a pair of two *different* instances of one family, and `topology_shortcut` becomes
family degree in the DDI graph rather than protein degree.

**COLLECT_BIAS** — Aggregates all per-attribute TSVs into a single interactive Plotly scatter plot (NMI vs detectability, colored by attribute, shaped by split).

**SIMILARITY_HEATMAP** — Plots a heatmap of pairwise BLASTp similarity between proteins in different splits, to visualize the degree of leakage.

**DDI_ATTRITION** — DDI mode only. One stacked bar per dataset accounting for every DDI that reached the splitting stage: discarded cross-cluster by the partitioner, removed by CD-HIT-2D, dropped because no domain-instance example was left for it, or kept. The counts are read back out of the DDI Partitioning and DDI Example Selection bars rather than re-derived, so the waterfall and the per-stage charts cannot disagree. The section is titled **DDI Attrition (within splitting)** because that is its scope: the four series sum to the DDIs handed to this pipeline in its `ppis` input, so DDIs a caller filtered out beforehand — by `--instance_tier human_only`, say — are outside the chart and remain the caller's to account for.

**MULTIQC** — Collects every dataset's `*_mqc.tsv`/`*_mqc.html` files into one combined report for the whole run (`results/multiqc/`). Per-attribute bias tables are excluded (the bias-scatter plot supersedes them); General Statistics, Classifier Performance, Positive vs Negative Pairs, and (in DDI mode) DDI Instance Expansion and DDI Attrition are merged across datasets (qualified sample names + an `ID` column); PPI/DDI Partitioning, DDI Example Selection, the similarity heatmap, and the bias-scatter plot remain one separate panel per dataset.

Before MultiQC runs, `relabel_mqc.py` applies the run's display order. It is the only step that sees the whole run, which is what the ordering needs: MultiQC sorts samples alphabetically and discovers custom-content sections in the order they are staged, so neither the datasets nor the sections can order themselves. It rewrites each file in place-of-copy and writes a `multiqc_config.yml` carrying `report_section_order`:

- **Sample rows** become `<NN>_<label>_<split>` — a 1-based index into (dataset order × `train, val, test, test_balanced, test_realistic, discarded`), zero-padded to the run's width. Rows that report a whole dataset rather than a split (Classifier Performance, DDI Attrition) become `<NN>_<label>`, padded to the dataset count alone. Bar-chart categories inside a single dataset's own chart keep their plain `1_train` form: only the intra-dataset order is at stake there.
- **Sections** are ordered cross-dataset summaries first (General Statistics, Positive vs Negative Pairs, Classifier Performance, DDI Instance Expansion, the ILP diagnostics, DDI Attrition), then one per-dataset kind at a time — Partitioning, DDI Example Selection, bias scatter, similarity heatmap — with the datasets inside each in the resolved order. Per-dataset section ids gain the same dataset index (`split_bar_02_minimal_leakage`); shared section ids are left exactly as they were, so external references to them survive.

`--mqc_order` sets the dataset order; `meta.mqc_labels` sets the labels. Both are report-only — see [Naming and ordering the report](#naming-and-ordering-the-report). Nothing here touches a published filename, a publish directory, or a `--split-name` value.

---

## Outputs

Every dataset from the samplesheet gets its own subtree under `--outdir`, named by its `id` column. Work shared across datasets (the deduplicated UniProt fetch and any embeddings shared by datasets requesting the same model) lives under a separate `_shared/` folder rather than being duplicated into every dataset's subtree. The combined MultiQC report for the whole run lives at the top level, `results/multiqc/`:

```
results/
├── _shared/
│   ├── data/                         # One deduplicated UniProt fetch batch (see Multiple datasets below)
│   │   └── sequences.fasta, go_annotations.tsv, species.tsv
│   │   └── instances.tsv, dropped_families.tsv   # DDI mode only (see --ddi_mode)
│   └── embeddings/
│       └── embeddings_<model>.npz    # One file per distinct embedding_model requested across datasets
├── multiqc/
│   ├── multiqc_report.html           # One combined report for the whole run
│   └── multiqc_report_data/          # MultiQC data folder
└── <id>/
    ├── multiqc/
    │   └── similarity_heatmap.html   # Standalone copy of this dataset's heatmap (also embedded in the combined report)
    ├── data/
    │   └── go_annotations.tsv        # GO annotations for this dataset's own proteins
    │   └── sequences.fasta           # FASTA for this dataset's own proteins
    │   └── species.tsv               # NCBI taxon IDs for this dataset's own proteins
    ├── similarities/
    │   └── all_vs_all.tsv            # BLAST evalue, bitscore and pident between this dataset's own proteins
    │   └── similarity.graph          # KaHIP input graph = all vs. all similarity graph (weighted edges) in METIS format
    │   └── node_mapping.tsv          # KaHIP just enumerates nodes, this maps them to protein IDs
    │   └── partitioned_proteome.txt  # KaHIP partitioned proteome (protein IDs) 
    ├── train.csv                     # Final labelled splits (positives + negatives)
    ├── val.csv
    ├── test_balanced.csv
    └── test_realistic.csv            # with 1:10 ratio of positives:negatives, negatives are uniformly sampled
```

`data/sequences.fasta` (and its `go_annotations.tsv`/`species.tsv`) is
always this dataset's own subset, even for datasets whose UniProt fetch was
folded into a shared batch with other datasets — this matters most for
`similarities/all_vs_all.tsv`, since BLAST's E-value/bitscore statistics
depend on exactly which proteins are in its search database, so it always
runs per-dataset even when the underlying sequences came from a shared fetch.

DDI mode keeps this layout and adds a few files to it — `data/instances.tsv`, an
`examples/` folder, and the instance-level `{split}_instances.csv` next to the
family-level `{split}.csv`; see [DDI mode](#ddi-mode) below.

A row asking for several negative sets suffixes each split with the set's name
(`train_ilp.csv`, `train_ilp_candidates.csv`, …) — see
[Several negative sets from one positive split](#several-negative-sets-from-one-positive-split-optional).

---

## Multiple datasets (samplesheet)

`--samplesheet` takes a CSV with one row per PPI dataset. Only `id` and `ppis`
are required — every other column may be left blank, in which case that
dataset falls back to the corresponding `nextflow.config` default. This makes
it possible to process several datasets with different parameters in a single
`nextflow run` invocation, all running in parallel:

```
id,ppis,split_method,negative_sampling_method,cdhit_identity,neg_ilp_lambda_jaccard
hippie,data/HIPPIE-current.csv,,,,
string,data/string.csv,ilp,ilp,0.5,0.5
```

| Column                                                                                                     | Overrides                   | Notes                                                                                                                                                                                                                                                   |
|------------------------------------------------------------------------------------------------------------|-----------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `id`                                                                                                       | —                           | Required. Used as the output subfolder name (`results/<id>/...`) and in logs.                                                                                                                                                                           |
| `ppis`                                                                                                     | —                           | Required. Path to this dataset's PPI CSV — or, under `--ddi_mode`, its DDI CSV (same two columns, Pfam accessions).                                                                                                                                      |
| `sequences`, `go_annotations`, `species`                                                                   | UniProt fetch step          | Supply all three to skip `FETCH_DATA` for this dataset. In DDI mode `go_annotations` is unused and the trio is `sequences`, `species`, `domain_instances` instead.                                                                                       |
| `domain_instances`                                                                                         | Pfam fetch step             | DDI mode only. A precomputed `instances.tsv`; supply it together with `sequences` and `species` to skip `FETCH_DOMAIN_META` for this dataset.                                                                                                            |
| `blast_results`                                                                                            | BLAST step                  | Supply to skip `RUN_BLAST` for this dataset (a precomputed `all_vs_all.tsv`).                                                                                                                                                                           |
| `candidate_network`                                                                                        | —                           | Candidate pool CSV (`protein1,protein2`) for the `ilp_candidates` negative set — high-confidence non-interacting *family* pairs in DDI mode, where it also feeds `SELECT_EXAMPLES`. Supplied per row only; it has no `nextflow.config` default. **Read only when the row lists `ilp_candidates`**: listing `ilp_candidates` without it is an error, and supplying it without listing `ilp_candidates` warns and ignores it everywhere. |
| `partition`, `node_mapping`                                                                                | Clustering step             | Only read when `--split_only` is set (see [Best-practice short run](#best-practice-short-run---split_only) below) — a precomputed `partitioned_proteome.txt`/`node_mapping.tsv` pair that lets that mode skip `CLUSTERING` entirely. Ignored otherwise. |
| `embedding_model`, `cdhit_identity`, `cdhit_wordsize`                                                      | `params.*` of the same name | Defaults: `embedding_model`: esm2, `cdhit_identity`: 0.4, `cdhit_wordsize`: 2                                                                                                                                                                           |
| `split_method`, `edge_weight`, `kahip_k`, `ilp_kahip_k`, `ilp_epsilon`                                     | `params.*` of the same name | Defaults: `split_method`: kahip (k=3), `edge_weight`: normalized_bitscore, `kahip_k`: 3, `ilp_kahip_k`: 100, `ilp_epsilon`: 0.05. `split_method=random` is a naive baseline, see below                                                                  |
| `train_split`, `val_split`, `test_split`                                                                   | `params.*` of the same name | Defaults: 0.8, 0.1, 0.1. Must sum to 1; any one of them may be `0` (that split then comes out empty). Only `ilp` and `random` read them — `kahip` sizes its splits from the KaHIP blocks themselves                                                       |
| `negative_sampling_method`                                                                                 | `params.*` of the same name | Defaults: default (alternatives: ilp, ilp_candidates, uniform — see below). A comma-separated list (`ilp,ilp_candidates`) asks for one negative set per entry over a single shared positive split; see [Several negative sets from one positive split](#several-negative-sets-from-one-positive-split-optional) |
| `neg_ilp_lambda_degree`, `neg_ilp_lambda_taxon_pair`, `neg_ilp_lambda_self_loop`, `neg_ilp_lambda_jaccard` | `params.*` of the same name | Only used when `negative_sampling_method` is `ilp`; see [Bias-aware ILP negative sampling](#bias-aware-ilp-negative-sampling-optional) below. Highly dataset-specific, so overridable per row rather than fixed run-wide. `--ddi_mode` forces `neg_ilp_lambda_jaccard` to 0 regardless of the row, since that term matches GO overlap and DDI mode has no GO annotations. |

Everything else (solver settings, Gurobi license, resource limits, seeds,
...) stays a run-wide default in
`nextflow.config` and is shared by every dataset in the samplesheet.

**Deduplicated work across datasets.** Datasets often overlap in which
proteins they contain, so two things are computed once per run rather than
once per dataset:
- **UniProt fetch**: every dataset that needs a fetch (doesn't supply
  `sequences`/`go_annotations`/`species`) has its unique proteins extracted
  and pooled with every other such dataset's, fetched together once, then
  split back out per dataset. BLAST still runs once per dataset on its own
  subset — see [Outputs](#outputs) above for why that matters. In DDI mode the
  same holds for the Pfam pass: one stream per run, not one per dataset.
- **Embeddings**: every dataset requesting the same `embedding_model` shares
  one embedding computation over the union of their sequences.

Both are purely a compute/storage optimization — a dataset in a mixed run
with others produces the same `train.csv`/`val.csv`/etc. it would if run
alone.

---

## Best-practice short run (`--split_only`)

For a quick, opinionated run that only produces the leakage-aware splits —
without embedding/training a baseline classifier or building the bias/MultiQC
report — skip straight to the ILP-based splitting core:

```bash
nextflow run main.nf --samplesheet samplesheet.csv --outdir results --split_only
```

`--split_only` runs only `SOLVE_ILP` → `CDHIT2D` → `REMOVE_REDUNDANT` →
`SAMPLE_NEGATIVES_ILP`, then stops — `FETCH_DATA`, `CLUSTERING` (`RUN_BLAST`/
`MAKE_METIS`/`RUN_KAHIP`), `TRAIN_BASELINE`, and `QC` never run. Because of
that, every samplesheet row must supply all of the following (no fallback to
fetching/computing them):

| Column           | Precomputed file                   |
|------------------|------------------------------------|
| `ppis`           | PPI CSV                            |
| `sequences`      | `sequences.fasta`                  |
| `go_annotations` | `go_annotations.tsv`               |
| `species`        | `species.tsv`                      |
| `partition`      | KaHIP's `partitioned_proteome.txt` |
| `node_mapping`   | `node_mapping.tsv`                 |

Under `--ddi_mode` the required set is the same except that `go_annotations`
(unused there) is replaced by `domain_instances`, and the partition is over clans
rather than proteins. `SELECT_EXAMPLES` still runs: in DDI mode the example tables
are the deliverable.

`candidate_network` remains optional, same as a normal run. `split_method`
and `negative_sampling_method` are forced to `ilp` regardless of what the
samplesheet says, since that's the only path this mode runs. Everything else
(`cdhit_identity`, `ilp_epsilon`, the `neg_ilp_*` weights, `--gurobi_license`/
`--ilp_solver`/`--neg_ilp_solver`, ...) still applies exactly as in a full run
— see [Bias-aware ILP negative sampling](#bias-aware-ilp-negative-sampling-optional)
above for those. The output is the same four files a full run produces —
`results/<id>/{train,val,test_balanced,test_realistic}.csv` — just without
the classifier/bias/MultiQC steps built on top of them.

---

## DDI mode

`--ddi_mode` splits **domain–domain interactions** — pairs of Pfam families —
instead of PPIs. Most of it is a relabelling: once the clustering node is a Pfam
**clan**, a DDI is an ordinary pairwise edge, so the partitioner, the ILP and the
negative samplers keep their exact meaning while proteins become clans and PPIs
become DDIs.

```bash
nextflow run main.nf --samplesheet samplesheet_ddi.csv --ddi_mode --outdir results -profile conda

# smoke test on the committed list of real Pfam family pairs. The first run queries
# Pfam live (~6.3 GB stream); the profile caches into <projectDir>/.pfam_cache, so
# every run after that is a stat and a read.
nextflow run main.nf -profile test_ddi,conda
```

The input is the two-column CSV from [Input PPI File](#input-ppi-file) with Pfam
accessions in it, and that is all a row needs: domain sequences, parent taxa and
the family → clan map all come from the Pfam pass in `FETCH_DOMAIN_META`.

### What it guarantees

- **One split per DDI, and clan-mates never separated.** Neither needs enforcing:
  a family sits in exactly one clan, and a clan in exactly one split.
- **No domain homology between splits** — the pipeline's usual homology barrier
  (KaHIP partition, then CD-HIT-2D), on domain instances. The CD-HIT verdict is taken per
  family and strictly: one redundant instance drops the family and every DDI
  touching it.
- **No parent protein in two splits.** `SELECT_EXAMPLES` applies this when it
  picks examples, rather than the partitioner applying it up front: as a
  partitioning constraint — "families sharing a parent must co-assign" — a handful
  of multi-domain hub proteins would chain families together transitively and pull
  most of the graph into one split.
- **Self-DDIs are kept**, positive and negative — dropping them while keeping
  self-DDI positives would make "same family ⇒ positive" a free shortcut.
- **No functional annotation.** Pfam families carry no GO terms, so the three
  `functional_relatedness_*` bias attributes and the negative sampler's Jaccard
  term go inactive, and `parent_degree` replaces them.

### Examples: `N` and `M`

A DDI's evidence is concrete domain-instance pairs, drawn purely from
co-occurrence — any instance of family A × any instance of family B, no PPI
network involved. Two counts govern that:

- **`M`** = `ddi_examples_pool_factor` × `N` instances sampled per family, in tier
  order (human-reviewed → human → reviewed → any, at random within a tier). This
  is all that BLAST, CD-HIT and the classifier ever see of a family. An empty tier
  is ordinary, not an error. `--instance_tier human_only` makes the last two tiers
  *ineligible* rather than merely least preferred (below).
- **`N`** = `ddi_examples_target` examples kept per DDI — a cap, not a quota. A
  DDI with fewer available keeps what it has; only one left with *zero* is
  dropped, and reported.

Instance ids are `family_protein_start_end`, e.g. `PF00069_P12345_10_250`, and
`data/instances.tsv` maps each one back to its family, clan, parent protein,
coordinates, taxon and source database.

`ddi_examples_pool_factor` is the single biggest cost driver in DDI mode, and not
only in `SELECT_EXAMPLES`. Raising it multiplies the instances per family, so
`FETCH_DOMAIN_META`'s reservoirs and `EMBED_SEQUENCES` grow linearly, `RUN_BLAST`
grows quadratically, and — because a protein carrying domains in two splits is
what makes it contested — the selection ILP's components grow denser. Its one
saving grace is that the per-DDI candidate pool stays capped at
`ddi_shortlist_factor` × `N`. Raise it one step at a time and check
`SELECT_EXAMPLES`'s component table and `RUN_BLAST`'s runtime each time; the
reserve of never-claimed proteins is also inert at factor 1 and live above it, so
factor ≥ 2 exercises code that factor 1 cannot reach.

### Restricting instances to human: `--instance_tier human_only`

By default (`any`) the four strata are a *preference* order: a family with no
human instance simply falls through to reviewed and then to anything. With
`--instance_tier human_only` the two non-human strata become **ineligible**, and
the consequence is the point of the option rather than a side effect:

- A family whose Pfam instances are all non-human keeps **zero** instances. Every
  DDI touching it then has no instance pair to represent it and drops out of the
  run. This is never an error and never aborts.
- Every family that *does* keep instances keeps exactly the instances it would
  have kept under `any` — the cascade fills top-down with the same room and each
  reservoir carries its own seed, so gating the lower strata cannot perturb the
  upper ones. `human_only`'s `instances.tsv` is the human subset of `any`'s, row
  for row.
- Both tiers get their own cache entry: the tier is part of
  `--interpro_cache`'s key, so a warm `any` cache cannot serve a `human_only` run.

`FETCH_DOMAIN_META` always writes **`_shared/data/dropped_families.tsv`**
(`family`, `reason`), one row per requested family that kept no instance:

| `reason` | meaning |
|---|---|
| `no_eligible_instances` | Pfam has the family, but nothing in a tier `--instance_tier` allows. Under `human_only` this is the expected bulk; under `any` it should be empty |
| `dead` | the accession is listed in `Pfam-A.dead` |
| `not_in_pfam` | no instances in `Pfam-A.fasta` and not listed as dead — usually a typo in the input |

The three are counted and warned separately on `FETCH_DOMAIN_META`'s stderr, so a
large `human_only` drop cannot hide among dead accessions. Note this compounds
with `ddi_examples_pool_factor`: asking for `M = 25` instances from human strata
alone will underfill most families, and fewer instances per family is *good* for
val/test survival (see the factor's effect above) but leaves smaller pools for
`SELECT_EXAMPLES`.

`species.tsv` and the `same_species` bias attribute are deliberately untouched —
under `human_only` that attribute goes constant and its NMI to ~0, which is the
intended sanity signal rather than something to suppress.

### Parameters

| Parameter                  | Default | Description                                                                                                                                     |
|----------------------------|---------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| `ddi_mode`                 | `false` | Interpret the interaction file's two columns as Pfam family accessions                                                                          |
| `ddi_examples_target`      | `5`     | `N`, the cap on examples kept per DDI                                                                                                           |
| `ddi_examples_pool_factor` | `5`     | `M` = this × `N`, the instances sampled per family                                                                                              |
| `ddi_select_max_sec`       | `300`   | Total `SELECT_EXAMPLES` ILP budget, shared across the independent components in proportion to their size. A stage that runs out keeps the best solution found; components reached after the whole budget is gone go to the greedy fallback. Both are warned and counted |
| `ddi_max_ilp_candidates`   | `200000`| A component with more candidate pairs than this skips the ILP for the greedy fallback, so one oversized component cannot exhaust memory during canonicalisation |
| `ddi_lambda_diversity`     | `0.1`   | How strongly a DDI's examples prefer distinct parents (`P1-P2, P3-P4` over `P1-P2, P1-P3`). Must stay below 0.5, so it never costs a DDI an example |
| `ddi_shortlist_factor`     | `4`     | Cap on a DDI's candidate pool before the ILP, as a multiple of `N`. A no-op at `M = N`, a guard for a larger pool                               |
| `ddi_candidate_factor`     | `4`     | Cap on `candidate_network` pairs per split, as a multiple of that split's DDI count                                                             |
| `ddi_select_verbose`       | `false` | Let the `SELECT_EXAMPLES` solver print its own log to `.command.err`. Off by default: one block per ILP solve, which is large at real DDI counts |
| `instance_tier`            | `any`   | `any` fills the four strata in preference order; `human_only` makes the two non-human strata ineligible, so a family with no human instance keeps zero instances and its DDIs drop out. Reported in `_shared/data/dropped_families.tsv` |
| `pfam_fasta`               | `null`  | Local `Pfam-A.fasta.gz`, skipping the ~6.3 GB download                                                                                          |
| `pfam_clans`               | `null`  | Local `Pfam-A.clans.tsv(.gz)`, skipping that download                                                                                           |
| `pfam_release`             | `null`  | Pin the release string (e.g. `38.2`) instead of downloading `Pfam.version`. That lookup happens before the cache directory is known, so it is the one fetch a warm cache cannot skip |
| `interpro_cache`           | `null`  | Directory for the cached downloads and sampled instances. **Must be an absolute path on a filesystem all compute nodes share** — it is resolved to one, but node-local scratch gives every task its own cold cache. A convenience only: a cold run produces identical output. `-profile test_ddi` sets it to `<projectDir>/.pfam_cache` |

`SELECT_EXAMPLES` reuses `--ilp_solver` and `--gurobi_license` rather than adding
its own; every other parameter (`cdhit_*`, `ilp_*`, `neg_ilp_*`, the split
fractions) keeps its usual meaning.

### Extra outputs

DDI mode publishes the same tree as PPI mode (see [Outputs](#outputs)) plus:

```
results/<id>/
├── data/
│   └── instances.tsv                            # instance -> family, clan, parent protein, coords, taxon
├── examples/
│   ├── {train,val,test}_sel.csv                 # the split's DDIs, zero-example ones removed
│   ├── {train,val,test}_examples.csv            # the selected instance pairs
│   ├── {train,val,test}_candidate_examples.csv  # the same, for candidate_network negatives
│   ├── {train,val,test}_universe.txt            # the parent proteins this split claimed; no other split may use them
│   ├── {train,val,test}_reserve.txt             # this split's share of the spare pool, weighted by the DDIs it kept
│   └── unclaimed.txt                            # parents no candidate example reached -- the spare pool, unpartitioned
├── {train,val,test_balanced,test_realistic}.csv            # family-level labelled pairs
└── {train,val,test_balanced,test_realistic}_instances.csv  # instance-level -- what the classifier trains on
```

A row asking for several negative sets suffixes both of those last two lines with
the set's name, except `test_realistic`, which stays one shared file — see
[Several negative sets from one positive split](#several-negative-sets-from-one-positive-split-optional).
Everything above them, `examples/` included, is produced once per row whatever the
negative sets are.

The instance-level files carry `protein1,protein2,label` plus `family1,family2`.
Those two extra columns are how `TRAIN_CLASSIFIER` and `BIAS_ANALYSIS` recognise
DDI mode — neither takes a flag for it.

### Reading the report

The MultiQC report gains "DDI Partitioning", "DDI Example Selection", "DDI
Example Selection ILP", "DDI Instance Expansion" and "DDI Attrition" — the last a
single stacked bar per dataset accounting for every input DDI (discarded
cross-cluster, removed by CD-HIT-2D, dropped for want of an example, or kept).
"DDI Example Selection ILP" describes the solve itself; the two numbers to check
there are **Greedy fallback components**, which should be 0, and **Largest
component (candidates)**, which is what decides whether the decomposition is
still doing its job as the input grows. Each test split gets two
classifier tables, one per example and one per DDI, the latter averaging a DDI's
example predictions before scoring. On the per-DDI table read **AUROC and
AUPRC**: averaging `N` near-chance probabilities pulls every DDI toward the same
mean, so a fixed 0.5 cut — and F1/MCC/precision/recall with it — says little until
the model clears chance.

`bin/other/check_ddi_invariants.py --results results/<outdir>` re-checks the whole
published tree: no parent protein and no family in two splits, no DDI over `N`,
and every id really an instance of the family it claims.

### Worth knowing

- `split_method=random` is the leaky baseline here too: it skips CD-HIT **and** the
  one-protein-per-split rule, so families and parent proteins may straddle splits.
  That is exactly the leakage the baseline exists to show.
- `cdhit_identity` 0.4 and `cdhit_wordsize` 2 are CD-HIT's own floor, so on ~100 aa
  domains no stricter homology cut is reachable through CD-HIT. The strict
  per-family verdict above is what compensates.
- BLAST runs as in PPI mode — no `-evalue`, hence `blastp`'s default of 10 — which
  keeps the two modes comparable but over-connects the clan graph on short
  sequences. That errs in the safe direction: KaHIP then separates harder, so more
  DDIs are discarded as cross-cluster and less leakage survives.
- Accessions Pfam has killed between releases are named in the fetch report and
  their DDIs drop out; one unresolvable accession never aborts a run.
- Positives and negatives are drawn the same way, from the same per-split protein
  universe, so DDI mode avoids PPI mode's asymmetry — there, positives come from
  real complexes and negatives from pairs merely not known to interact, which is
  itself learnable. `parent_degree` is the check: a nonzero NMI means one class
  reuses its parents more than the other.

---

## Naive baseline: the topology shortcut (optional)

Many PPI-splitting publications randomly split pairs 80/10/10 and sample
negatives uniformly at random. This inflates reported performance: positive
degree follows a power law while uniform-random negatives don't, so a
protein's *training* degree — specifically what fraction of its training
interactions are positive, `pos_degree_in_training(p) /
(pos_degree_in_training(p) + neg_degree_in_training(p))` — becomes
predictive of the label by itself, especially for proteins at the extremes
of the degree distribution. A model can exploit this "topology shortcut"
instead of learning any real interaction signal, because a naive random
split lets the same protein land in both training and test.

Reproduce this deliberately-bad setup with:

```bash
nextflow run main.nf --split_method random --negative_sampling_method uniform
```

(or per-dataset via the samplesheet's `split_method`/`negative_sampling_method`
columns, see [Multiple datasets](#multiple-datasets-samplesheet) above).
`split_method=random` (`SPLIT_RANDOM`) shuffles PPIs and slices them by
proportion only — no homology- or topology-aware partitioning — and,
critically, **skips `CDHIT2D`/`REMOVE_REDUNDANT` entirely**: since those
steps treat a protein shared between train and test as trivially
~100%-self-similar and strip it back out of test, running them here would
erase almost all of the train/test protein overlap this baseline exists to
demonstrate. `negative_sampling_method=uniform` draws negatives fully
uniformly at random for every split (not just the realistic test set),
matching how the "many publications" this baseline is modeled on sample
negatives.

The `topology_shortcut` bias attribute (see the Attributes table above)
quantifies exactly this effect: it's only computed for pairs where at least
one endpoint occurred in training, and the whole attribute is skipped (no
MultiQC output at all) when *no* val/test/test_balanced/test_realistic pair
qualifies — which is the normal, expected case for the real `kahip`/`ilp`
splits, since those keep train/val/test protein sets disjoint by design.
Run a `random`-split dataset alongside a `kahip`/`ilp`-split dataset on the
same PPI data in one samplesheet to see the difference directly: the naive
baseline's report will show a `topology_shortcut` table with elevated
NMI/detectability that the leakage-aware splits' reports won't show at all.

---

## Bias-aware ILP negative sampling (optional)

`SAMPLE_NEGATIVES_ILP` is an opt-in alternative to `SAMPLE_NEGATIVES` that chooses
the negative set by solving a mixed-integer linear program (CVXPY), rather than
sampling at random. It matches per-protein degree, per-taxon-pair interaction
counts, self-interaction counts, and mean GO-BP Jaccard similarity between the
positive and negative sets. `bin/sample_negatives_ilp.py` samples exactly one
split per invocation; the process runs once per split (train, val,
test_balanced) and Nextflow executes all three in parallel. `test_realistic`
is deliberately excluded — same as under `negative_sampling_method=default`,
it always gets uniform-at-random negatives (via `SAMPLE_NEGATIVES_DEGREE`)
regardless of `negative_sampling_method`, since the point of that split is to
simulate an uncontrolled random screen, and bias-matching its negatives to
the positives would defeat that purpose. Together they produce the same four 
output files (`train.csv`, `val.csv`, `test_balanced.csv`, 
`test_realistic.csv`) as the default sampler, so all downstream steps are 
unaffected. 

Enable it with:

```bash
nextflow run main.nf --negative_sampling_method ilp
```

The `--neg_ilp_lambda_*` weights are highly
dataset-specific (they depend on each dataset's degree distribution, taxon
composition, and GO annotation coverage), so when running multiple datasets
via `--samplesheet` they are set **per row**, not as a single run-wide value —
see the samplesheet column reference in
[Multiple datasets (samplesheet)](#multiple-datasets-samplesheet) above. The
values below are just the `nextflow.config` fallback used for any row that
leaves them blank.

| Parameter                   | Default   | Description                                                                                                                                                              |
|-----------------------------|-----------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `negative_sampling_method`  | `default` | `default` (random), `ilp`, `ilp_candidates` (`ilp` restricted to `candidate_network`), or `uniform`. A comma-separated list asks for one negative set per entry — see [Several negative sets from one positive split](#several-negative-sets-from-one-positive-split-optional) |
| `candidate_network`         | —         | CSV (`protein1,protein2`) restricting the candidate pool, e.g. a Negatome database or a topology-driven pool. Samplesheet column only — there is no `nextflow.config` default. Read only by the `ilp_candidates` negative set; recommended for large protein universes (see below). |
| `neg_ilp_lambda_degree`     | `1.0`     | Weight of the per-protein aggregate degree-matching term                                                                                                                 |
| `neg_ilp_lambda_taxon_pair` | `1.0`     | Weight of the global taxon-pair matching term                                                                                                                            |
| `neg_ilp_lambda_self_loop`  | `1.0`     | Weight of the self-interaction count matching term                                                                                                                       |
| `neg_ilp_lambda_jaccard`    | `1.0`     | Weight of the mean GO-BP Jaccard matching term. Forced to 0 under `--ddi_mode` (no GO annotations), which skips the term outright rather than matching an all-zero one     |
| `neg_ilp_solver`            | `auto`    | `auto`, `gurobi`, `scip`, or `highs`. `auto` tries Gurobi first, then falls back to an open-source solver.                                                               |
| `neg_ilp_time_limit`        | `7200`    | Solver time limit in seconds                                                                                                                                             |
| `neg_ilp_mip_gap`           | `0.01`    | Solver MIP gap tolerance                                                                                                                                                 |
| `gurobi_license`            | `"gurobi.lic"` | Path to a Gurobi license file, only used if the `gurobi` solver is selected. It is staged with `checkIfExists`, so a nonexistent path fails the run at channel construction |

For Gurobi, install it yourself and point `--gurobi_license` at your license
file (`pip install gurobipy` is already pulled in by `environment.yml`). With
no usable Gurobi license, `auto` falls back to SCIP or HiGHS, both installed
by `environment.yml`.

The default candidate pool is the full upper-triangle complement of the
positive set, which is quadratic in the number of proteins (~1.1×10⁸ pairs for
15k proteins). If it exceeds the process's candidate cap, the pool is randomly
subsampled instead of built in full — weighted toward each protein's own
negative-degree cap rather than uniformly, so a low-degree protein doesn't end
up with far more candidates than it can ever use. That cap is not a parameter:
`SAMPLE_NEGATIVES_ILP` sets it to 4× the split's positive count (raised on each
retry attempt) and passes it to `bin/sample_negatives_ilp.py --max-candidates`.
For large datasets, supply `candidate_network` to restrict the pool to a
deliberately-curated set instead (e.g. a Negatome database or a topology-driven
pool); if that network is itself larger than the cap, it's capped down the same
degree-weighted way rather than used in full.

Each split's process writes its own `<split>_mqc.tsv` diagnostics row and,
optionally, `<split>_residuals_mqc.tsv` (per-protein degree residuals); MultiQC
picks up all of them.

In DDI mode both samplers run exactly as described here, one level up: the pairs
they draw are Pfam family pairs on the family DDI graph, and `EXPAND_NEGATIVES`
then turns them into domain-instance pairs.

---

## Several negative sets from one positive split (optional)

`negative_sampling_method` takes a comma-separated list, one entry per negative
set wanted:

```bash
nextflow run main.nf --negative_sampling_method ilp,ilp_candidates
```

Every entry names a sampler — `default`, `ilp`, `ilp_candidates`, or `uniform` —
and the entries share **one** positive split. Splitting, redundancy removal and
(in DDI mode) `SELECT_EXAMPLES` all still run once per row, so the positive rows
of the resulting datasets are identical by construction, not by luck: the sets
differ only in their negatives. That is the point of the feature — comparing two
negative-sampling strategies with the positives held fixed.

`ilp_candidates` is `SAMPLE_NEGATIVES_ILP` restricted to the row's
`candidate_network` pool; `ilp` is the same sampler unrestricted, and does not
read the network even when one is supplied. Validation, at channel construction:

- `ilp_candidates` listed with no `candidate_network` → **error**.
- `candidate_network` supplied but `ilp_candidates` not listed → **warning**, and
  the network is ignored everywhere, `SELECT_EXAMPLES` included.
- the same method listed twice → **error**.

**Sizing the candidate network.** `ilp_candidates` can only draw pairs whose *both*
endpoints are in the split it is sampling, because nodes are split-exclusive. Survival
is therefore quadratic in a split's share: a network of `P` pairs spread evenly over the
node set leaves roughly `P · f²` usable pairs in a split holding fraction `f` of the
nodes. `SAMPLE_NEGATIVES_ILP` errors out rather than under-fill —

```
RuntimeError: val: need 11 negatives but only 7 candidate pairs are available.
Supply a larger --candidate-network or lower the negative ratio.
```

— so size the network by the *smallest* split, not by the dataset. A split with `N`
positives needs on the order of `sqrt(4N)` distinct nodes represented in the network to
cover itself. `ilp` is unaffected: it draws from the full non-positive complement.

**Filenames.** One entry leaves every output name exactly as it was
(`train.csv`, `test_realistic_instances.csv`, …), so existing samplesheets are
unaffected. Several entries suffix each split with the set's name:

```
results/<id>/
├── train_ilp.csv                 ├── train_ilp_candidates.csv
├── val_ilp.csv                   ├── val_ilp_candidates.csv
├── test_balanced_ilp.csv         ├── test_balanced_ilp_candidates.csv
└── test_realistic.csv            (one shared file, see below)
```

`test_realistic` is the one exception. Its negatives are drawn uniformly at
random for *every* method (the ILP path excludes that split by design, and
`uniform`/`default` both sample it uniformly), so over one shared positive split
with one seed its content cannot depend on the negative set. It is therefore
sampled — and in DDI mode expanded — once, published unsuffixed, and reported as
belonging to every set. The pipeline logs one line per affected dataset saying
so.

**MultiQC.** With several sets, the negative-sampling, classifier and bias
sections tag their samples `<id>_<negset>` so the sets do not overlay each other
in one chart; the splitting-stage sections and the DDI attrition waterfall stay
per-dataset, because that stage runs once per row — and take the *first* set's
label, since a row produces one of them and several labels. `test_realistic`'s
negative-sampling row likewise appears once, under that same first-set label; the
file on disk stays the unsuffixed `test_realistic.csv` either way. Each of a row's
sets is a separate entry in `--mqc_order`, and `meta.mqc_labels` can give each of
them its own display name; see
[Naming and ordering the report](#naming-and-ordering-the-report).

`bin/other/check_ddi_invariants.py` understands the suffixed layout and adds one
check for it: every negative set of a row must carry the same positive rows in
every split.

---

## Naming and ordering the report

Two knobs control how the MultiQC report reads. Both affect the report only: no
published filename, no publish directory and no `--split-name` value depends on
either, so a downstream consumer that joins on paths is unaffected.

### `--mqc_order`: which dataset comes first

MultiQC sorts samples alphabetically, which puts the datasets in alphabetical
order of their id — never the order a reader wants. `--mqc_order` names them
instead:

```bash
nextflow run main.nf --samplesheet s.csv \
    --mqc_order 'random,minimal_leakage,external_test'
```

In a config file it is a list (`mqc_order = ['random', 'minimal_leakage']`); on
the command line, a comma-separated string. The names are **display labels**, one
per (samplesheet row, negative set) — so a row asking for two negative sets
contributes two names. Every rule here is advisory, warned about on stderr and in
`.nextflow.log`, and never fatal:

- a name matching no dataset is ignored, and the labels that *do* exist are listed
  so a typo is obvious;
- a name given twice keeps its first position;
- a dataset the list does not mention is appended after the listed ones,
  alphabetically — so a partial order is still a total one;
- the default, `[]`, is alphabetical.

`bin/relabel_mqc.py` then turns the resolved order into the numeric prefixes
MultiQC sorts on, and into the `report_section_order` config that fixes the
section order — see the **MULTIQC** step above for the exact shapes.

### `meta.mqc_labels`: what each dataset is called

A row asking for several negative sets gets a `_<negset>` suffix on its display
label, so the report says `minimal_leakage_ilp_candidates` where the pipeline that
drove the run may call that dataset something else entirely. `meta.mqc_labels` is
an optional per-dataset map from negative-set name to display label:

```groovy
[ id: 'minimal_leakage',
  negative_sampling_method: 'ilp,ilp_candidates',
  mqc_labels: [ 'ilp': 'minimal_leakage', 'ilp_candidates': 'minimal_leakage_hcni' ],
  /* ...every other meta key... */ ]
```

It is only reachable from an including pipeline — there is no samplesheet column
for it, which is exactly why a standalone run's labels, and therefore every task
hash, are unchanged. An absent key, a non-map value, or a negative set the map
does not cover falls back to `<id><negset suffix>`; a map whose keys do not match
the row's negative sets is warned about and the uncovered sets fall back
individually. Like every other `meta` key it must be set when `meta` is built and
never mutated, because every `join()`/`combine(by: 0)` in the pipeline keys on the
whole map.

Artefacts that belong to the whole row rather than to one negative set — the
partitioning bars, the DDI example-selection tables, the similarity heatmap, the
DDI attrition waterfall — run once per row and take the **first** negative set's
label.

Putting the two together, for a run of three rows and five negative sets:

```groovy
mqc_labels: ['uniform': 'random']                                                  // row 1
mqc_labels: ['ilp': 'minimal_leakage', 'ilp_candidates': 'minimal_leakage_hcni']    // row 2
mqc_labels: ['ilp': 'external_test',   'ilp_candidates': 'external_test_hcni']      // row 3

mqc_order = ['random', 'minimal_leakage', 'minimal_leakage_hcni',
             'external_test', 'external_test_hcni']
```

5 labels × 6 splits = 30, so the width is 2 and the rows run
`01_random_train` … `30_external_test_hcni_discarded`.

---

## Standalone STRING channel analysis

To investigate which STRING evidence channels explain classifier performance differences between datasets, use the standalone script (not part of the Nextflow pipeline):

```bash
python bin/analyse_string_channels.py \
    --train      results/<id>/train.csv \
    --test       results/<id>/test_balanced.csv \
    --embeddings results/_shared/embeddings/embeddings_<model>.npz \
    --out        string_channel_analysis.tsv
```

This fits a Ridge regressor (on positive pairs only) to predict each STRING evidence channel score from pair embeddings, and reports train and test Spearman ρ per channel. `combined_score` is excluded since it is derived from the individual channels.

---

## Embedding this pipeline in another Nextflow pipeline

The whole pipeline is the named workflow `PPI_SPLITTING` in `main.nf`; the anonymous
`workflow { }` entry is a thin caller that builds the dataset channel from
`--samplesheet` and nothing else. An including pipeline therefore skips the
samplesheet entirely and builds the channel itself:

```groovy
include { PPI_SPLITTING } from './subworkflows/external/ppi-splitting/main.nf'

// tuple(meta, filesMap) -- one item per dataset
datasets_ch = channel.of(
    tuple(
        [ id: 'minimal_leakage', split_method: 'ilp',
          // one entry per negative set wanted out of this row's single positive split
          negative_sampling_method: 'ilp,ilp_candidates',
          train_split: 0.7, val_split: 0.1, test_split: 0.2, /* ...every other meta key... */ ],
        [ ppis: file('3did.csv'), sequences: file('sequences.fasta'),
          species: file('species.tsv'), domain_instances: file('instances.tsv'),
          candidate_network: file('candidate_network.csv'),
          go_annotations: [], blast_results: [],
          partition: [], node_mapping: [] ]
    )
)

out = PPI_SPLITTING(datasets_ch)
```

Three things to get right:

1. **`meta` must carry every key** `buildDatasetsChannel()` sets — the subworkflows
   read `meta.split_method`, `meta.cdhit_identity` and the rest directly, and a
   missing key surfaces as a null in a rendered command line, not as an error.
   `meta` must also not be mutated afterwards: every `join()`/`combine(by: 0)` in
   the pipeline keys on the whole map. One key is available *only* here and is
   optional: `mqc_labels`, which names the MultiQC rows in your own vocabulary —
   see [Naming and ordering the report](#naming-and-ordering-the-report).
2. **`filesMap` must have all nine keys.** An absent optional file is `[]`, never
   `null` — a `path` input accepts `[]` as "no file" and `null` breaks staging.
3. **Include `conf/params.config`**, before your own `params { }` block:

   ```groovy
   includeConfig 'subworkflows/external/ppi-splitting/conf/params.config'
   ```

   Nextflow reads only the root project's `nextflow.config`, so without this every
   `params.*` this pipeline reads is undefined. Order matters because a later
   assignment wins and `outdir` is defined on both sides (`seed` too, with the same
   meaning). Those two are the only collisions.

### What it emits

| emit | shape |
|---|---|
| `instances` | `tuple(meta, instances.tsv)` — `tuple(meta, [])` in PPI mode |
| `sequences` | `tuple(meta, sequences.fasta)` |
| `labelled` | `tuple(meta, negset, label, csv)` — labelled pairs at node level (Pfam family in DDI mode) |
| `labelled_inst` | the same at domain-instance level; empty in PPI mode |
| `multiqc_report` | `multiqc_report.html`; empty under `--split_only` |

`negset` is one of the negative-sampling methods the row asked for, carried as its
own tuple field rather than inside `meta` — putting it in `meta` would rekey every
`join()`/`combine(by: 0)` in the pipeline. A row whose
`negative_sampling_method` lists several methods emits one item per
`(negset, label)`, all over the same positive rows; see
[Several negative sets from one positive split](#several-negative-sets-from-one-positive-split-optional),
and note that `test_realistic` is one shared file emitted once per `negset`.
`label` is one of `train`, `val`, `test_balanced`, `test_realistic`. There is no
`versions` channel — this pipeline does not capture tool versions anywhere.

`negative_sampling_method` is validated inside `PPI_SPLITTING`, not in
`buildDatasetsChannel()`, so a channel you build in Groovy is held to the same
rules (unknown or duplicated method, `ilp_candidates` without a
`candidate_network`) rather than failing later inside a task.

Everything is still published to `--outdir` exactly as in a standalone run; the
emits exist so an including pipeline can ingest or re-publish without knowing this
pipeline's layout.

### Which profile

Use **`-profile docker`**. Nextflow puts only the *root* project's `bin/` on
`PATH`, and when this pipeline is included the root project is yours, so
`sample_negatives.py` and friends would not resolve. The image bakes `bin/` in, so
processes work identically standalone and embedded — which is also why **the image
tag and the submodule tag have to move together**: a `bin/` change with a stale
image is a silently wrong run, not a failure. `-profile conda` is supported for
standalone runs only, for exactly that `$projectDir/bin` reason.

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.10
- Conda (for the environment) — or install the packages in `environment.yml` manually, or use `-profile docker` and the image built from `Dockerfile`
- Internet access for the initial UniProt fetch (subsequent runs use cached Nextflow work directories)
- A GPU is recommended but not required for `esm2` and `prot_t5` embedding models. It is required under `-profile gpu`, which turns an unusable CUDA device into an error instead of a slow CPU run. The wheel's CUDA major version must not exceed the driver's (minors are compatible). `Dockerfile` pins `torch==2.10.0+cu128` for that reason; **`environment.yml` does not pin torch at all**, so a conda env built from it takes whatever pip resolves that day and can reproduce the silent-CPU failure `-profile gpu` exists to catch
- For DDI mode, internet access for the Pfam pass — one ~6.3 GB transfer per run, or none if `--pfam_fasta`/`--pfam_clans` point at local copies. `--interpro_cache <abs-dir>` makes repeat runs a stat and a read (plus one small `Pfam.version` request, which `--pfam_release` removes)
