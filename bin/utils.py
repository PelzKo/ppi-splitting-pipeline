#!/usr/bin/env python3
"""Shared I/O utilities for PPI pipeline scripts."""

import csv
import heapq
import sys
from collections import defaultdict

import numpy as np

# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------


def read_fasta(path):
    """Return {protein_id: sequence} from a FASTA file."""
    seqs = {}
    acc, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(parts)
                acc, parts = line[1:].split()[0], []
            elif line:
                parts.append(line)
    if acc and parts:
        seqs[acc] = "".join(parts)
    return seqs


def write_fasta(seqs, proteins, path):
    """Write proteins from seqs to a FASTA file, sorted by ID."""
    with open(path, "w") as fh:
        for p in sorted(proteins):
            if p in seqs:
                fh.write(f">{p}\n{seqs[p]}\n")


# ---------------------------------------------------------------------------
# PPI CSV I/O
# ---------------------------------------------------------------------------


def read_ppis(path):
    """Return list of row dicts from a PPI CSV, stripping protein ID whitespace."""
    seen, rows = set(), []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            row["protein1"] = row["protein1"].strip()
            row["protein2"] = row["protein2"].strip()
            p1, p2 = row["protein1"], row["protein2"]
            key = (min(p1, p2), max(p1, p2))
            if key not in seen:
                seen.add(key)
                rows.append(row)
    print(f"  {len(rows):,} unique PPIs", file=sys.stderr)
    return rows


def write_ppi_csv(rows, path):
    """Write PPI row dicts to CSV, preserving all columns from the input."""
    fieldnames = list(rows[0].keys()) if rows else ["protein1", "protein2"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# Domain instance table I/O (DDI mode)
# ---------------------------------------------------------------------------

INSTANCE_COLUMNS = ["instance_id", "family", "clan", "protein_id", "start", "end", "taxon_id", "source_db"]


def read_instances(path):
    """Return list of row dicts from an instances.tsv (bin/fetch_domains.py's output).

    One table serves every DDI-mode consumer: MAKE_METIS reads
    instance_id -> clan, the splitters invert to family -> {instance_id} and
    read family -> clan, SELECT_EXAMPLES reads protein_id. Callers build the
    dict they need -- keeping a single reader is what stops two views of the
    same file from disagreeing.
    """
    rows = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in INSTANCE_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path} is missing instances.tsv column(s): {', '.join(missing)}")
        for row in reader:
            rows.append({c: row[c].strip() for c in INSTANCE_COLUMNS})
    return rows


def instances_by_family(rows):
    """Invert read_instances() output to {family: {instance_id}}."""
    members = {}
    for row in rows:
        members.setdefault(row["family"], set()).add(row["instance_id"])
    return members


def instance_parent(instance_id):
    """Parent protein accession of a domain instance id, or the id itself if it does not parse.

    bin/fetch_domains.py builds instance ids as f"{family}_{protein}_{start}_{end}"
    -- underscore-separated rather than pipe-separated on purpose, since a "|"
    would be parsed as an NCBI defline by `makeblastdb -parse_seqids` and then
    undone by the split("|") in make_metis.py. Neither a Pfam accession nor a
    UniProt accession contains an underscore, so a right-anchored split recovers
    all four fields. This is the only place that parse lives; every consumer that
    has the instances.tsv to hand should read its protein_id column instead.
    """
    parts = instance_id.rsplit("_", 3)
    return parts[1] if len(parts) == 4 else instance_id


EMPTY = frozenset()


def pair_candidates(fam1, fam2, members, available):
    """Every instance pair that could represent the interaction (fam1, fam2).

    Pure co-occurrence, per the design: any instance of fam1 against any
    instance of fam2, with no PPI network consulted. `available` restricts to
    instances that actually reached this split's FASTA -- instances.tsv can list
    a family member whose sequence never arrived, and an example with no
    sequence cannot be embedded.

    A self-interaction (fam1 == fam2) yields unordered pairs *including* i == i,
    which is the direct analogue of PPI mode's (P, P) rows and is what keeps a
    single-instance self-DDI representable at all. The diversity term
    deprioritises them, since such a pair spends its parent's budget twice.

    Shared by SELECT_EXAMPLES (positives, candidate_network negatives) and
    EXPAND_NEGATIVES (sampled negatives), so both levels of DDI mode draw their
    instance pairs by exactly the same rule.
    """
    # The default has to be a set, not (): a family can reach a split CSV without
    # ever appearing in instances.tsv -- a dead or mistyped Pfam accession that
    # FETCH_DOMAIN_META could not resolve is the ordinary case, and the random
    # splitter passes it straight through where solve_ilp/sort_ppis drop it for
    # having no clan. Such a family simply has no candidates, which the caller
    # already reports as a DDI dropped for want of an example.
    la = sorted(members.get(fam1, EMPTY) & available)
    if fam1 == fam2:
        return [(a, b) for i, a in enumerate(la) for b in la[i:]]
    lb = sorted(members.get(fam2, EMPTY) & available)
    return [(a, b) for a in la for b in lb]


def diverse_pick(pairs, k, parent_of, rng):
    """Pick <= k pairs, preferring ones whose parent proteins are still unused.

    Realises the design's diversity preference -- P1-P2, P3-P4, P5-P6 over
    P1-P2, P1-P3, P1-P4 -- and doubles as SELECT_EXAMPLES' shortlist trimmer.
    Deterministic given `rng`.

    Each pick takes the lowest (reuse score, shuffled index), where the score is
    used[parent_a] + used[parent_b]. A plain min() over the remaining pool costs
    O(k * n) with two dict lookups per comparison, which is 20 * 225 comparisons
    per DDI at --ddi_examples_pool_factor 3 and, over ~90k units, minutes of
    single-threaded Python before the solver starts. A pick only changes the score
    of pairs sharing one of the two parents just used, so the scan is replaced by a
    lazy heap: push a pair's score when it changes, skip entries whose score is
    stale on pop. Scores only ever rise, so a stale entry always sorts before its
    replacement and is discarded rather than hiding it -- the (score, index)
    ordering, and therefore the output, is identical to the min() scan at
    O((n + k * touched) log n).
    """
    if k <= 0 or not pairs:
        return []
    order = sorted(pairs)
    rng.shuffle(order)

    touching = defaultdict(list)
    for i, (a, b) in enumerate(order):
        pa, pb = parent_of[a], parent_of[b]
        touching[pa].append(i)
        if pb != pa:
            touching[pb].append(i)

    score = [0] * len(order)
    taken = [False] * len(order)
    heap = [(0, i) for i in range(len(order))]  # already in heap order
    used, chosen = defaultdict(int), []
    while len(chosen) < k and heap:
        s, i = heapq.heappop(heap)
        if taken[i] or s != score[i]:
            continue
        taken[i] = True
        a, b = order[i]
        chosen.append(order[i])
        # A pair whose two instances share a parent spends that parent twice,
        # which is the double reuse it is -- so increment per endpoint, not per
        # distinct parent.
        used[parent_of[a]] += 1
        used[parent_of[b]] += 1
        for p in {parent_of[a], parent_of[b]}:
            for j in touching[p]:
                if not taken[j]:
                    score[j] = used[parent_of[order[j][0]]] + used[parent_of[order[j][1]]]
                    heapq.heappush(heap, (score[j], j))
    return chosen


def expand_members(groups, members):
    """Map group ids to the union of their members; identity when `members` is falsy.

    The split CSVs and the split FASTAs speak two different vocabularies in DDI
    mode: the CSV's protein1/protein2 hold Pfam *families* while the FASTA is
    keyed by domain *instance*. write_fasta() skips ids it does not know
    silently, so handing it families produces three empty FASTAs and no
    warning -- this is the conversion that prevents it. In PPI mode `members`
    is None, the two vocabularies coincide and this is a plain copy.
    """
    if not members:
        return set(groups)
    return {m for g in groups for m in members.get(g, ())}


def strict_survivors(present, keep, members):
    """Groups all of whose PRESENT members survived; identity when `members` is falsy.

    `present` are the sequence ids in one split's FASTA and `keep` the subset
    CD-HIT-2D reported as non-redundant, so in PPI mode this is exactly the
    `present & keep` the redundancy filter has always applied. In DDI mode a
    family is kept only if every instance of it that reached this split
    survived: STRICT rather than a majority vote, which is what keeps the CSV
    and the FASTA consistent for free -- a kept family has no lost instance, so
    expand_members() can never reintroduce a sequence CD-HIT flagged.
    """
    if not members:
        return present & keep
    lost = present - keep
    return {g for g, ms in members.items() if (ms & present) and not (ms & lost)}


# ---------------------------------------------------------------------------
# KaHIP partition I/O
# ---------------------------------------------------------------------------


def read_node_mapping(path):
    """Return {node_id (int, 1-indexed): protein_id}."""
    mapping = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mapping[int(row["node_id"])] = row["protein_id"]
    return mapping


def read_partition(path):
    """Return list where element i is the partition of node (i+1)."""
    partitions = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                partitions.append(int(line))
    return partitions


def read_labelled_csv(path):
    """Return (pairs, labels) from a CSV with protein1, protein2, label columns."""
    pairs, labels = [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            pairs.append((row["protein1"].strip(), row["protein2"].strip()))
            labels.append(int(row["label"]))
    return pairs, np.array(labels)


def read_family_pairs(path):
    """Per-row order-independent family pair, or None if the CSV has no family columns.

    DDI mode's labelled CSVs (EXPAND_NEGATIVES' output) are domain-*instance*
    pairs carrying the family pair they represent in family1/family2, so several
    rows belong to one DDI. read_labelled_csv() drops those columns, and this is
    what lets a caller aggregate its per-row predictions back to DDI level. Rows
    come back in file order and neither reader dedupes, so the two line up
    index for index.

    None -- rather than an empty list -- distinguishes a PPI-mode CSV, which has
    no DDI level at all, from a DDI-mode split that is simply empty.
    """
    with open(path) as fh:
        reader = csv.DictReader(fh)
        cols = reader.fieldnames or []
        if "family1" not in cols or "family2" not in cols:
            return None
        pairs = []
        for row in reader:
            f1, f2 = row["family1"].strip(), row["family2"].strip()
            pairs.append((min(f1, f2), max(f1, f2)))
    return pairs


# ---------------------------------------------------------------------------
# Embeddings
# ---------------------------------------------------------------------------


def load_embeddings(path):
    """Return {protein_id: embedding_array} from an NPZ file."""
    raw = np.load(path, allow_pickle=False)
    return {k: raw[k] for k in raw.files}


# ---------------------------------------------------------------------------
# MultiQC sample naming
# ---------------------------------------------------------------------------

# The intra-dataset split order. MultiQC sorts rows alphabetically by Sample
# name, so a plain f"{label}_{split}" would show test_balanced,
# test_realistic, train, val -- the numeric prefix these helpers add is what
# forces train -> val -> test -> ... instead.
#
# A list rather than the {name: key} dict this used to be, so the order is
# stated once and the keys cannot drift from it. Every name the dict held is
# still here, including "test" and "discarded".
SPLIT_ORDER = ["train", "val", "test", "test_balanced", "test_realistic", "discarded"]

# One stderr line per unknown name per script run, not per row.
_INDEX_WARNED = set()


def _index(seq, value, what="name"):
    """Position of `value` in `seq`, or one past the end when it is unknown.

    Unknown names sort after every known one instead of displacing them, which
    is what makes an unrecognised split label or dataset a cosmetic problem and
    not a reordered report. An unknown name is reported once to stderr.
    """
    if value in seq:
        return seq.index(value)
    key = (what, value)
    if key not in _INDEX_WARNED:
        _INDEX_WARNED.add(key)
        print(f"Warning: unknown {what} {value!r} -- it will sort after every known one", file=sys.stderr)
    return len(seq)


def mqc_sample(label, split):
    """MultiQC Sample name for one (dataset, split), ordered within the dataset.

    The global, run-wide numbering lives in bin/relabel_mqc.py instead: it needs
    the run's whole dataset order, and threading that into every task would put
    the run's dataset set into the hash of tasks as expensive as SELECT_EXAMPLES
    and SAMPLE_NEGATIVES_ILP -- so adding one dataset would re-solve every other.
    relabel_mqc.py parses the name back out of the collected *_mqc.tsv files and
    rewrites it as "<NN>_<label>_<split>" just before MultiQC runs.
    """
    return f"{label}_{_index(SPLIT_ORDER, split, 'split') + 1}_{split}"


def mqc_dataset(label):
    """MultiQC Sample name for a whole dataset -- one row, no split.

    Used by the classifier-performance tables and the DDI attrition bar, which
    report one row per dataset. Identity today; it exists so those writers say
    which of the two Sample vocabularies they are using, which is what
    relabel_mqc.py keys its dataset-level vs split-level rewrite on.
    """
    return label


def mqc_category(name):
    """Bar-chart category inside one dataset's own chart -- split order only.

    No dataset component: these charts are already scoped to one dataset, so
    only the intra-dataset order is at stake and a global index would be noise.
    relabel_mqc.py leaves these alone for exactly that reason.
    """
    return f"{_index(SPLIT_ORDER, name, 'split') + 1}_{name}"
