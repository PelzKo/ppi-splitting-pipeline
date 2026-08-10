#!/usr/bin/env python3
"""Shared I/O utilities for PPI pipeline scripts."""

import csv
import sys

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

# MultiQC sorts rows alphabetically by Sample name, so a plain f"{id}_{split}"
# would show test_balanced, test_realistic, train, val. This numeric prefix
# forces train/val/test/test_balanced/test_realistic order per dataset.
SPLIT_SORT_KEY = {"train": 1, "val": 2, "test": 3, "test_balanced": 4, "test_realistic": 5, "discarded": 6}


def mqc_sample(id_, split):
    """Build a MultiQC Sample name for one dataset+split that sorts correctly."""
    return f"{id_}_{SPLIT_SORT_KEY.get(split, 9)}_{split}"


def mqc_category(name):
    """Same sort-key prefix as mqc_sample, for bar-chart categories that are
    already scoped to one dataset's own chart (no id_ qualification needed)."""
    return f"{SPLIT_SORT_KEY.get(name, 9)}_{name}"
