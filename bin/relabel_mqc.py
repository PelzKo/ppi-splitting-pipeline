#!/usr/bin/env python3
"""
Apply the run-wide MultiQC display order, just before MultiQC runs.

Every *_mqc.tsv / *_mqc.html this pipeline emits is written by a task that sees
one dataset. The global numbering cannot be: it depends on how many datasets the
run has and in which order the caller wants them. Threading that into every
emitting task would put the run's whole dataset set into the hash of tasks as
expensive as SELECT_EXAMPLES and SAMPLE_NEGATIVES_ILP, so adding one dataset
would re-solve every other one on -resume. This script does it in one place
instead: it reads the collected files, rewrites the two things MultiQC sorts on,
and writes a config that pins the section order.

What it rewrites
----------------
Sample column (*_mqc.tsv body rows only)
    split-level    "<label>_<k>_<split>[_<tail>]"  ->  "<NN>_<label>_<split>[_<tail>]"
    dataset-level  "<label>"                       ->  "<NN>_<label>"
    Split-level NN indexes the (dataset x split) cross product, dataset-level NN
    indexes the datasets alone, each zero-padded to the width of the largest
    index the run can produce -- so MultiQC's alphabetical sort reproduces the
    intended order and the padding does not change halfway through a report.
    utils.mqc_category() values ("1_train") carry no dataset component on
    purpose -- they are categories inside one dataset's own chart -- and are left
    alone.

Section id (the "# id:" / "id:" header line)
    "<kind>_<label>"  ->  "<kind>_<NN>_<label>"   for the per-dataset kinds
    Shared sections (pos_neg_bar, ddi_attrition, classifier_metrics_*, ...) have
    no dataset component and keep their ids verbatim, so external references to
    them survive. They are ordered by the generated config alone.

What it writes
--------------
<out>/...                 the rewritten tree, relative paths preserved
<config>                  report_section_order: shared summaries first, then the
                          per-dataset kinds one at a time with the datasets in
                          the resolved order. MultiQC renders higher `order`
                          first; General Statistics is always first regardless.
"""

import argparse
import os
import re
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from utils import SPLIT_ORDER, _index  # noqa: E402

# Section-id prefixes that are followed by a dataset display label. Longest
# match wins, so ddi_examples_bar is not read as ddi_examples.
DATASET_KINDS = [
    "split_bar",
    "ddi_examples_bar",
    "ddi_examples_stats",
    "ddi_examples_ilp",
    "bias_scatter",
    "similarity_heatmap",
]

# Report order. The shared, cross-dataset summaries come first -- they are the
# ones a reader compares datasets in -- then each per-dataset kind in turn with
# the datasets inside it in the resolved order. Partitioning deliberately follows
# the attrition waterfall, which is the summary of it.
#
# A "kind" here is a section id with its "_<NN>_<label>" tail removed, so both
# shared and per-dataset ids classify the same way. An id matching none of these
# is placed after all of them, sorted, rather than dropped.
KIND_ORDER = [
    # shared, one section merging every dataset's rows
    "neg_generalstats",
    "pos_neg_bar",
    "classifier_metrics",
    "ddi_expand_bar",
    "ddi_expand_stats",
    "neg_ilp_diagnostics",
    "neg_ilp_residuals",
    "ddi_attrition",
    # per-dataset, one section each
    "split_bar",
    "ddi_examples_bar",
    "ddi_examples_stats",
    "ddi_examples_ilp",
    "bias_scatter",
    "similarity_heatmap",
]

# Longest first: "test" would otherwise swallow "test_balanced".
_SPLITS_LONGEST_FIRST = sorted(SPLIT_ORDER, key=len, reverse=True)

# utils.mqc_category() output: a category inside one dataset's own chart, which
# carries no dataset component by design. Recognised so it is not reported as an
# unmatched Sample name.
_CATEGORY = re.compile(r"^\d+_(?:" + "|".join(re.escape(s) for s in SPLIT_ORDER) + r")$")

_TSV_SECTION_ID = re.compile(r"^# id:\s*(['\"])(?P<id>.*?)\1\s*$")
_HTML_SECTION_ID = re.compile(r"^id:\s*(['\"])(?P<id>.*?)\1\s*$")


class Renamer:
    """Holds the resolved order and every rewrite derived from it."""

    def __init__(self, order):
        self.order = order
        self.n_datasets = max(len(order), 1)
        self.dataset_width = len(str(self.n_datasets))
        self.split_width = len(str(self.n_datasets * len(SPLIT_ORDER)))
        # Longest first so "minimal_leakage" cannot shadow "minimal_leakage_hcni".
        self.labels_longest_first = sorted(order, key=len, reverse=True)
        self.unmatched = set()
        self.section_ids = []

    # -- indices ----------------------------------------------------------
    def dataset_index(self, label):
        return f"{_index(self.order, label, 'dataset') + 1:0{self.dataset_width}d}"

    def split_index(self, label, split):
        idx = _index(self.order, label, "dataset") * len(SPLIT_ORDER) + _index(SPLIT_ORDER, split, "split") + 1
        return f"{idx:0{self.split_width}d}"

    # -- Sample column ----------------------------------------------------
    def _split_label(self, sample):
        """Return (label, rest) if `sample` starts with a known display label."""
        for label in self.labels_longest_first:
            if sample == label:
                return label, None
            if sample.startswith(label + "_"):
                return label, sample[len(label) + 1 :]
        return None, None

    def sample(self, value):
        """Rewrite one Sample cell. Anything unrecognised is returned unchanged."""
        label, rest = self._split_label(value)
        if label is None:
            if not _CATEGORY.match(value):
                self.unmatched.add(value)
            return value
        if rest is None:
            # Dataset-level: utils.mqc_dataset(), one row per dataset.
            return f"{self.dataset_index(label)}_{label}"
        # Split-level: utils.mqc_sample() wrote "<label>_<k>_<split>", possibly
        # with a further "_<tail>" (sample_negatives_ilp.py's per-protein rows).
        m = re.match(r"^(\d+)_(?P<tail>.+)$", rest)
        if not m:
            self.unmatched.add(value)
            return value
        tail = m.group("tail")
        for split in _SPLITS_LONGEST_FIRST:
            if tail == split:
                return f"{self.split_index(label, split)}_{label}_{split}"
            if tail.startswith(split + "_"):
                return f"{self.split_index(label, split)}_{label}_{tail}"
        self.unmatched.add(value)
        return value

    # -- section id -------------------------------------------------------
    def section_id(self, sid):
        """Insert the dataset index into a per-dataset section id."""
        kinds = sorted((k for k in DATASET_KINDS if sid.startswith(k + "_")), key=len, reverse=True)
        for kind in kinds:
            label = sid[len(kind) + 1 :]
            if label in self.order:
                return f"{kind}_{self.dataset_index(label)}_{label}"
            # A per-dataset kind whose label is not one of the run's -- report it
            # rather than inventing an index for it.
            self.unmatched.add(sid)
            return sid
        return sid


def kind_of(sid, order, dataset_width):
    """The section id with its "_<NN>_<label>" tail removed, for KIND_ORDER."""
    for kind in sorted(DATASET_KINDS, key=len, reverse=True):
        prefix = kind + "_"
        if sid.startswith(prefix):
            rest = sid[len(prefix) :]
            if re.match(rf"^\d{{{dataset_width}}}_", rest):
                return kind
    for kind in sorted(KIND_ORDER, key=len, reverse=True):
        if sid == kind or sid.startswith(kind + "_"):
            return kind
    return None


def rewrite_tsv(lines, renamer):
    """Rewrite a custom-content TSV: the section id, then the Sample column.

    MultiQC's TSV custom content is a "#"-prefixed YAML header, then one column
    header row, then the body -- so the first non-"#" line is skipped and every
    line after it has its first field rewritten.
    """
    out = []
    seen_header_row = False
    for line in lines:
        if line.startswith("#"):
            m = _TSV_SECTION_ID.match(line.rstrip("\n"))
            if m:
                sid = m.group("id")
                new = renamer.section_id(sid)
                renamer.section_ids.append(new)
                out.append(f"# id: '{new}'\n")
                continue
            out.append(line)
            continue
        if not seen_header_row:
            # The column-header row. Its first cell is the literal "Sample".
            seen_header_row = True
            out.append(line)
            continue
        if not line.strip():
            out.append(line)
            continue
        fields = line.rstrip("\n").split("\t")
        fields[0] = renamer.sample(fields[0])
        out.append("\t".join(fields) + "\n")
    return out


def rewrite_html(lines, renamer):
    """Rewrite the leading HTML-comment header of a custom-content HTML file."""
    out = []
    in_header = False
    done = False
    for line in lines:
        stripped = line.strip()
        if not done:
            if stripped.startswith("<!--"):
                in_header = True
            if in_header:
                m = _HTML_SECTION_ID.match(stripped)
                if m:
                    sid = m.group("id")
                    new = renamer.section_id(sid)
                    renamer.section_ids.append(new)
                    out.append(f"id: '{new}'\n")
                    continue
                if stripped.endswith("-->"):
                    in_header = False
                    done = True
        out.append(line)
    return out


def write_config(path, section_ids, order, dataset_width):
    """Pin the section order. Higher `order` renders first in MultiQC.

    Every id present is listed, so nothing falls back to MultiQC's discovery
    order -- which is what made the section order arbitrary per run, since
    MULTIQC stages its inputs into directories numbered by channel arrival.
    """
    ids = sorted(set(section_ids))
    by_kind = {}
    unknown = []
    for sid in ids:
        kind = kind_of(sid, order, dataset_width)
        if kind is None:
            unknown.append(sid)
        else:
            by_kind.setdefault(kind, []).append(sid)

    ranked = []
    for kind in KIND_ORDER:
        # Within a per-dataset kind the ids already carry the zero-padded
        # dataset index, so a plain sort is the resolved order.
        ranked.extend(sorted(by_kind.get(kind, [])))
    ranked.extend(sorted(unknown))

    # MultiQC numbers sections from 10 (bottom) upwards in steps of 10, and a
    # higher number renders nearer the top. Every discovered section is listed, so
    # none keeps a default -- but the base is lifted well clear of that range
    # anyway, so a section a future writer adds and this script does not recognise
    # cannot tie with one it placed deliberately.
    #
    # general_stats is pinned above everything. MultiQC renders the General
    # Statistics table first on its own, so the entry is belt-and-braces; an
    # unrecognised key here is ignored rather than an error.
    step = 10
    base = 100000
    with open(path, "w") as fh:
        fh.write("# Generated per run by bin/relabel_mqc.py -- do not edit.\n")
        fh.write("# Higher order renders nearer the top. Shared cross-dataset summaries\n")
        fh.write("# first, then each per-dataset kind with the datasets in the resolved\n")
        fh.write("# --mqc-order.\n")
        fh.write("report_section_order:\n")
        fh.write(f"  general_stats:\n    order: {base + step}\n")
        for i, sid in enumerate(ranked):
            fh.write(f"  {sid}:\n    order: {base - i * step}\n")
    return ranked


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source", help="directory holding the collected *_mqc.tsv / *_mqc.html files")
    ap.add_argument("--out", required=True, help="directory to write the rewritten tree into")
    ap.add_argument(
        "--mqc-order",
        default="",
        help="comma-joined display labels, in report order (main.nf's resolved order)",
    )
    ap.add_argument("--config", required=True, help="MultiQC config to write (report_section_order)")
    args = ap.parse_args()

    order = [s.strip() for s in args.mqc_order.split(",") if s.strip()]
    renamer = Renamer(order)

    out_root = os.path.abspath(args.out)
    src_root = os.path.abspath(args.source)
    os.makedirs(out_root, exist_ok=True)

    n_files = 0
    for dirpath, dirnames, filenames in os.walk(src_root):
        # Never descend into our own output, and skip Nextflow's task files.
        dirnames[:] = [d for d in dirnames if os.path.join(dirpath, d) != out_root and not d.startswith(".")]
        for name in sorted(filenames):
            if not (name.endswith("_mqc.tsv") or name.endswith("_mqc.html")):
                continue
            src = os.path.join(dirpath, name)
            rel = os.path.relpath(src, src_root)
            dst = os.path.join(out_root, rel)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            with open(src) as fh:
                lines = fh.readlines()
            if name.endswith("_mqc.tsv"):
                lines = rewrite_tsv(lines, renamer)
            else:
                lines = rewrite_html(lines, renamer)
            with open(dst, "w") as fh:
                fh.writelines(lines)
            n_files += 1

    if n_files == 0:
        print(f"Warning: no *_mqc.tsv / *_mqc.html files found under {src_root}", file=sys.stderr)

    ranked = write_config(args.config, renamer.section_ids, order, renamer.dataset_width)

    print(
        f"Relabelled {n_files} MultiQC file(s), {len(ranked)} section(s), "
        f"{len(order)} dataset(s): {', '.join(order) if order else '(none)'}",
        file=sys.stderr,
    )
    if renamer.unmatched:
        # Not fatal: an unrecognised Sample or section id keeps its old name and
        # sorts after everything else, which is a cosmetic problem only.
        shown = sorted(renamer.unmatched)[:10]
        more = "" if len(renamer.unmatched) <= 10 else f" (and {len(renamer.unmatched) - 10} more)"
        print(
            "Warning: left unchanged because they match no display label of this run: " f"{', '.join(shown)}{more}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
