#!/usr/bin/env python3
"""
DDI mode: one attrition waterfall per dataset, from the splitting stage's own bars.

Every input DDI ends up in exactly one of four buckets, and each is already
counted by a MultiQC bar written upstream -- so this reduces those bars rather
than re-deriving anything, which is what stops the waterfall and the per-stage
charts from disagreeing:

    discarded (cross-cluster)   the partitioner put the DDI's two families in
                                different clusters   [sort_ppis.py / remove_redundant.py]
    discarded (CD-HIT-2D)       one family lost an instance to redundancy
                                filtering            [remove_redundant.py]
    dropped (no example)        the DDI survived splitting but reached zero
                                domain-instance examples -- every candidate pair
                                blocked because Barrier B gave its parent
                                proteins to another split, or the family had no
                                usable instance      [select_examples.py]
    kept                        made it into a split with >= 1 example

The two source bars are identified by their column headers rather than their
filenames, so this keeps working if a splitter's output file is renamed:

    Sample  Kept  Discarded (KaHIP/ILP)  Discarded (CD-HIT-2D)   <- partitioning
    Sample  Full N examples  Partial  Dropped (0 examples)       <- example selection

Input
-----
Every *_mqc.tsv one dataset's splitting stage produced (extras are ignored).

Output
------
ddi_attrition_mqc.tsv  -- one bargraph row for this dataset; MultiQC merges the
                          section run-wide, so the rows line up as one chart.
"""

import argparse
import sys

_PARTITION_MARKER = "Discarded (KaHIP/ILP)"
_SELECTION_MARKER = "Dropped (0 examples)"


def read_bar(path):
    """Return (header, [(sample, [int, ...]), ...]) from a bargraph MultiQC TSV.

    None if the file has no tab-separated header row -- e.g. an HTML section or
    one of the stats tables that share this channel.
    """
    header, rows = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            try:
                rows.append((parts[0], [int(float(v)) for v in parts[1:]]))
            except ValueError:
                # A non-numeric body column means this is not a bargraph.
                return None
    if header is None or len(header) < 2:
        return None
    return header, rows


def find_bar(paths, marker):
    """The first parsed bar whose header contains `marker`, or None."""
    for path in paths:
        bar = read_bar(path)
        if bar and marker in bar[0]:
            return bar
    return None


def column(header, rows, name):
    """Sum one named column over every row; 0 if the column is absent."""
    if name not in header:
        return 0
    i = header.index(name) - 1
    return sum(vals[i] for _, vals in rows if i < len(vals))


def split_rows(rows):
    """The per-split rows of the partitioning bar, i.e. everything but its
    single stacked "discarded" row. mqc_category() prefixes a numeric sort key,
    so match on the suffix rather than on equality."""
    return [(s, v) for s, v in rows if not s.endswith("discarded")]


def collect(paths):
    """Return the four bucket counts plus what was and wasn't found."""
    partition = find_bar(paths, _PARTITION_MARKER)
    selection = find_bar(paths, _SELECTION_MARKER)

    if partition is None:
        sys.exit(
            f"none of {len(paths)} input file(s) is the PPI/DDI Partitioning bar "
            f"(no '{_PARTITION_MARKER}' column) -- cannot account for the input DDIs"
        )

    header, rows = partition
    after_split = column(header, split_rows(rows), "Kept")
    cross = column(header, rows, _PARTITION_MARKER)
    cdhit = column(header, rows, "Discarded (CD-HIT-2D)")

    if selection is None:
        # split_method=random reaches SELECT_EXAMPLES too, so this should not
        # happen in DDI mode -- report it rather than silently showing 0 drops.
        print(
            f"Warning: no example-selection bar (no '{_SELECTION_MARKER}' column) among the "
            "inputs -- reporting 0 DDIs dropped for want of an example",
            file=sys.stderr,
        )
        dropped = 0
    else:
        dropped = column(selection[0], selection[1], _SELECTION_MARKER)

    return {
        "kept": after_split - dropped,
        "cross": cross,
        "cdhit": cdhit,
        "dropped": dropped,
    }, (after_split, selection is not None)


def write_mqc(stats, id_):
    with open("ddi_attrition_mqc.tsv", "w") as fh:
        fh.write(
            "# id: 'ddi_attrition'\n"
            "# section_name: 'DDI Attrition'\n"
            "# description: 'Where every input DDI went, one bar per dataset. Discarded "
            "cross-cluster means the partitioner put the two Pfam families in different "
            "clusters; discarded by CD-HIT-2D means a family lost an instance to redundancy "
            "filtering; dropped (no example) means the DDI survived splitting but reached zero "
            "domain-instance examples, because Barrier B gave its parent proteins to another "
            "split or the family had no usable instance. The four series sum to the input DDI "
            "count.'\n"
            "# plot_type: 'bargraph'\n"
            "# pconfig:\n"
            "#     id: 'ddi_attrition_plot'\n"
            "#     title: 'DDI Attrition: where the input DDIs went'\n"
            "#     ylab: '# DDIs'\n"
            "Sample\tKept\tDiscarded (cross-cluster)\tDiscarded (CD-HIT-2D)\tDropped (no example)\n"
            f"{id_}\t{stats['kept']}\t{stats['cross']}\t{stats['cdhit']}\t{stats['dropped']}\n"
        )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsvs", nargs="+", help="the splitting stage's *_mqc.tsv files for one dataset")
    ap.add_argument("--id", required=True, help="Dataset ID, for MultiQC tagging")
    args = ap.parse_args()

    stats, (after_split, had_selection) = collect(args.tsvs)
    total = stats["kept"] + stats["cross"] + stats["cdhit"] + stats["dropped"]

    print(
        f"{total} input DDIs: {stats['cross']} cross-cluster, {stats['cdhit']} CD-HIT-2D, "
        f"{stats['dropped']} no example, {stats['kept']} kept",
        file=sys.stderr,
    )
    if stats["kept"] < 0:
        # More DDIs dropped for want of an example than survived splitting: the
        # two bars are counting different populations, so the waterfall is wrong.
        sys.exit(
            f"example selection dropped {stats['dropped']} DDIs but only {after_split} survived "
            "splitting -- the partitioning and selection bars disagree"
        )
    if not had_selection:
        print("  (example-selection stage not accounted for -- see warning above)", file=sys.stderr)

    write_mqc(stats, args.id)


if __name__ == "__main__":
    main()
