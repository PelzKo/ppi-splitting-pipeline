#!/usr/bin/env python3
"""Check the invariants a DDI-mode run must satisfy, over a published results tree.

Not part of the pipeline -- run by hand after a run, so the end-to-end checks are
one command rather than the ad-hoc shell the implementation plan used per step.

    python bin/other/check_ddi_invariants.py --results results/results_test_ddi

Point --results at either a single dataset directory or the whole outdir, in
which case every dataset under it is checked (`_shared` is skipped). Exits 1 if
any strict dataset fails.

What is checked, per dataset:

  families   no Pfam family, and no clan, appears in more than one split
  parents    no parent protein is used by more than one split
  vocabulary every id in the *_instances.csv files is a real instance of the
             family that row declares, and its parent is in that split's own
             universe (or in the never-claimed reserve EXPAND_NEGATIVES may draw
             from when a negative pair runs short)
  cap    no (family pair, label) carries more than N examples
  labels no instance pair appears as both a positive and a negative
  universes  the three protein universes are pairwise disjoint, and disjoint
             from unclaimed.txt
  positives  when a row asked for several negative sets, every set carries the
             same positive rows -- that shared positive split is the whole point
             of the fan-out

A row asking for several negative sets (negative_sampling_method as a
comma-separated list) publishes one file set per name --
train_ilp_instances.csv, train_ilp_candidates_instances.csv, ... -- except
test_realistic, which cannot differ between them and is published once,
unsuffixed, and checked as part of every set. Each set is checked on its own and
named in the section header; a single-negative-set run keeps the historic
unsuffixed filenames and prints no name.

split_method=random is the deliberate exception: it puts the same family in
several splits so the naive baseline shows the leak, and SELECT_EXAMPLES gets
--allow-shared-parents there. Datasets whose id matches --leaky-ids are checked the
other way round -- the leak-shaped invariants must FAIL, because a random
baseline that came back clean would mean the one-split-per-parent rule was
applied where it should not have been. The vocabulary and cap checks still apply to them.
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils import instance_parent, read_instances

LABELS = ["train", "val", "test_balanced", "test_realistic"]
SUFFIX = "_instances.csv"


def split_of(label):
    """The SELECT_EXAMPLES split a sampling label draws its universe from."""
    return "test" if label.startswith("test") else label


def discover_negsets(path):
    """Return {negset: {label: csv path}} for one dataset directory.

    A single-negative-set run publishes the historic unsuffixed filenames and
    comes back under the key "". Several sets suffix every file with the set's own
    name -- except test_realistic, whose content cannot depend on the negative set,
    so it is published once, unsuffixed, and attached here to every set found. The
    label is matched longest-first, since test_balanced_ilp starts with neither
    "test" nor any label but its own.
    """
    by_negset, shared_realistic = defaultdict(dict), None
    longest_first = sorted(LABELS, key=len, reverse=True)
    for name in sorted(os.listdir(path)):
        if not name.endswith(SUFFIX):
            continue
        stem = name[: -len(SUFFIX)]
        label = next((lbl for lbl in longest_first if stem == lbl or stem.startswith(lbl + "_")), None)
        if label is None:
            continue
        negset = "" if stem == label else stem[len(label) + 1 :]
        if label == "test_realistic" and not negset:
            shared_realistic = os.path.join(path, name)
        by_negset[negset][label] = os.path.join(path, name)
    if shared_realistic is not None:
        for files in by_negset.values():
            files.setdefault("test_realistic", shared_realistic)
        # With named sets present, "" is not a set of its own -- it is only the
        # shared test_realistic file, already handed to each of them above.
        if len(by_negset) > 1 and set(by_negset.get("", {})) == {"test_realistic"}:
            del by_negset[""]
    return dict(by_negset)


def find_datasets(root):
    """Return [(dataset_id, path)], accepting either one dataset dir or a whole outdir."""
    if discover_negsets(root):
        return [(os.path.basename(os.path.abspath(root)), root)]
    found = []
    for name in sorted(os.listdir(root)):
        path = os.path.join(root, name)
        if name == "_shared" or not os.path.isdir(path):
            continue
        if discover_negsets(path):
            found.append((name, path))
    return found


def positive_pairs(rows):
    """The positive instance pairs of one labelled table, as a comparable set."""
    return {
        (row["family1"], row["family2"], tuple(sorted((row["protein1"], row["protein2"]))))
        for row in rows
        if row.get("label") == "1"
    }


def find_instances(dataset_dir, explicit):
    if explicit:
        return explicit
    candidates = [
        os.path.join(dataset_dir, "data", "instances.tsv"),
        os.path.join(os.path.dirname(os.path.abspath(dataset_dir)), "_shared", "data", "instances.tsv"),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    return None


def read_rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def read_ids(path):
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


class Report:
    """Collects check outcomes so one dataset reports all its failures, not just the first."""

    def __init__(self, dataset):
        self.dataset = dataset
        self.failures = 0

    def check(self, ok, name, detail=""):
        status = "PASS" if ok else "FAIL"
        if not ok:
            self.failures += 1
        print(f"  [{status}] {name}{(' -- ' + detail) if detail else ''}")

    def note(self, text):
        print(f"  [note] {text}")


def check_dataset(dataset, path, instances_path, target, leaky, negset="", files=None):
    """Check one negative set of one dataset; returns (failures, {label: rows})."""
    shown = f"{dataset} [{negset}]" if negset else dataset
    print(f"\n== {shown}{'  (leaky baseline: leak-shaped checks are inverted)' if leaky else ''}")
    rep = Report(dataset)
    if files is None:
        files = {lbl: os.path.join(path, f"{lbl}{SUFFIX}") for lbl in LABELS}

    rows = read_instances(instances_path)
    inst_family = {r["instance_id"]: r["family"] for r in rows}
    inst_parent = {r["instance_id"]: r["protein_id"] for r in rows}
    fam_clan = {r["family"]: r["clan"] for r in rows}

    labelled, missing = {}, []
    for lbl in LABELS:
        f = files.get(lbl)
        if f and os.path.exists(f):
            labelled[lbl] = read_rows(f)
        else:
            missing.append(lbl)
    if missing:
        rep.check(False, "all four labelled instance tables present", f"missing {', '.join(missing)}")
        return rep.failures, labelled
    rep.check(True, "all four labelled instance tables present")

    universes, unclaimed = {}, read_ids(os.path.join(path, "examples", "unclaimed.txt"))
    for split in ("train", "val", "test"):
        ids = read_ids(os.path.join(path, "examples", f"{split}_universe.txt"))
        if ids is not None:
            universes[split] = ids
    if not universes:
        rep.note("no examples/ directory published -- universe checks skipped")

    bad_family, bad_parent, over, both = [], [], [], []
    fam_by_split, clan_by_split, parents_by_split = defaultdict(set), defaultdict(set), defaultdict(set)
    all_pairs = set()

    for lbl, table in labelled.items():
        split = split_of(lbl)
        # Counted per table, not across all four: test_balanced and test_realistic
        # both expand the same test-split positives, so a shared counter would
        # report every test DDI at 2N and fail the cap check on a correct run.
        per_pair, pairs_by_label = defaultdict(int), defaultdict(set)
        for row in table:
            if not {"protein1", "protein2", "label", "family1", "family2"} <= set(row):
                rep.check(False, f"{lbl}: instance table carries the DDI columns", "expected family1/family2/label")
                return rep.failures, labelled
            for col, fam_col in (("protein1", "family1"), ("protein2", "family2")):
                iid = row[col]
                if inst_family.get(iid) != row[fam_col]:
                    bad_family.append(f"{lbl}:{iid}")
                    continue
                parent = inst_parent.get(iid, instance_parent(iid))
                parents_by_split[split].add(parent)
                allowed = universes.get(split)
                if allowed is not None and parent not in allowed and parent not in (unclaimed or set()):
                    bad_parent.append(f"{lbl}:{parent}")
                fam_by_split[split].add(row[fam_col])
                if row[fam_col] in fam_clan:
                    clan_by_split[split].add(fam_clan[row[fam_col]])
            per_pair[(row["family1"], row["family2"], row["label"])] += 1
            pairs_by_label[row["label"]].add(tuple(sorted((row["protein1"], row["protein2"]))))
        over += [f"{lbl}:{f1}-{f2}(label {y})={n}" for (f1, f2, y), n in per_pair.items() if n > target]
        both += [f"{lbl}:{p}" for p in pairs_by_label.get("1", set()) & pairs_by_label.get("0", set())]
        all_pairs |= {(f1, f2) for f1, f2, _ in per_pair}

    rep.check(not bad_family, "every id is an instance of the family its row declares", f"{len(bad_family)} bad")
    if universes:
        rep.check(
            not bad_parent,
            "every parent is in its own split's universe or the unclaimed reserve",
            f"{len(bad_parent)} outside",
        )
    rep.check(not over, f"no (family pair, label) carries more than N={target} examples", ", ".join(over[:4]))
    rep.check(not both, "no instance pair carries both labels", f"{len(both)} do")

    counts = ", ".join(f"{lbl}={len(labelled[lbl])}" for lbl in LABELS)
    rep.note(f"rows: {counts}; distinct family pairs: {len(all_pairs)}")

    # The leak-shaped ones. Under split_method=random these are supposed to fail,
    # so the assertion flips rather than being skipped -- a random baseline that
    # came back clean would mean the one-split-per-parent rule was applied where it must not be.
    shared = {}
    for a, b in (("train", "val"), ("train", "test"), ("val", "test")):
        shared[("family", a, b)] = fam_by_split[a] & fam_by_split[b]
        shared[("clan", a, b)] = clan_by_split[a] & clan_by_split[b]
        shared[("parent", a, b)] = parents_by_split[a] & parents_by_split[b]
        if universes:
            shared[("universe", a, b)] = universes.get(a, set()) & universes.get(b, set())

    for kind in ("family", "clan", "parent", "universe"):
        overlaps = {k: v for k, v in shared.items() if k[0] == kind}
        if not overlaps:
            continue
        total = sum(len(v) for v in overlaps.values())
        detail = ", ".join(f"{a}/{b}: {len(v)}" for (_, a, b), v in sorted(overlaps.items()))
        name = {
            "family": "no Pfam family in two splits",
            "clan": "no clan in two splits",
            "parent": "no parent protein in two splits",
            "universe": "protein universes pairwise disjoint",
        }[kind]
        if leaky:
            rep.check(total > 0, f"leaky baseline still leaks ({kind})", detail)
        else:
            rep.check(total == 0, name, detail)

    if unclaimed is not None and universes and not leaky:
        bleed = {s: len(u & unclaimed) for s, u in universes.items()}
        rep.check(
            not any(bleed.values()),
            "unclaimed.txt disjoint from every universe",
            ", ".join(f"{s}: {n}" for s, n in sorted(bleed.items())),
        )
        rep.note(f"{len(unclaimed)} proteins never claimed by any split")

    return rep.failures, labelled


def check_shared_positives(dataset, tables):
    """The point of asking for several negative sets is that they sit on one shared
    positive split, so every set must carry the same positive rows."""
    if len(tables) < 2:
        return 0
    print(f"\n== {dataset}  (across {len(tables)} negative sets)")
    rep = Report(dataset)
    ref_negset = sorted(tables)[0]
    for lbl in LABELS:
        ref = positive_pairs(tables[ref_negset].get(lbl, []))
        differing = []
        for ns, rows in sorted(tables.items()):
            if ns == ref_negset:
                continue
            got = positive_pairs(rows.get(lbl, []))
            if got != ref:
                differing.append(f"{ns}: {len(got ^ ref)} differ")
        rep.check(
            not differing, f"{lbl}: same positives in every negative set ({len(ref)} pairs)", "; ".join(differing)
        )
    return rep.failures


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, help="a dataset directory, or the whole outdir")
    ap.add_argument("--instances", help="instances.tsv; found under the results tree when omitted")
    ap.add_argument("--examples-target", type=int, default=5, help="N, the per-DDI example cap (default 5)")
    ap.add_argument(
        "--leaky-ids",
        default=r"(?i)leak|random",
        help="regex for dataset ids run with split_method=random, whose leak-shaped invariants are inverted",
    )
    args = ap.parse_args()

    datasets = find_datasets(args.results)
    if not datasets:
        sys.exit(f"No dataset with *_instances.csv found under {args.results} -- is this a DDI-mode run?")

    leaky_re = re.compile(args.leaky_ids)
    failures = 0
    for dataset, path in datasets:
        instances_path = find_instances(path, args.instances)
        if not instances_path:
            print(f"\n== {dataset}\n  [FAIL] instances.tsv not found -- pass --instances")
            failures += 1
            continue
        leaky = bool(leaky_re.search(dataset))
        negsets = discover_negsets(path)
        tables = {}
        for negset in sorted(negsets):
            n_failed, labelled = check_dataset(
                dataset, path, instances_path, args.examples_target, leaky, negset, negsets[negset]
            )
            failures += n_failed
            tables[negset] = labelled
        failures += check_shared_positives(dataset, tables)

    print(f"\n{len(datasets)} dataset(s) checked, {failures} failed check(s)")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
