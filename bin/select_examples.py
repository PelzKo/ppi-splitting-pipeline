#!/usr/bin/env python3
"""
DDI mode: turn family-pair DDIs into domain-instance-pair examples (Barrier B).

The split CSVs handed over by SOLVE_ILP/SORT_PPIS/SPLIT_RANDOM (and filtered by
REMOVE_REDUNDANT) hold Pfam *family* pairs, but a row a classifier can train on
is a pair of concrete domain *instances* -- and the parent proteins of those two
instances must not turn up in another split, or the parent's other domains carry
the interaction across the split boundary. That is Barrier B.

Selection and Barrier B are one problem, not two: claiming a parent protein for
one split removes it from every other split's options, so which examples a split
can still reach depends on what the other splits took. Greedy order would decide
the outcome, so it is solved as an ILP over all three splits jointly.

Two reductions keep that ILP small:

  * Only *contested* parents -- proteins carrying candidates in more than one
    split -- get claim variables. A parent in play for a single split can be
    claimed for free, so every DDI whose candidates touch no contested parent is
    decoupled from the rest of the problem and is settled by a local diversity
    greedy instead. This is a decomposition, not an approximation.
  * A per-DDI shortlist caps the candidate pool at K = shortlist_factor * N. At
    the default pool size (M = N = 5, so 25 candidates) it trims almost nothing
    and exists as a guard for a larger --ddi_examples_pool_factor.

The objective is lexicographic, solved as one bounded stage per level, so the
constants stay at the scale of the example counts instead of the products a
single weighted objective would need:

  1. keep as many positive DDIs as possible at >= 1 example (a DDI that reaches
     zero is dropped outright, so this level is what stops the solver starving
     one DDI to complete another),
  2. then maximise the positive DDIs reaching the full N, then total positive
     examples, then example diversity (distinct parents),
  3. then the same for the candidate_network negatives, which ride in the same
     claim accounting so their examples inherit Barrier B by construction --
     but strictly below the positives, so a negative can never take a parent a
     positive needed.

--no-barrier-b turns the whole claim mechanism off, which is what split_method=
random needs: that path deliberately puts a node in more than one split so the
naive baseline shows the leak, and enforcing Barrier B on top would repair part
of it. Every unit is then decoupled and the ILP is skipped entirely.

Note the caveat that follows from step 3: a candidate pair SAMPLE_NEGATIVES does
not ultimately select will have claimed proteins for nothing, and purely random
negatives stay outside this ILP altogether -- so Barrier B is exact for
high-confidence negatives and holds via the per-split protein universe for the
rest.
"""

import argparse
import math
import os
import random
import sys
from collections import defaultdict
from dataclasses import dataclass, field

import cvxpy as cp
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import (  # noqa: E402
    diverse_pick,
    instances_by_family,
    mqc_sample,
    pair_candidates,
    read_fasta,
    read_instances,
    read_ppis,
    write_ppi_csv,
)
from solve_ilp import _seed_option  # noqa: E402

SPLITS = ["train", "val", "test"]
EXAMPLE_COLUMNS = ["protein1", "protein2", "family1", "family2"]


@dataclass
class Unit:
    """One DDI (positive) or candidate_network pair, with its example candidates."""

    split: str
    kind: str  # "pos" | "cand"
    row: dict  # the source CSV row; protein1/protein2 hold the two families
    fam1: str
    fam2: str
    cands: list  # [(instance_a, instance_b)] after shortlisting
    n_raw: int  # candidates before the shortlist trim
    picked: list = field(default_factory=list)


# ---------------------------------------------------------------------------
# Candidate enumeration
# ---------------------------------------------------------------------------
#
# pair_candidates() and diverse_pick() live in utils.py: EXPAND_NEGATIVES draws
# the sampled negatives' instance pairs by the same two rules, and it must not
# import this module (which would pull cvxpy into a process that needs no solver).


def build_units(kind, rows_by_split, members, split_ids, parent_of, shortlist_k, seed):
    """Wrap each interaction row as a Unit, enumerating and shortlisting candidates."""
    units = []
    for split in SPLITS:
        for row in rows_by_split[split]:
            fam1, fam2 = row["protein1"], row["protein2"]
            cands = pair_candidates(fam1, fam2, members, split_ids[split])
            n_raw = len(cands)
            if n_raw > shortlist_k:
                # Seeded per unit, not from one global RNG, so the shortlist does
                # not depend on how many units were processed before this one.
                rng = random.Random(f"{seed}:{kind}:{split}:{fam1}:{fam2}")
                cands = diverse_pick(cands, shortlist_k, parent_of, rng)
            units.append(Unit(split, kind, row, fam1, fam2, cands, n_raw))
    return units


# ---------------------------------------------------------------------------
# The selection ILP
# ---------------------------------------------------------------------------


def _sum(var, idx):
    """cp.sum over an index list, tolerating the empty list."""
    return cp.sum(var[idx]) if idx else cp.Constant(0.0)


def _run_stage(problem, label, solver, seed, secs, verbose):
    """Solve one lexicographic stage; return its objective value or None."""
    kwargs = dict(time_limit=secs, verbose=verbose)
    if solver:
        seed_opt = _seed_option(solver, seed)
        if not seed_opt:
            print(
                f"Warning: no seed option is known for solver {solver}, so its tie-breaking is "
                f"unseeded and this selection is not reproducible run to run.",
                file=sys.stderr,
            )
        problem.solve(solver=solver, **kwargs, **seed_opt)
    else:
        print(
            "Warning: no --solver given, so CVXPY picks one and its internal randomisation stays "
            "unseeded; pass --solver for a reproducible selection.",
            file=sys.stderr,
        )
        problem.solve(**kwargs)

    if problem.status not in cp.settings.SOLUTION_PRESENT:
        print(f"  stage '{label}': solver status {problem.status}", file=sys.stderr)
        return None
    if problem.status == cp.settings.USER_LIMIT:
        print(
            f"  stage '{label}': hit the {secs}s limit before proving optimality; "
            f"using the best incumbent found (suboptimal).",
            file=sys.stderr,
        )
    print(f"  stage '{label}': objective {problem.value:,.4g} ({problem.status})", file=sys.stderr)
    return float(problem.value)


def solve_selection(units, parent_of, contested, n, lam, max_sec, solver, seed, verbose):
    """Choose <= n examples per unit subject to Barrier B. Fills Unit.picked.

    Variables
        y[e]   in {0,1}  candidate example e is selected
        nz[d]  in {0,1}  unit d ends up with >= 1 example
        r[d]   in {0,1}  unit d reaches the full n
        c[p,s] in {0,1}  contested parent p is claimed by split s
        o[d,p] >= 0      how often unit d reuses parent p beyond the first time

    Constraints
        (1) sum_{e in d} y[e] <= min(n, |C_d|)          per-unit cap (R5: n is a cap)
        (2) sum_{e in d} y[e] >= nz[d]
        (3) sum_{e in d} y[e] >= n * r[d]
        (4) sum_s c[p,s] <= 1                           Barrier B
        (5) y[e] <= c[parent(e), split(d(e))]           both parents, contested only
        (6) o[d,p] >= (uses of p by d's selected examples) - 1
    """
    cand_index, unit_rows, unit_cols, caps = [], [], [], []
    for ui, u in enumerate(units):
        for pair in u.cands:
            unit_rows.append(ui)
            unit_cols.append(len(cand_index))
            cand_index.append((ui, pair))
        caps.append(min(n, len(u.cands)))
    n_units, n_cand = len(units), len(cand_index)
    if n_cand == 0:
        return

    # (p, split) claim variables exist only for contested parents -- an
    # uncontested parent is in play for one split only, so claiming it costs
    # nothing and constraint 5 would be slack.
    claim_index = {}
    for ui, u in enumerate(units):
        for a, b in u.cands:
            for p in (parent_of[a], parent_of[b]):
                if p in contested:
                    claim_index.setdefault((p, u.split), len(claim_index))

    y = cp.Variable(n_cand, boolean=True)
    nz = cp.Variable(n_units, boolean=True)
    r = cp.Variable(n_units, boolean=True)

    per_unit = sp.coo_matrix((np.ones(n_cand), (unit_rows, unit_cols)), shape=(n_units, n_cand)).tocsr()
    unit_examples = per_unit @ y
    cons = [unit_examples <= np.array(caps, dtype=float), unit_examples >= nz, unit_examples >= n * r]

    if claim_index:
        c = cp.Variable(len(claim_index), boolean=True)
        by_prot = defaultdict(list)
        for (p, _), idx in claim_index.items():
            by_prot[p].append(idx)
        rows, cols = [], []
        for pi, idxs in enumerate(sorted(by_prot.values())):
            for idx in idxs:
                rows.append(pi)
                cols.append(idx)
        barrier = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(by_prot), len(claim_index))).tocsr()
        cons.append(barrier @ c <= np.ones(len(by_prot)))

        e_idx, cl_idx = [], []
        for k, (ui, (a, b)) in enumerate(cand_index):
            for p in {parent_of[a], parent_of[b]}:
                if p in contested:
                    e_idx.append(k)
                    cl_idx.append(claim_index[(p, units[ui].split)])
        if e_idx:
            cons.append(y[e_idx] <= c[cl_idx])

    # Diversity: one overflow variable per (unit, parent). A pair whose two
    # instances share a parent contributes 2 to that parent's count -- coo_matrix
    # sums the duplicate entries -- which is exactly the double spend it is.
    over_index, o_rows, o_cols = {}, [], []
    for k, (ui, (a, b)) in enumerate(cand_index):
        for p in (parent_of[a], parent_of[b]):
            key = (ui, p)
            if key not in over_index:
                over_index[key] = len(over_index)
            o_rows.append(over_index[key])
            o_cols.append(k)
    o = cp.Variable(len(over_index), nonneg=True)
    per_parent = sp.coo_matrix((np.ones(len(o_rows)), (o_rows, o_cols)), shape=(len(over_index), n_cand)).tocsr()
    cons.append(per_parent @ y - 1.0 <= o)

    pos_u = [i for i, u in enumerate(units) if u.kind == "pos"]
    cand_u = [i for i, u in enumerate(units) if u.kind == "cand"]
    pos_y = [k for k, (ui, _) in enumerate(cand_index) if units[ui].kind == "pos"]
    cand_y = [k for k, (ui, _) in enumerate(cand_index) if units[ui].kind == "cand"]
    pos_o = [idx for (ui, _), idx in over_index.items() if units[ui].kind == "pos"]
    cand_o = [idx for (ui, _), idx in over_index.items() if units[ui].kind == "cand"]

    # Within a stage, `big` only has to outrank that stage's own example and
    # diversity tail, so it stays at the scale of the example count rather than
    # the product a single all-levels objective would need.
    stages = []
    if pos_u:
        stages.append(("keep every positive DDI", _sum(nz, pos_u), 0.2))
        big = float(sum(caps[i] for i in pos_u)) + 1.0
        stages.append(("fill positives to N", big * _sum(r, pos_u) + _sum(y, pos_y) - lam * _sum(o, pos_o), 0.5))
    if cand_u:
        stages.append(("expand candidate negatives", _sum(y, cand_y) - lam * _sum(o, cand_o), 0.3))

    total_share = sum(s for _, _, s in stages)
    for label, expr, share in stages:
        secs = max(10, int(max_sec * share / total_share))
        value = _run_stage(cp.Problem(cp.Maximize(expr), cons), label, solver, seed, secs, verbose)
        if value is None:
            print(
                f"SELECT_EXAMPLES could not solve the '{label}' stage. The model is always feasible "
                f"(selecting nothing satisfies every constraint), so this is a solver failure rather "
                f"than an over-constrained problem -- check the solver's own log above.",
                file=sys.stderr,
            )
            sys.exit(1)
        # Freeze this level before optimising the next. An incumbent from a
        # timed-out stage is still achievable, so the bound stays satisfiable.
        cons = cons + [expr >= value - 1e-4]

    y_val = np.asarray(y.value).ravel()
    for k, (ui, pair) in enumerate(cand_index):
        if y_val[k] > 0.5:
            units[ui].picked.append(pair)
    for ui, u in enumerate(units):
        # Trim defensively: a solver returning 0.5+eps on more candidates than
        # the cap would otherwise leak an extra example into the output.
        u.picked = sorted(u.picked)[: caps[ui]]


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def example_rows(units):
    """Instance-pair rows for one split, carrying the family pair and provenance.

    protein1/protein2 hold the two instance ids -- the same column names every
    downstream consumer already reads -- with family1/family2 added so
    train_classifier.py can aggregate example predictions back to the DDI, and
    any extra column the input CSV carried (source, confidence, ...) copied
    through untouched.
    """
    rows = []
    for u in units:
        extra = {k: v for k, v in u.row.items() if k not in ("protein1", "protein2")}
        for a, b in u.picked:
            rows.append({"protein1": a, "protein2": b, "family1": u.fam1, "family2": u.fam2, **extra})
    return rows


def write_examples(units, path):
    """write_ppi_csv, but with the DDI header for an empty split.

    write_ppi_csv falls back to protein1,protein2 when there are no rows, which
    would leave a zero-fraction split's file without family1/family2 and make it
    read as a PPI-mode file to the next stage.
    """
    rows = example_rows(units)
    if rows:
        write_ppi_csv(rows, path)
        return len(rows)
    with open(path, "w", newline="") as fh:
        fh.write(",".join(EXAMPLE_COLUMNS) + "\n")
    return 0


def write_ids(ids, path):
    with open(path, "w") as fh:
        for i in sorted(ids):
            fh.write(f"{i}\n")


def write_mqc(stats, id_):
    """Two MultiQC sections: the per-split DDI outcome bar and a stats table."""
    with open("select_examples_bar_mqc.tsv", "w") as fh:
        fh.write(
            f"# id: 'ddi_examples_bar_{id_}'\n"
            f"# section_name: 'DDI Example Selection: {id_}'\n"
            "# description: 'Positive DDIs per split, coloured by how many domain-instance "
            "examples each one ended up with. A DDI that reached zero examples -- every "
            "candidate blocked because Barrier B gave its parent proteins to another split, "
            "or the family had no usable instance -- is dropped from the split.'\n"
            "# plot_type: 'bargraph'\n"
            "# pconfig:\n"
            f"#     id: 'ddi_examples_bar_plot_{id_}'\n"
            f"#     title: 'DDI Example Selection: DDIs per outcome ({id_})'\n"
            "#     ylab: '# DDIs'\n"
            "Sample\tFull N examples\tPartial\tDropped (0 examples)\n"
        )
        for split in SPLITS:
            st = stats[split]
            fh.write(f"{mqc_sample(id_, split)}\t{st['full']}\t{st['partial']}\t{st['dropped']}\n")

    with open("select_examples_stats_mqc.tsv", "w") as fh:
        fh.write(
            f"# id: 'ddi_examples_stats_{id_}'\n"
            f"# section_name: 'DDI Example Selection Stats: {id_}'\n"
            "# description: 'Per split: the DDIs that survived selection, the instance-pair "
            "examples written for them, the parent proteins claimed for this split (its protein "
            "universe, which the negative sampler draws from), and how many DDIs had their "
            "candidate pool shortlisted before the ILP saw it.'\n"
            "# plot_type: 'table'\n"
            "# pconfig:\n"
            f"#     id: 'ddi_examples_stats_table_{id_}'\n"
            f"#     title: 'DDI Example Selection Stats ({id_})'\n"
            "Sample\tDDIs in\tDDIs kept\tExamples\tExamples per DDI\tProtein universe\t"
            "Contested parents\tShortlisted DDIs\tCandidate pairs\tCandidate examples\n"
        )
        for split in SPLITS:
            st = stats[split]
            per = st["examples"] / st["kept"] if st["kept"] else 0.0
            fh.write(
                f"{mqc_sample(id_, split)}\t{st['ddis_in']}\t{st['kept']}\t{st['examples']}\t"
                f"{per:.2f}\t{st['universe']}\t{st['contested']}\t{st['shortlisted']}\t"
                f"{st['cand_pairs']}\t{st['cand_examples']}\n"
            )


# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--train_ppis", required=True, help="train split CSV (Pfam family pairs)")
    ap.add_argument("--val_ppis", required=True)
    ap.add_argument("--test_ppis", required=True)
    ap.add_argument("--train_fasta", required=True, help="train split FASTA (domain instances)")
    ap.add_argument("--val_fasta", required=True)
    ap.add_argument("--test_fasta", required=True)
    ap.add_argument("--instances", required=True, help="instances.tsv: family/clan/instance/parent table")
    ap.add_argument(
        "--candidate-network",
        default=None,
        help="high-confidence negative family pairs. Their parents claim proteins in the same ILP, "
        "so their examples satisfy Barrier B too -- see this module's docstring for the caveat.",
    )
    ap.add_argument("--examples-target", type=int, default=5, help="N: cap on examples per DDI (default 5)")
    ap.add_argument(
        "--shortlist-factor",
        type=int,
        default=4,
        help="cap a DDI's candidate pool at this multiple of N before the ILP (default 4)",
    )
    ap.add_argument(
        "--candidate-factor",
        type=float,
        default=4.0,
        help="cap candidate_network pairs per split at this multiple of the split's DDI count "
        "(default 4), mirroring SAMPLE_NEGATIVES_ILP's own candidate cap",
    )
    ap.add_argument(
        "--lambda-diversity",
        type=float,
        default=0.1,
        help="weight on parent reuse. Must stay below 0.5: one extra example raises the overflow "
        "sum by at most 2, so 2*lambda < 1 is what keeps diversity from ever costing an example.",
    )
    ap.add_argument(
        "--no-barrier-b",
        action="store_true",
        help="let a parent protein be claimed by several splits at once. For split_method=random, "
        "which deliberately puts the same node in more than one split so the baseline shows the "
        "leak: enforcing Barrier B there would repair part of that leak and blunt the very "
        "comparison the naive baseline exists to make.",
    )
    ap.add_argument("--max-sec", type=int, default=300, help="total ILP time budget in seconds (default 300)")
    ap.add_argument("--solver", default=None, help="CVXPY solver name, e.g. GUROBI, SCIP (default: auto)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--id", required=True, help="Dataset ID, for MultiQC tagging")
    ap.add_argument("--verbose", action="store_true", help="let the solver print its own log")
    args = ap.parse_args()

    n = args.examples_target
    if n < 1:
        sys.exit("--examples-target must be at least 1")
    if args.lambda_diversity >= 0.5:
        sys.exit(
            f"--lambda-diversity must be < 0.5 (got {args.lambda_diversity}): at 0.5 or above the "
            f"diversity penalty can outweigh an example, which inverts the objective's priorities."
        )

    inst_rows = read_instances(args.instances)
    members = instances_by_family(inst_rows)
    parent_of = {r["instance_id"]: r["protein_id"] for r in inst_rows}
    all_parents = set(parent_of.values())
    print(
        f"{len(inst_rows):,} instances over {len(members):,} families and {len(all_parents):,} parent proteins",
        file=sys.stderr,
    )

    ppi_paths = {"train": args.train_ppis, "val": args.val_ppis, "test": args.test_ppis}
    fasta_paths = {"train": args.train_fasta, "val": args.val_fasta, "test": args.test_fasta}
    rows_by_split = {s: read_ppis(p) for s, p in ppi_paths.items()}
    split_ids = {s: set(read_fasta(p)) for s, p in fasta_paths.items()}

    shortlist_k = max(n, args.shortlist_factor * n)
    units = build_units("pos", rows_by_split, members, split_ids, parent_of, shortlist_k, args.seed)

    # candidate_network pairs are family pairs like the positives, and families
    # are split-exclusive, so each pair belongs to at most one split.
    cand_by_split = {s: [] for s in SPLITS}
    if args.candidate_network:
        # A set per family, not a single split: leakage-aware splits make families
        # exclusive, but split_method=random deliberately does not, and a pair whose
        # two families share more than one split belongs in each of them.
        fam_splits = defaultdict(set)
        for s in SPLITS:
            for row in rows_by_split[s]:
                fam_splits[row["protein1"]].add(s)
                fam_splits[row["protein2"]].add(s)
        n_unplaced = 0
        for row in read_ppis(args.candidate_network):
            shared = fam_splits.get(row["protein1"], set()) & fam_splits.get(row["protein2"], set())
            for s in shared:
                cand_by_split[s].append(row)
            if not shared:
                n_unplaced += 1
        for s in SPLITS:
            limit = int(math.ceil(args.candidate_factor * len(rows_by_split[s])))
            if len(cand_by_split[s]) > limit:
                rng = random.Random(f"{args.seed}:candidate_network:{s}")
                keep = rng.sample(range(len(cand_by_split[s])), limit)
                cand_by_split[s] = [cand_by_split[s][i] for i in sorted(keep)]
                print(
                    f"  {s}: candidate_network capped to {limit:,} pairs "
                    f"(--candidate-factor {args.candidate_factor})",
                    file=sys.stderr,
                )
        print(
            f"  candidate_network: {sum(len(v) for v in cand_by_split.values()):,} pairs placed in a "
            f"split, {n_unplaced:,} skipped (families in different splits or in none)",
            file=sys.stderr,
        )
        units += build_units("cand", cand_by_split, members, split_ids, parent_of, shortlist_k, args.seed)

    # A parent is contested when it carries candidates in more than one split;
    # only those need claim variables. splits_of also records every parent some
    # candidate reaches, which is what makes "unclaimed" mean genuinely never in play.
    splits_of = defaultdict(set)
    for u in units:
        for a, b in u.cands:
            splits_of[parent_of[a]].add(u.split)
            splits_of[parent_of[b]].add(u.split)
    contested = set() if args.no_barrier_b else {p for p, ss in splits_of.items() if len(ss) > 1}
    if args.no_barrier_b:
        print(
            "--no-barrier-b: parent proteins may be claimed by several splits at once. Every unit "
            "is therefore decoupled and the ILP is skipped.",
            file=sys.stderr,
        )
    print(
        f"{len(units):,} units ({sum(1 for u in units if u.kind == 'pos'):,} positive); "
        f"{len(splits_of):,} parent proteins in play, {len(contested):,} contested",
        file=sys.stderr,
    )

    free = [
        u for u in units if all(parent_of[a] not in contested and parent_of[b] not in contested for a, b in u.cands)
    ]
    # Split by identity, not equality: Unit is a plain dataclass, so two rows of
    # the same DDI in different splits would compare equal and cross-contaminate.
    free_ids = {id(u) for u in free}
    ilp_units = [u for u in units if id(u) not in free_ids]
    print(f"  {len(free):,} decoupled (local greedy), {len(ilp_units):,} in the ILP", file=sys.stderr)

    for u in free:
        rng = random.Random(f"{args.seed}:free:{u.kind}:{u.split}:{u.fam1}:{u.fam2}")
        u.picked = sorted(diverse_pick(u.cands, min(n, len(u.cands)), parent_of, rng))

    if ilp_units:
        print("Solving the selection ILP …", file=sys.stderr)
        solve_selection(
            ilp_units,
            parent_of,
            contested,
            n,
            args.lambda_diversity,
            args.max_sec,
            args.solver,
            args.seed,
            args.verbose,
        )
    else:
        why = "--no-barrier-b" if args.no_barrier_b else "no contested parent protein"
        print(f"ILP skipped ({why}): every unit's examples are decided locally.", file=sys.stderr)

    # The per-split protein universe: every uncontested parent in play for this
    # split (safe -- no other split can want it) plus the contested parents its
    # selected examples actually used. Deriving it from the selection rather than
    # from c[p,s] keeps it tight: c carries no objective weight, so the solver is
    # free to claim a contested parent it never uses, which would deny it to a
    # split that could.
    universe = {s: set() for s in SPLITS}
    for p, ss in splits_of.items():
        if p in contested:
            continue
        # Exactly one split, unless --no-barrier-b let a parent stay in several.
        for s in ss:
            universe[s].add(p)
    for u in units:
        for a, b in u.picked:
            for p in (parent_of[a], parent_of[b]):
                if p in contested:
                    universe[u.split].add(p)

    stats = {}
    for split in SPLITS:
        pos = [u for u in units if u.kind == "pos" and u.split == split]
        cand = [u for u in units if u.kind == "cand" and u.split == split]
        kept = [u for u in pos if u.picked]
        write_ppi_csv([u.row for u in kept], f"{split}_sel.csv")
        n_examples = write_examples(kept, f"{split}_examples.csv")
        n_cand_examples = write_examples([u for u in cand if u.picked], f"{split}_candidate_examples.csv")
        write_ids(universe[split], f"{split}_universe.txt")
        stats[split] = {
            "ddis_in": len(pos),
            "kept": len(kept),
            "dropped": len(pos) - len(kept),
            "full": sum(1 for u in pos if len(u.picked) >= n),
            "partial": sum(1 for u in pos if 0 < len(u.picked) < n),
            "examples": n_examples,
            "universe": len(universe[split]),
            "contested": len(universe[split] & contested),
            "shortlisted": sum(1 for u in pos if u.n_raw > len(u.cands)),
            "cand_pairs": len(cand),
            "cand_examples": n_cand_examples,
        }
        st = stats[split]
        print(
            f"  {split}: {st['kept']:,}/{st['ddis_in']:,} DDIs kept "
            f"({st['full']:,} at full N, {st['partial']:,} partial, {st['dropped']:,} dropped), "
            f"{st['examples']:,} examples, {st['universe']:,} proteins claimed",
            file=sys.stderr,
        )

    # Parents no candidate ever reached. A contested parent the ILP left
    # unclaimed is deliberately in neither list -- handing it to a split now
    # would reintroduce the leak Barrier B just prevented.
    write_ids(all_parents - set(splits_of), "unclaimed.txt")
    print(
        f"{len(all_parents - set(splits_of)):,} parent proteins never in play "
        f"(free for the negative sampler to extend a split's universe with)",
        file=sys.stderr,
    )

    write_mqc(stats, args.id)


if __name__ == "__main__":
    main()
