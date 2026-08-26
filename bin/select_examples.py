#!/usr/bin/env python3
"""
DDI mode: turn family-pair DDIs into domain-instance-pair examples, each parent
protein used by at most one split.

The split CSVs handed over by SOLVE_ILP/SORT_PPIS/SPLIT_RANDOM (and filtered by
REMOVE_REDUNDANT) hold Pfam *family* pairs, but a row a classifier can train on
is a pair of concrete domain *instances* -- and the parent proteins of those two
instances must not turn up in another split, or the parent's other domains carry
the interaction across the split boundary.

Selection and that one-split-per-parent rule are one problem, not two: claiming a
parent protein for one split removes it from every other split's options, so which
examples a split can still reach depends on what the other splits took. Greedy
order would decide the outcome, so it is solved as an ILP over all three splits
jointly.

Three reductions keep that ILP small, none of them an approximation:

  * Only *contested* parents -- proteins carrying candidates in more than one
    split -- get claim variables. A parent in play for a single split can be
    claimed for free, so every DDI whose candidates touch no contested parent is
    decoupled from the rest of the problem and is settled by a local diversity
    greedy instead.
  * What is left is then cut into connected components over the contested parents
    they share (component_partition), and each component is solved on its own.
    Every constraint is either per-unit or per-parent and the objective is a sum
    of per-unit terms, so units sharing no contested parent share nothing -- this
    is the same reduction as the point above, carried one step further from
    "touches no contested parent" to "shares no contested parent, transitively".
    It is what makes the target scale (~90k surviving DDIs, so millions of
    variables in one model) tractable at all.
  * A per-DDI shortlist caps the candidate pool at K = shortlist_factor * N. At
    the default pool size (M = N = 5, so 25 candidates) it trims almost nothing
    and exists as a guard for a larger --ddi_examples_pool_factor.

Within a component the objective is lexicographic, solved as one bounded stage
per level, so every coefficient stays 1 or lambda rather than a constant large
enough to outrank its own tail:

  1. keep as many positive DDIs as possible at >= 1 example (a DDI that reaches
     zero is dropped outright, so this level is what stops the solver starving
     one DDI to complete another),
  2. then maximise the positive DDIs reaching the full N,
  3. then total positive examples, less the diversity penalty on parent reuse,
  4. then the same for the candidate_network negatives, which ride in the same
     claim accounting so their examples inherit the one-split-per-parent rule by
     construction -- but strictly below the positives, so a negative can never
     take a parent a positive needed.

A component too large for --max-ilp-candidates, or reached after the --max-sec
budget is spent, drops to a deterministic greedy (greedy_component) that keeps
one-split-per-parent exactly and loses only optimality. That path is counted and
reported, never silent.

--allow-shared-parents turns the whole claim mechanism off, which is what
split_method=random needs: that path deliberately puts a node in more than one
split so the naive baseline shows the leak, and holding parents to one split each
would repair part of it. Every unit is then decoupled and the ILP is skipped
entirely.

Note the caveat that follows from step 3: a candidate pair SAMPLE_NEGATIVES does
not ultimately select will have claimed proteins for nothing, and purely random
negatives stay outside this ILP altogether -- so the rule is exact for
high-confidence negatives and holds via the per-split protein universe for the
rest.
"""

import argparse
import math
import os
import random
import statistics
import sys
import time
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
from solve_ilp import _solver_options  # noqa: E402

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


_WARNED = set()


def _warn_once(msg, key=None):
    """Print a warning the first time only, with the detail of that first time.

    The component decomposition turns what used to be three solves per task into
    four per component, so a per-solve warning would run to thousands of near
    identical lines in .command.err. `key` is what makes them one warning rather
    than one per component: the message carries the numbers of the case that
    tripped it, the key carries only its kind. The totals are reported at the end
    and in the DDI Example Selection ILP MultiQC table.
    """
    key = msg if key is None else key
    if key not in _WARNED:
        _WARNED.add(key)
        print(msg, file=sys.stderr)


@dataclass
class SolverEnv:
    """One solver environment reused across every component solve.

    cvxpy builds a fresh Gurobi environment per solve unless it is handed one, and
    each build re-reads and re-validates the licence. That was paid three times per
    task before the decomposition and would be paid ~4x per component after it.

    Self-disabling: if the installed cvxpy or solver does not take an `env`
    argument the first solve raises, the environment is dropped, and every later
    solve runs the default way. The open solvers have no such parameter and are
    unaffected either way.
    """

    env: object = None

    def kwargs(self):
        return {"env": self.env} if self.env is not None else {}

    def disable(self, exc):
        print(
            f"Note: the solver did not accept a shared environment ({exc!r}); building one per "
            f"solve instead. Slower at many components, identical in result.",
            file=sys.stderr,
        )
        self.env = None


def make_solver_env(solver):
    """A shared gurobipy.Env when the solver is Gurobi, an inert holder otherwise."""
    if (solver or "").upper() != cp.GUROBI:
        return SolverEnv(None)
    try:
        import gurobipy

        return SolverEnv(gurobipy.Env())
    except Exception as exc:  # gurobipy absent, or the licence check failed here
        print(f"Note: no shared Gurobi environment ({exc!r}); each solve builds its own.", file=sys.stderr)
        return SolverEnv(None)


@dataclass
class SolveContext:
    """Solver settings and run-level counters, threaded through every stage solve."""

    solver: str = None
    seed: int = 42
    verbose: bool = False
    env: SolverEnv = field(default_factory=SolverEnv)
    stages: int = 0  # stage solves attempted
    timeouts: int = 0  # of those, the ones that returned an incumbent at the time limit


def _run_stage(problem, label, secs, ctx):
    """Solve one lexicographic stage; return its objective value or None.

    None means the solver came back without a usable solution -- with a per-component
    share of the budget the ordinary cause is a limit reached before any incumbent --
    and the caller sends that component to the greedy fallback. A solver that is
    missing or misconfigured raises instead, which still stops the run rather than
    quietly degrading every component.
    """
    # solve_ilp._solver_options maps both the time limit and the seed to the names
    # the chosen solver actually takes -- time_limit is CVXPY's HiGHS-specific
    # kwarg, Gurobi wants TimeLimit and SCIP wants scip_params={'limits/time': ...},
    # and SCIP rejects the bare kwarg outright rather than ignoring it. The generic
    # name is only used as a last resort for a solver it does not know.
    kwargs = {"verbose": ctx.verbose}
    if ctx.solver:
        solver_opts = _solver_options(ctx.solver, ctx.seed, secs)
        if not solver_opts:
            _warn_once(
                f"Warning: no solver-specific options are known for solver {ctx.solver}, so its time "
                f"limit is passed under CVXPY's generic name (may be rejected) and its tie-breaking is "
                f"unseeded, so this selection is not reproducible run to run."
            )
            solver_opts = {"time_limit": secs}
        kwargs.update(solver_opts)
        kwargs["solver"] = ctx.solver
    else:
        _warn_once(
            "Warning: no --solver given, so CVXPY picks one and its internal randomisation stays "
            "unseeded; pass --solver for a reproducible selection."
        )
        kwargs["time_limit"] = secs

    ctx.stages += 1
    try:
        problem.solve(**kwargs, **ctx.env.kwargs())
    except Exception as exc:
        if not ctx.env.kwargs():
            raise
        ctx.env.disable(exc)
        problem.solve(**kwargs)

    if problem.status not in cp.settings.SOLUTION_PRESENT:
        _warn_once(
            f"Warning: stage '{label}' came back with solver status {problem.status} and no usable "
            f"solution at a {secs:.0f}s limit; that component falls back to the greedy. Raise "
            f"--max-sec if this affects many components.",
            key=f"nosolution:{label}",
        )
        return None
    if problem.status == cp.settings.USER_LIMIT:
        ctx.timeouts += 1
        _warn_once(
            f"Warning: stage '{label}' hit its {secs:.0f}s limit before proving optimality and is "
            f"using the best incumbent found. Counted in the DDI Example Selection ILP table.",
            key=f"userlimit:{label}",
        )
    if ctx.verbose:
        print(f"  stage '{label}': objective {problem.value:,.4g} ({problem.status})", file=sys.stderr)
    return float(problem.value)


def solve_selection(units, parent_of, contested, n, lam, max_sec, ctx):
    """Choose <= n examples per unit, one split per parent protein. Fills Unit.picked.

    Variables
        y[e]   in {0,1}  candidate example e is selected
        nz[d]  in {0,1}  unit d ends up with >= 1 example
        r[d]   in {0,1}  unit d reaches the full n
        c[p,s] in {0,1}  contested parent p is claimed by split s
        o[d,p] >= 0      how often unit d reuses parent p beyond the first time

    Constraints
        (1) sum_{e in d} y[e] <= min(n, |C_d|)          per-unit cap (n is a cap, not a quota)
        (2) sum_{e in d} y[e] >= nz[d]
        (3) sum_{e in d} y[e] >= n * r[d]
        (4) sum_s c[p,s] <= 1                           one split per parent protein
        (5) y[e] <= c[parent(e), split(d(e))]           both parents, contested only
        (6) o[d,p] >= (uses of p by d's selected examples) - 1

    Returns True when the selection was solved, False when a stage came back with
    no usable solution and the caller should fall back to the greedy. Unit.picked
    is only written on the True path, so a False leaves the units untouched.
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
        return True

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
        for pi, (_p, idxs) in enumerate(sorted(by_prot.items())):
            for idx in idxs:
                rows.append(pi)
                cols.append(idx)
        one_split_per_parent = sp.coo_matrix(
            (np.ones(len(rows)), (rows, cols)), shape=(len(by_prot), len(claim_index))
        ).tocsr()
        cons.append(one_split_per_parent @ c <= np.ones(len(by_prot)))

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
    #
    # A key with a single entry is dropped: its row reads y_k - 1 <= o with
    # y_k in {0,1}, so the left side is never positive and o = 0 at every optimum.
    # The row and its variable are unconditionally slack, so removing them cannot
    # change the solution -- and they are the common case for a lopsided DDI (a
    # 1 x 5 family pair leaves five of its six keys droppable).
    key_count = defaultdict(int)
    for ui, (a, b) in cand_index:
        for p in (parent_of[a], parent_of[b]):
            key_count[(ui, p)] += 1
    over_index, o_rows, o_cols = {}, [], []
    for k, (ui, (a, b)) in enumerate(cand_index):
        for p in (parent_of[a], parent_of[b]):
            key = (ui, p)
            if key_count[key] < 2:
                continue
            if key not in over_index:
                over_index[key] = len(over_index)
            o_rows.append(over_index[key])
            o_cols.append(k)
    o = None
    if over_index:
        o = cp.Variable(len(over_index), nonneg=True)
        per_parent = sp.coo_matrix((np.ones(len(o_rows)), (o_rows, o_cols)), shape=(len(over_index), n_cand)).tocsr()
        cons.append(per_parent @ y - 1.0 <= o)

    pos_u = [i for i, u in enumerate(units) if u.kind == "pos"]
    cand_u = [i for i, u in enumerate(units) if u.kind == "cand"]
    pos_y = [k for k, (ui, _) in enumerate(cand_index) if units[ui].kind == "pos"]
    cand_y = [k for k, (ui, _) in enumerate(cand_index) if units[ui].kind == "cand"]
    pos_o = [idx for (ui, _), idx in over_index.items() if units[ui].kind == "pos"]
    cand_o = [idx for (ui, _), idx in over_index.items() if units[ui].kind == "cand"]

    # One objective level per stage, every coefficient 1 or lam. An earlier
    # formulation folded levels 2 and 3 into one expression scaled by
    # big = sum(caps) + 1, which is ~50,001 at 10k DDIs and N = 5 -- putting the
    # objective around 5e8 while the freeze tolerance below stayed absolute, five
    # orders of magnitude under the solver's own feasibility slack at that
    # magnitude. Splitting the stage removes the constant rather than tuning it.
    stages = []
    if pos_u:
        stages.append(("keep every positive DDI", _sum(nz, pos_u), 0.15))
        stages.append(("fill positives to N", _sum(r, pos_u), 0.35))
        stages.append(("maximise positive examples", _sum(y, pos_y) - lam * _sum(o, pos_o), 0.25))
    if cand_u:
        stages.append(("expand candidate negatives", _sum(y, cand_y) - lam * _sum(o, cand_o), 0.25))

    total_share = sum(s for _, _, s in stages)
    for label, expr, share in stages:
        secs = max(1.0, max_sec * share / total_share)
        value = _run_stage(cp.Problem(cp.Maximize(expr), cons), label, secs, ctx)
        if value is None:
            return False
        # Freeze this level before optimising the next. An incumbent from a
        # timed-out stage is still achievable, so the bound stays satisfiable.
        # Relative, so the slack tracks the objective's own magnitude.
        cons = cons + [expr >= value - max(1e-4, 1e-6 * abs(value))]

    y_val = np.asarray(y.value).ravel()
    for k, (ui, pair) in enumerate(cand_index):
        if y_val[k] > 0.5:
            units[ui].picked.append(pair)
    for ui, u in enumerate(units):
        # Trim defensively: a solver returning 0.5+eps on more candidates than
        # the cap would otherwise leak an extra example into the output.
        u.picked = sorted(u.picked)[: caps[ui]]
    return True


# ---------------------------------------------------------------------------
# Decomposition, fallback and the driver over both
# ---------------------------------------------------------------------------


def component_partition(units, parent_of, contested):
    """Cut the ILP units into independent sub-problems over the parents they share.

    Nodes are units and contested parents, with an edge wherever a unit has a
    candidate pair on that parent; the connected components of that graph are
    sub-problems that share no variable and no constraint. Every constraint is
    either per-unit (the cap, the >= nz and >= n*r rows, the overflow rows) or
    per-parent (one split per parent, and y <= c), and the objective is a sum of
    per-unit terms -- so maximising each component's stage separately and freezing
    each component's own optimum is identical to doing it globally.

    A contested parent therefore lives in exactly one component, which is what lets
    the greedy fallback keep its claim bookkeeping local.

    Returns lists of Units, largest candidate count first, so the components that
    dominate the cost are solved before any stop-loss can fire.
    """
    uf = list(range(len(units)))

    def find(x):
        while uf[x] != x:
            uf[x] = uf[uf[x]]  # path halving
            x = uf[x]
        return x

    first_seen = {}
    for ui, u in enumerate(units):
        for a, b in u.cands:
            for p in (parent_of[a], parent_of[b]):
                if p not in contested:
                    continue
                if p in first_seen:
                    ra, rb = find(ui), find(first_seen[p])
                    if ra != rb:
                        uf[ra] = rb
                else:
                    first_seen[p] = ui

    groups = defaultdict(list)
    for ui in range(len(units)):
        groups[find(ui)].append(ui)
    # Sorted by size then by the lowest unit index, so the order does not depend on
    # dict iteration and two runs of the same input solve the same problems in the
    # same sequence -- which matters once a wall-clock stop-loss decides where the
    # ILP stops and the greedy starts.
    keyed = sorted((-sum(len(units[i].cands) for i in idxs), idxs[0], idxs) for idxs in groups.values())
    return [[units[i] for i in idxs] for _, _, idxs in keyed]


def greedy_component(units, parent_of, contested, n, seed, tag):
    """Deterministic fallback for a component the ILP cannot be run on.

    Positives before candidate_network pairs -- the objective's one hard priority --
    and a seeded shuffle within each kind, which interleaves the three splits so no
    split gets first refusal on every contested parent. Ordering by split size
    instead would hand the largest split every contested parent it can use rather
    than its share of them.

    Each unit then picks from the candidates whose parents no *other* split has
    claimed yet, and claims what it took. Claims stay local to the component because
    component_partition puts every unit touching a given contested parent in one
    component.

    One split per parent protein is preserved exactly; what is lost is optimality,
    which is the whole reason the ILP exists -- so every call must be counted and
    reported, never silently substituted.
    """
    rng = random.Random(f"{seed}:greedy:{tag}")
    order = []
    for kind in ("pos", "cand"):
        group = [u for u in units if u.kind == kind]
        rng.shuffle(group)
        order += group

    claimed = {}
    for u in order:
        usable = [
            (a, b)
            for a, b in u.cands
            if claimed.get(parent_of[a], u.split) == u.split and claimed.get(parent_of[b], u.split) == u.split
        ]
        # Seeded from the unit itself, like the decoupled path, so a unit's own
        # choice among what is still available does not depend on the shuffle.
        pick_rng = random.Random(f"{seed}:greedy:{u.kind}:{u.split}:{u.fam1}:{u.fam2}")
        u.picked = sorted(diverse_pick(usable, min(n, len(usable)), parent_of, pick_rng))
        for a, b in u.picked:
            for p in (parent_of[a], parent_of[b]):
                if p in contested:
                    claimed[p] = u.split


def run_selection(units, parent_of, contested, n, lam, max_sec, max_cand, ctx):
    """Solve every component, with a candidate cap and a wall-clock stop-loss.

    The time budget is shared out in proportion to each component's candidate
    count, recomputed against what is left after each one, so a component that
    finishes early hands its unused seconds to the rest. Once the budget is spent
    the remaining components go to the greedy rather than the task overrunning its
    allocation -- at real DDI counts --max-sec (params.ddi_select_max_sec) is the
    knob to raise if that starts happening.

    Returns the numbers the DDI Selection ILP MultiQC table reports.
    """
    comps = component_partition(units, parent_of, contested)
    sizes = [sum(len(u.cands) for u in c) for c in comps]
    n_units = [len(c) for c in comps]
    print(
        f"  {len(comps):,} independent component(s); units min/median/max "
        f"{min(n_units):,}/{statistics.median(n_units):,.0f}/{max(n_units):,}, candidates "
        f"{min(sizes):,}/{statistics.median(sizes):,.0f}/{max(sizes):,}. Largest: "
        f"{sizes[0]:,} candidates over {n_units[0]:,} units.",
        file=sys.stderr,
    )

    started = time.monotonic()
    fallback, fallback_units, remaining_cand = 0, 0, sum(sizes)
    for ci, comp in enumerate(comps):
        elapsed = time.monotonic() - started
        kind, reason = None, None
        if sizes[ci] > max_cand:
            kind = "candidate cap"
            reason = f"{sizes[ci]:,} candidates is over --max-ilp-candidates ({max_cand:,})"
        elif elapsed >= max_sec:
            kind = "time budget"
            reason = f"the {max_sec}s ILP budget ran out after {ci:,} of {len(comps):,} components"
        if kind is None:
            share = (max_sec - elapsed) * sizes[ci] / remaining_cand if remaining_cand else max_sec
            if not solve_selection(comp, parent_of, contested, n, lam, share, ctx):
                kind = "no solution"
                reason = "the solver returned no usable solution"
        remaining_cand -= sizes[ci]
        if kind is not None:
            fallback += 1
            fallback_units += len(comp)
            _warn_once(
                f"WARNING: at least one component fell back to the greedy selection ({reason}). "
                f"One split per parent protein still holds, but that component's selection is no "
                f"longer optimal. The DDI Example Selection ILP MultiQC table reports the total.",
                key=f"fallback:{kind}",
            )
            greedy_component(comp, parent_of, contested, n, ctx.seed, str(ci))

    info = {
        "units": len(units),
        "components": len(comps),
        "largest_units": n_units[0],
        "largest_cands": sizes[0],
        "fallback": fallback,
        "fallback_units": fallback_units,
        "timeouts": ctx.timeouts,
        "seconds": time.monotonic() - started,
    }
    print(
        f"  ILP done in {info['seconds']:,.1f}s over {ctx.stages:,} stage solves "
        f"({info['timeouts']:,} hit their time limit); {fallback:,} component(s) "
        f"covering {fallback_units:,} units fell back to the greedy.",
        file=sys.stderr,
    )
    return info


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
        # family1/family2 are excluded as well as the id columns: an input CSV --
        # a --candidate-network file in particular, since those come from outside
        # the pipeline -- that happens to carry a family1 column would otherwise
        # overwrite the computed family pair, poisoning utils.read_family_pairs()
        # and with it train_classifier.py's family aggregation and
        # bias_analysis.py's node_pairs. Matches expand_negatives.write_expanded.
        extra = {k: v for k, v in u.row.items() if k not in ("protein1", "protein2", "family1", "family2")}
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


def partition_reserve(unclaimed, weights, seed, allow_shared):
    """Hand each never-in-play parent protein to exactly one split.

    Partitioned here rather than in each EXPAND_NEGATIVES task because this is
    the one place that knows all three splits: weighting the draw by the DDIs each
    split actually kept is correct on every split method, where weighting it by
    train_split/val_split/test_split was not -- split_method=kahip never reads
    those fractions (sort_ppis.py ranks the KaHIP blocks and hands largest ->
    train), so an 80/10/10 reserve could sit on top of a realised 50/30/20 split.

    Doing it once also drops the requirement that four independent tasks reach the
    same partition from the same file and seed. That was true, but it was a
    property that could quietly stop being true.

    Under --allow-shared-parents (split_method=random) every split gets the whole
    pool, since overlapping universes are the point on that path.
    """
    if allow_shared:
        return {s: set(unclaimed) for s in SPLITS}
    out = {s: set() for s in SPLITS}
    active = [s for s in SPLITS if weights.get(s, 0) > 0]
    if not active:
        return out
    total = float(sum(weights[s] for s in active))
    rng = random.Random(f"{seed}:reserve")
    for p in sorted(unclaimed):
        x = rng.random() * total
        acc = 0.0
        for s in active:
            acc += weights[s]
            if x < acc:
                out[s].add(p)
                break
    return out


def write_mqc(stats, ilp, id_):
    """Three MultiQC sections: the per-split DDI outcome bar, a per-split stats
    table, and one row per dataset describing the selection ILP itself.

    The ILP numbers are deliberately not columns on the per-split table: a
    component is cut over contested parents, which span splits by definition, so
    there is no per-split value to report there.
    """
    with open("select_examples_bar_mqc.tsv", "w") as fh:
        fh.write(
            f"# id: 'ddi_examples_bar_{id_}'\n"
            f"# section_name: 'DDI Example Selection: {id_}'\n"
            "# description: 'Positive DDIs per split, coloured by how many domain-instance "
            "examples each one ended up with. A DDI that reached zero examples -- every "
            "candidate blocked because another split already had its parent proteins, or the "
            "family had no usable instance -- is dropped from the split.'\n"
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

    with open("select_examples_ilp_mqc.tsv", "w") as fh:
        fh.write(
            # Its own section per dataset, exactly like the two tables above: custom
            # content sharing an id across files would have to be merged by MultiQC,
            # and a section that fails to parse costs the whole report, not one table.
            f"# id: 'ddi_examples_ilp_{id_}'\n"
            f"# section_name: 'DDI Example Selection ILP: {id_}'\n"
            "# description: 'Shape of the selection ILP, one row per dataset. Units are the "
            "interactions left after the reduction that settles anything touching no contested "
            "parent locally; those are cut into independent components over the contested parents "
            "they share, and each component is solved on its own. A component that exceeds "
            "--max-ilp-candidates, or that is reached after the time budget is spent, falls back to "
            "a deterministic greedy: one split per parent protein still holds there, but that "
            "component is no longer optimal, so a nonzero fallback count is worth acting on.'\n"
            "# plot_type: 'table'\n"
            "# pconfig:\n"
            f"#     id: 'ddi_examples_ilp_table_{id_}'\n"
            f"#     title: 'DDI Example Selection ILP ({id_})'\n"
            "Sample\tUnits in ILP\tComponents\tLargest component (units)\tLargest component (candidates)\t"
            "Greedy fallback components\tGreedy fallback units\tStages at time limit\tILP seconds\n"
        )
        fh.write(
            f"{id_}\t{ilp['units']}\t{ilp['components']}\t{ilp['largest_units']}\t{ilp['largest_cands']}\t"
            f"{ilp['fallback']}\t{ilp['fallback_units']}\t{ilp['timeouts']}\t{ilp['seconds']:.1f}\n"
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
        "so their examples respect the one-split-per-parent rule too -- see this module's "
        "docstring for the caveat.",
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
        "--allow-shared-parents",
        action="store_true",
        help="let a parent protein be claimed by several splits at once. For split_method=random, "
        "which deliberately puts the same node in more than one split so the baseline shows the "
        "leak: holding parents to one split each there would repair part of that leak and blunt "
        "the very comparison the naive baseline exists to make.",
    )
    ap.add_argument(
        "--max-sec",
        type=int,
        default=300,
        help="total ILP time budget in seconds (default 300), shared out across the components in "
        "proportion to their candidate counts. Components reached after it is spent fall back to "
        "the greedy, so raise it rather than let that happen at scale.",
    )
    ap.add_argument(
        "--max-ilp-candidates",
        type=int,
        default=200000,
        help="skip the ILP for any component with more candidate pairs than this and use the "
        "deterministic greedy instead (default 200000). A guard against a single component that "
        "would exhaust memory during canonicalisation; the fallback keeps one split per parent "
        "protein, it only loses optimality.",
    )
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
    contested = set() if args.allow_shared_parents else {p for p, ss in splits_of.items() if len(ss) > 1}
    if args.allow_shared_parents:
        print(
            "--allow-shared-parents: parent proteins may be claimed by several splits at once. "
            "Every unit is therefore decoupled and the ILP is skipped.",
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

    ilp = dict.fromkeys(
        ("units", "components", "largest_units", "largest_cands", "fallback", "fallback_units", "timeouts"), 0
    )
    ilp["seconds"] = 0.0
    if ilp_units:
        print("Solving the selection ILP …", file=sys.stderr)
        ctx = SolveContext(solver=args.solver, seed=args.seed, verbose=args.verbose, env=make_solver_env(args.solver))
        ilp = run_selection(
            ilp_units,
            parent_of,
            contested,
            n,
            args.lambda_diversity,
            args.max_sec,
            args.max_ilp_candidates,
            ctx,
        )
    else:
        why = "--allow-shared-parents" if args.allow_shared_parents else "no contested parent protein"
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
        # Exactly one split, unless --allow-shared-parents let a parent stay in several.
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

    # Parents no candidate ever reached, published whole as a diagnostic and split
    # three ways for the negative sampler. A contested parent the ILP left
    # unclaimed is deliberately in neither -- handing it to a split now
    # would reintroduce the leak the one-split-per-parent rule just prevented.
    #
    # This reserve exists so a negative family pair whose own split universe
    # cannot reach N examples has somewhere else to draw from without breaking
    # one-split-per-parent: a parent here belongs to no split, so giving it to one
    # takes it from none. Each reserve protein goes to exactly one split, weighted
    # by the DDIs that split actually kept -- see partition_reserve().
    #
    # It is all but inert at --examples-pool-factor 1 (M = N). A reserve parent is
    # only usable downstream if one of its instances is also in the target split's
    # FASTA, and every instance in that FASTA already has its parent recorded in
    # splits_of: the instance's family is in the split CSV, and pair_candidates()
    # enumerates the full cross product of both families' available instances, so
    # every available instance of every family in the split turns up in some
    # candidate pair. The only escape is a family whose every DDI had a partner
    # family with zero available instances.
    #
    # Above factor 1 the mechanism becomes live: a family then holds M = N * factor
    # instances while selection claims at most N per DDI, so genuinely unclaimed
    # parents exist and expand_negatives' extension path fires. The target
    # configuration is factor 2-3, so this is not dead code -- it is code whose
    # only exercise so far has been at the one setting where it cannot trigger.
    # Treat its first factor >= 2 run as untested: watch the "Extended with
    # unclaimed" MQC column go non-zero and confirm check_ddi_invariants.py still
    # passes on that run.
    unclaimed = all_parents - set(splits_of)
    write_ids(unclaimed, "unclaimed.txt")
    reserve = partition_reserve(unclaimed, {s: stats[s]["kept"] for s in SPLITS}, args.seed, args.allow_shared_parents)
    for split in SPLITS:
        write_ids(reserve[split], f"{split}_reserve.txt")
    print(
        f"{len(unclaimed):,} parent proteins never in play, split "
        + ", ".join(f"{len(reserve[s]):,} {s}" for s in SPLITS)
        + " (free for the negative sampler to extend a split's universe with)",
        file=sys.stderr,
    )

    write_mqc(stats, ilp, args.id)


if __name__ == "__main__":
    main()
