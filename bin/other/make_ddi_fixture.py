#!/usr/bin/env python3
"""Generate `data/test_ddis.csv`, the committed DDI smoke-test input (`-profile test_ddi`).

Not part of the pipeline -- run by hand when the fixture needs regenerating.

The fixture is a list of **real Pfam family pairs**, and nothing else: the
test_ddi profile supplies no sequences, no species table and no instances.tsv,
so a run exercises FETCH_DOMAIN_META for real -- the Pfam-A.fasta.gz stream, the
clan table, the tier cascade and the instance sampling included. That costs a
~6.3 GB download per cold run, which is the price of testing the fetch at all.

Accessions come from Pfam's own clan table rather than from memory, so a typo or
a dead family cannot creep in. The composition is deliberate:

  * `--clanned-clans` clans contributing `--per-clan` families each. Families
    sharing a clan are genuinely homologous, so BLAST finds them, MAKE_METIS
    contracts them onto one node, and its same-clan skip is actually exercised.
    Clans are taken evenly spaced down the size ranking, so the set spans large
    and small ones instead of only the giants.
  * `--unclanned` families with no clan at all. 16,201 of Pfam's 30,134 families
    are unclanned and each becomes its own singleton clan, so this is the common
    case, not an edge case -- and it is what Risk 4 (clan dust) is about.
  * One deliberately nonexistent accession, so the missing-family report has
    something to report and a run proves it does not abort on one.

The pairing has structure on purpose: a same-clan DDI lives inside a single
contracted node and therefore always survives cluster assignment, while a
cross-clan one only survives if the ILP co-assigns the two clusters. A purely
random pairing would be almost entirely cross-cluster and the attrition numbers
would say more about the fixture than about the pipeline.

Usage:
    python bin/other/make_ddi_fixture.py --out data/test_ddis.csv
    python bin/other/make_ddi_fixture.py --clans /path/to/Pfam-A.clans.tsv.gz    # offline
"""

import argparse
import csv
import gzip
import io
import os
import random
import sys
import urllib.request

CLANS_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz"

# Not a Pfam accession, and not in Pfam-A.dead either: the "typo" branch of
# fetch_domains.py's missing-family report rather than the "killed family" one.
MISSING_ACCESSION = "PF99999"


def load_clans(path):
    """Return (clanned {clan: [family]}, unclanned [family]) from Pfam-A.clans.tsv(.gz)."""
    if path:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            text = fh.read()
    else:
        print(f"Downloading {CLANS_URL} ...", file=sys.stderr)
        with urllib.request.urlopen(CLANS_URL, timeout=300) as resp:
            raw = resp.read()
        text = gzip.GzipFile(fileobj=io.BytesIO(raw)).read().decode("utf-8", "replace")

    clanned, unclanned = {}, []
    for line in text.splitlines():
        cols = line.split("\t")
        if len(cols) < 2 or not cols[0].startswith("PF"):
            continue
        family, clan = cols[0].strip(), cols[1].strip()
        if clan:
            clanned.setdefault(clan, []).append(family)
        else:
            unclanned.append(family)
    return clanned, unclanned


def pick_families(clanned, unclanned, n_clans, per_clan, n_unclanned):
    """Deterministic selection -- no RNG, so it cannot drift with the Python version."""
    eligible = sorted((c for c, fams in clanned.items() if len(fams) >= per_clan), key=lambda c: (-len(clanned[c]), c))
    if len(eligible) < n_clans:
        sys.exit(f"only {len(eligible)} clans have >= {per_clan} families")
    # Evenly spaced down the size ranking rather than the top N, so the set is
    # not exclusively P-loop-NTPase-sized superclans.
    step = len(eligible) / n_clans
    clans = [eligible[int(i * step)] for i in range(n_clans)]

    families, clan_of = [], {}
    for clan in clans:
        for family in sorted(clanned[clan])[:per_clan]:
            families.append(family)
            clan_of[family] = clan

    singles = sorted(unclanned)
    step = len(singles) / n_unclanned
    for i in range(n_unclanned):
        family = singles[int(i * step)]
        families.append(family)
        clan_of[family] = f"CL_{family}"  # what fetch_domains.py will synthesise
    return families, clan_of


def build_ddis(rng, families, clan_of, n_ddis, frac_same_clan, frac_self, include_missing):
    """Return sorted, deduplicated (family1, family2) pairs."""
    by_clan = {}
    for family in families:
        by_clan.setdefault(clan_of[family], []).append(family)
    multi = [c for c, fams in by_clan.items() if len(fams) > 1]

    pairs, seen = [], set()
    guard = 0
    while len(pairs) < n_ddis and guard < n_ddis * 50:
        guard += 1
        roll = rng.random()
        if roll < frac_self:
            a = b = rng.choice(families)
        elif roll < frac_self + frac_same_clan and multi:
            # Same contracted node: survives cluster assignment unconditionally.
            a, b = rng.sample(by_clan[rng.choice(multi)], 2)
        else:
            a, b = rng.sample(families, 2)
        key = tuple(sorted((a, b)))
        if key in seen:
            continue
        seen.add(key)
        pairs.append(key)

    if include_missing:
        pairs.append(tuple(sorted((MISSING_ACCESSION, families[0]))))
    return sorted(pairs)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="data/test_ddis.csv")
    ap.add_argument("--clans", help="local Pfam-A.clans.tsv(.gz); downloaded when omitted")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--clanned-clans", type=int, default=20, help="clans contributing families (default 20)")
    ap.add_argument("--per-clan", type=int, default=4, help="families taken from each clan (default 4)")
    ap.add_argument("--unclanned", type=int, default=12, help="unclanned families, i.e. singleton clans")
    ap.add_argument("--ddis", type=int, default=250)
    ap.add_argument("--same-clan-fraction", type=float, default=0.45)
    ap.add_argument("--self-fraction", type=float, default=0.06)
    ap.add_argument("--no-missing", action="store_true", help="omit the nonexistent-accession DDI")
    args = ap.parse_args()

    clanned, unclanned = load_clans(args.clans)
    families, clan_of = pick_families(clanned, unclanned, args.clanned_clans, args.per_clan, args.unclanned)
    ddis = build_ddis(
        random.Random(args.seed),
        families,
        clan_of,
        args.ddis,
        args.same_clan_fraction,
        args.self_fraction,
        not args.no_missing,
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        # LF, not csv's excel-dialect CRLF: this is a committed text fixture and
        # data/test_ppis.csv next to it is LF.
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["protein1", "protein2"])
        writer.writerows(ddis)

    same_clan = sum(1 for a, b in ddis if a != b and clan_of.get(a) and clan_of.get(a) == clan_of.get(b))
    print(
        f"{len(families)} families over {len(set(clan_of.values()))} clans "
        f"({args.unclanned} unclanned, i.e. singleton clans)",
        file=sys.stderr,
    )
    print(
        f"{len(ddis)} DDIs: {same_clan} same-clan, {sum(1 for a, b in ddis if a == b)} self, "
        f"{0 if args.no_missing else 1} with a nonexistent family",
        file=sys.stderr,
    )
    print(f"Written {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
