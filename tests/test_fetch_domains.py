#!/usr/bin/env python3
"""Checks for bin/fetch_domains.py's regions-based instance sampling.

The reason this file exists: instances used to come from ``Pfam-A.fasta.gz``, which
Pfam publishes 90 % non-redundant. For a conserved family the human member is
>=90 % identical to another organism's and is dropped, so ``--tier human_only``
reported ``no_eligible_instances`` for families that human demonstrably carries --
10 566 of 12 769 on a real request. ``Pfam-A.regions.tsv.gz`` has no redundancy
reduction, and the Swiss-Prot flat file supplies the sequence and taxon it lacks.

What is asserted here is the part that can silently drift:

* the five output files' formats, because ppi-splitting's own splitters and an
  embedding pipeline downstream both parse them positionally-by-name;
* ``instance_id`` construction and the 1-based inclusive slice, because embedding
  HDF5 keys are built from it and a one-off error is invisible;
* the human-identical invariant: a family's human instances under ``--tier any``
  are exactly the ones it keeps under ``human_only``;
* that a region whose parent is outside the universe is skipped, not guessed at;
* that filling a family prefers distinct parent proteins, which is the quantity
  SELECT_EXAMPLES and an external test set are short of.

Runs end to end with local fixtures and no network at all. Run directly
(``python3 tests/test_fetch_domains.py``) or via pytest.
"""

import csv
import gzip
import os
import subprocess
import sys
import tempfile
from collections import Counter

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)

FETCH = os.path.join(BIN, "fetch_domains.py")

from fetch_domains import (  # noqa: E402
    INSTANCE_COLUMNS,
    TIER_ELIGIBILITY,
    eligible_taxa,
    parse_uniprot_dat,
    tier_of,
)

ENV = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))

REGIONS_HEADER = [
    "pfamseq_acc", "seq_version", "crc64", "md5",
    "pfamA_acc", "seq_start", "seq_end", "ali_start", "ali_end",
]

# PF00001: one human protein with two regions, plus a second human protein -- so
# "prefer distinct parents" has something to prefer. PF00002: human + mouse, which
# is what separates human_only from any. PF00003: mouse only, so human_only must
# drop it with no_eligible_instances rather than fall back. PF00009: a region whose
# parent is not in the DAT at all.
REGIONS = [
    ("P11111", "PF00001", 10, 20, 11, 19),
    ("P11111", "PF00001", 50, 62, 51, 60),
    ("P22222", "PF00001", 5, 15, 6, 14),
    ("P33333", "PF00002", 1, 12, 2, 11),
    ("Q44444", "PF00002", 3, 14, 4, 13),
    ("Q55555", "PF00003", 7, 17, 8, 16),
    ("P99999", "PF00009", 1, 10, 1, 10),   # P99999 absent from the DAT
    ("P66666", "PF00004", 900, 999, 901, 990),  # coordinates past the sequence end
] + [
    # PF00005 is the repeat-domain shape that motivates the per-parent cap: one
    # protein carrying ten copies of the family, another carrying one. A plain
    # reservoir would fill the family almost entirely from P77777.
    ("P77777", "PF00005", i * 3 + 1, i * 3 + 3, i * 3 + 1, i * 3 + 3) for i in range(10)
] + [
    ("P88888", "PF00005", 1, 6, 1, 6),
]

# accession -> (taxon, sequence). Sequences are long enough for every slice above.
PROTEINS = {
    "P11111": ("9606", "".join(f"{i % 10}" for i in range(80)).replace("0", "A")),
    "P22222": ("9606", "M" * 40),
    "P33333": ("9606", "K" * 40),
    "Q44444": ("10090", "L" * 40),
    "Q55555": ("10090", "V" * 40),
    "P66666": ("9606", "W" * 40),
    "P77777": ("9606", "R" * 40),
    "P88888": ("9606", "S" * 40),
}

FAMILIES = ["PF00001", "PF00002", "PF00003", "PF00004", "PF00005", "PF00009"]


def write_fixtures(tmp):
    regions = os.path.join(tmp, "Pfam-A.regions.tsv.gz")
    with gzip.open(regions, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(REGIONS_HEADER)
        for acc, fam, s_start, s_end, a_start, a_end in REGIONS:
            w.writerow([acc, "1", "CRC", "MD5", fam, s_start, s_end, a_start, a_end])

    dat = os.path.join(tmp, "uniprot_sprot.dat.gz")
    with gzip.open(dat, "wt") as fh:
        for acc, (taxon, seq) in PROTEINS.items():
            fh.write(f"ID   {acc}_TEST   Reviewed;   {len(seq)} AA.\n")
            fh.write(f"AC   {acc};\n")
            fh.write(f"OX   NCBI_TaxID={taxon};\n")
            fh.write(f"SQ   SEQUENCE   {len(seq)} AA;  0 MW;  0 CRC64;\n")
            for i in range(0, len(seq), 60):
                fh.write(f"     {seq[i:i + 60]}\n")
            fh.write("//\n")

    clans = os.path.join(tmp, "Pfam-A.clans.tsv.gz")
    with gzip.open(clans, "wt") as fh:
        fh.write("PF00001\tCL0001\tname\tshort\tdesc\n")   # PF00002+ have no clan row

    families = os.path.join(tmp, "families.txt")
    with open(families, "w") as fh:
        fh.write("\n".join(FAMILIES) + "\n")
    return regions, dat, clans, families


def run(tmp, tier, pool_size=25):
    regions, dat, clans, families = write_fixtures(tmp)
    out = os.path.join(tmp, tier)
    os.makedirs(out, exist_ok=True)
    proc = subprocess.run(
        [sys.executable, FETCH,
         "--families", families,
         "--pool-size", str(pool_size),
         "--seed", "42",
         "--tier", tier,
         "--clans", clans,
         "--pfam-regions", regions,
         "--uniprot-dat", dat,
         "--pfam-release", "38.2"],
        cwd=out, env=ENV, capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise AssertionError(f"fetch_domains.py failed:\n{proc.stderr}")
    return out, proc.stderr


def read_instances(out):
    with open(os.path.join(out, "instances.tsv"), newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta(out):
    seqs, key = {}, None
    with open(os.path.join(out, "sequences.fasta")) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                key = line[1:]
            elif key:
                seqs[key] = seqs.get(key, "") + line
    return seqs


def test_unit_helpers():
    assert eligible_taxa(TIER_ELIGIBILITY["human_only"]) == {"9606"}
    assert eligible_taxa(TIER_ELIGIBILITY["any"]) is None
    # Taxon, not organism mnemonic, and reviewed-ness still splits the ladder.
    assert tier_of("9606", "P1", {"P1"}) == 0
    assert tier_of("9606", "P1", set()) == 1
    assert tier_of("10090", "P1", {"P1"}) == 2
    assert tier_of("10090", "P1", set()) == 3

    with tempfile.TemporaryDirectory() as tmp:
        _regions, dat, _clans, _families = write_fixtures(tmp)
        with gzip.open(dat, "rb") as fh:
            everything = parse_uniprot_dat(fh)
        assert set(everything) == set(PROTEINS), sorted(everything)
        assert everything["P11111"][0] == "9606"
        assert everything["P11111"][1] == PROTEINS["P11111"][1]
        with gzip.open(dat, "rb") as fh:
            human = parse_uniprot_dat(fh, {"9606"})
        assert set(human) == {a for a, (t, _s) in PROTEINS.items() if t == "9606"}, sorted(human)
        assert "Q44444" not in human and "Q55555" not in human, sorted(human)


def test_formats_and_coordinates():
    with tempfile.TemporaryDirectory() as tmp:
        out, log = run(tmp, "any")

        rows = read_instances(out)
        with open(os.path.join(out, "instances.tsv"), newline="") as fh:
            header = fh.readline().rstrip("\r\n").split("\t")
        assert header == list(INSTANCE_COLUMNS), header

        by_id = {r["instance_id"]: r for r in rows}
        # instance_id is {family}_{accession}_{ali_start}_{ali_end} ...
        assert "PF00001_P11111_11_19" in by_id, sorted(by_id)
        # ... and the sequence is the 1-based inclusive slice of the parent.
        parent = PROTEINS["P11111"][1]
        seqs = read_fasta(out)
        assert seqs["PF00001_P11111_11_19"] == parent[10:19]
        assert len(seqs["PF00001_P11111_11_19"]) == 9

        rec = by_id["PF00001_P11111_11_19"]
        assert (rec["family"], rec["protein_id"], rec["start"], rec["end"]) == ("PF00001", "P11111", "11", "19")
        assert rec["taxon_id"] == "9606"
        assert rec["source_db"] == "reviewed"
        assert rec["clan"] == "CL0001"
        # A family with no clan row gets the synthetic singleton, never a blank --
        # MAKE_METIS contracts the graph on this column.
        assert by_id["PF00002_P33333_2_11"]["clan"] == "CL_PF00002"

        # The mouse-only family survives under `any`, via tier 2.
        assert any(r["family"] == "PF00003" for r in rows), sorted(by_id)

        # Every FASTA key is an instance id and every instance has a sequence.
        assert set(seqs) == set(by_id), (sorted(set(seqs) ^ set(by_id)))

        # species.tsv: instance rows plus one row per family, whose taxon is the
        # modal taxon of its instances.
        with open(os.path.join(out, "species.tsv"), newline="") as fh:
            species = {r["protein_id"]: r["taxon_id"] for r in csv.DictReader(fh, delimiter="\t")}
        assert species["PF00001"] == "9606"
        assert species["PF00003"] == "10090"
        assert species["PF00001_P11111_11_19"] == "9606"

        # dropped_families.tsv is always written and its reasons are load-bearing
        # downstream. PF00009's only parent is outside the universe; PF00004's only
        # region is out of range.
        with open(os.path.join(out, "dropped_families.tsv"), newline="") as fh:
            dropped = {r["family"]: r["reason"] for r in csv.DictReader(fh, delimiter="\t")}
        assert set(dropped) == {"PF00004", "PF00009"}, dropped
        assert set(dropped.values()) == {"no_eligible_instances"}, dropped

        # go_annotations.tsv exists and is header-only: SUBSET_DOMAIN_DATA requires
        # the file, DDI mode has no family-level GO to put in it.
        with open(os.path.join(out, "go_annotations.tsv")) as fh:
            go = fh.read().splitlines()
        assert len(go) == 1 and go[0].split("\t")[0] == "protein_id", go

        assert "region rows skipped" in log, log
        assert "coordinates fall outside" in log, log


def test_human_only_is_a_subset_and_keeps_the_same_human_instances():
    with tempfile.TemporaryDirectory() as tmp:
        any_out, _ = run(tmp, "any")
        human_out, human_log = run(tmp, "human_only")

        any_rows = read_instances(any_out)
        any_ids = {r["instance_id"] for r in any_rows}
        human_ids = {r["instance_id"] for r in read_instances(human_out)}

        assert human_ids < any_ids, sorted(human_ids - any_ids)
        # The invariant that matters, and it is stronger than "subset": the human
        # instances a family keeps are byte-identical between the two tiers. Only
        # the non-human ones disappear -- both Q44444's PF00002 region, whose family
        # survives on its human sibling, and all of mouse-only PF00003.
        assert human_ids == {r["instance_id"] for r in any_rows if r["taxon_id"] == "9606"}
        assert all(r["taxon_id"] == "9606" for r in read_instances(human_out))

        # A family with only non-human regions is dropped rather than filled from
        # a non-human stratum -- that is what human_only means.
        with open(os.path.join(human_out, "dropped_families.tsv"), newline="") as fh:
            dropped = {r["family"]: r["reason"] for r in csv.DictReader(fh, delimiter="\t")}
        assert dropped["PF00003"] == "no_eligible_instances", dropped
        assert "9606" in human_log, human_log


def test_filling_prefers_distinct_parent_proteins():
    with tempfile.TemporaryDirectory() as tmp:
        # Room for two of PF00001's three regions. P11111 owns two of them, so a
        # uniform draw would sometimes take both and leave one parent; the
        # round-robin must take one from each parent instead.
        out, _ = run(tmp, "any", pool_size=2)
        parents = [r["protein_id"] for r in read_instances(out) if r["family"] == "PF00001"]
        assert len(parents) == 2, parents
        assert len(set(parents)) == 2, parents


def test_a_repeat_domain_family_is_not_filled_from_one_protein():
    with tempfile.TemporaryDirectory() as tmp:
        # pool_size 25 -> cap 5. P77777 carries ten PF00005 regions, P88888 one.
        out, log = run(tmp, "any", pool_size=25)
        rows = [r for r in read_instances(out) if r["family"] == "PF00005"]
        parents = Counter(r["protein_id"] for r in rows)
        assert parents["P77777"] == 5, parents
        assert parents["P88888"] == 1, parents
        assert "per-family cap of 5" in log, log

        # And the cap scales with pool_size rather than being a constant: at
        # pool_size 5 the cap is 1, so the family reduces to one region per parent.
        out, _ = run(tmp, "any", pool_size=5)
        rows = [r for r in read_instances(out) if r["family"] == "PF00005"]
        parents = Counter(r["protein_id"] for r in rows)
        assert parents == {"P77777": 1, "P88888": 1}, parents


if __name__ == "__main__":
    test_unit_helpers()
    test_formats_and_coordinates()
    test_human_only_is_a_subset_and_keeps_the_same_human_instances()
    test_filling_prefers_distinct_parent_proteins()
    test_a_repeat_domain_family_is_not_filled_from_one_protein()
    print("ok: regions-based sampling satisfies the output contract")
