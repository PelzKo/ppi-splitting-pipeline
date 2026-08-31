#!/usr/bin/env python3
"""Fetch Pfam domain instances for all families in a DDI CSV (DDI mode's FETCH_DATA).

In DDI mode the interaction file's protein1/protein2 columns hold Pfam family
accessions rather than UniProt accessions. This script turns that family list
into the same four artefacts DATA_PREP produces in PPI mode -- sequences.fasta,
species.tsv, go_annotations.tsv -- plus instances.tsv, the single table that
maps every sampled domain instance back to its family, clan and parent protein,
and dropped_families.tsv, one row per requested family that kept no instance.

Everything comes from three bulk downloads and one streaming pass, with no
per-family API requests: the per-family InterPro endpoint costs ~29 s, so the
full 30,134-family set would be ~242 sequential hours against a public service.

    Pfam-A.regions.tsv.gz every Pfam region of every reference-proteome protein:
                          family, parent accession and the alignment coordinates
    uniprot_sprot.dat.gz  the Swiss-Prot flat file: parent sequence (SQ) and
                          taxon (OX) per accession, and -- being Swiss-Prot --
                          the reviewed accession set itself
    Pfam-A.clans.tsv.gz   family -> clan; a blank clan column means singleton

**Why regions and not Pfam-A.fasta.** Pfam-A.fasta is "90 % non-redundant" (Pfam's
own release notes). For a well-conserved family the human member is >=90 %
identical to some other organism's and is dropped, so the family keeps plenty of
records with no _HUMAN among them and --tier human_only correctly reports
`no_eligible_instances` for it. Measured over a 12,769-family human request:
10,566 families dropped that way, leaving 2,140. Pfam-A.regions.tsv.gz is the same
reference-proteome basis with no redundancy reduction, so every human region is
present. It carries no sequence and no taxon, which is why the Swiss-Prot flat
file joins it here: the domain is cut as sequence[ali_start-1:ali_end].

Coordinates are the *alignment* region (ali_start/ali_end), not the envelope
(seq_start/seq_end) -- the same region Pfam-A.fasta's deflines carried, so
instance ids stay comparable with what earlier releases of this script produced.

Instance ids use '_' as the separator, never '|':

    PF00069_P12345_10_250

make_metis.py, bias_analysis.py and plot_similarity_heatmap.py all undo
`makeblastdb -parse_seqids` with `s.split("|")[1] if "|" in s else s`. On a
pipe-delimited instance id that yields the bare parent accession, which is
absent from `lengths`, so every hit would be dropped and the similarity graph
would come out empty with no error at all. Parse the id back with
`rsplit("_", 3)`.
"""

import argparse
import csv
import gzip
import hashlib
import io
import os
import random
import re
import shutil
import sys
import time
import urllib.request
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fetch_data import _ssl_context, write_go_tsv, write_species_tsv
from utils import INSTANCE_COLUMNS, read_ppis, write_fasta

PFAM_BASE = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
PFAM_REGIONS_URL = f"{PFAM_BASE}/Pfam-A.regions.tsv.gz"
PFAM_CLANS_URL = f"{PFAM_BASE}/Pfam-A.clans.tsv.gz"
PFAM_DEAD_URL = f"{PFAM_BASE}/Pfam-A.dead.gz"
PFAM_VERSION_URL = f"{PFAM_BASE}/Pfam.version.gz"
SPROT_URL = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/"
    "complete/uniprot_sprot.dat.gz"
)

# Bumped whenever the sampling or output format changes, so a cache written by
# an older version of this script is ignored rather than silently reused. 2:
# instances now come from Pfam-A.regions.tsv.gz rather than the 90 %-non-redundant
# Pfam-A.fasta, so a v1 cache holds a *different* -- and for --tier human_only a
# much smaller -- instance set for the same key inputs.
FORMAT_VERSION = 2

# The columns Pfam-A.regions.tsv.gz is read by name, not by position: the file has
# gained columns before (ali_start/ali_end are newer than seq_start/seq_end) and a
# positional read would silently shift.
REGIONS_FAMILY_COL = "pfamA_acc"
REGIONS_PROTEIN_COL = "pfamseq_acc"
REGIONS_START_COL = "ali_start"
REGIONS_END_COL = "ali_end"
REGIONS_REQUIRED = (REGIONS_FAMILY_COL, REGIONS_PROTEIN_COL, REGIONS_START_COL, REGIONS_END_COL)

HUMAN_TAXON = "9606"

# At most pool_size // PARENT_CAP_DIVISOR regions of one parent protein are offered
# to a (family, tier) reservoir, so a family ends up spread over several parents.
# Without it a repeat-domain family fills every slot from one protein -- one human
# protein carries 30+ zf-C2H2 or WD40 regions -- and then no instance *pair* with
# distinct parents exists, which is what SELECT_EXAMPLES and the external test set
# both require. The reservoir itself cannot prevent this: it samples uniformly over
# regions, and those regions genuinely are mostly one protein's.
PARENT_CAP_DIVISOR = 5

# OX   NCBI_TaxID=9606;
DAT_TAXON_RE = re.compile(r"NCBI_TaxID=(\d+)")

DEAD_AC_RE = re.compile(r"^#=GF\s+AC\s+(PF\d+)")

# Disjoint strata, in fill order: human-reviewed -> human -> reviewed -> any.
# Disjoint rather than nested so an instance lands in exactly one reservoir and
# the cascade needs no deduplication at the family boundary.
TIERS = ("human_reviewed", "human", "reviewed", "any")
N_TIERS = len(TIERS)

# --tier restricts which strata may be filled at all, as opposed to which are
# preferred. Under "human_only" the two non-human strata are never offered a
# record, so a family whose instances are all non-human ends with *zero*
# instances rather than falling back to them -- which is the point: domainsplit
# ingests these instances into a human-only database and aborts on any other
# taxon. Because the cascade fills top-down with the same room and every
# reservoir carries its own seed, the human instances a family keeps under
# "human_only" are identical to the ones it keeps under "any".
TIER_ELIGIBILITY = {
    "any": frozenset(range(N_TIERS)),
    "human_only": frozenset({0, 1}),
}

# One row per requested family that reaches zero instances, with the reason.
# Published rather than made a MultiQC section: the fetch runs once per *run*
# under the synthetic "_shared" meta, so it could only ever be a run-wide
# section, and QC does not run at all under --split_only.
DROPPED_FAMILIES_COLUMNS = ("family", "reason")
DROP_REASONS = {
    "no_eligible_instances": "have Pfam records but none in a tier --tier allows",
    "dead": "are dead in Pfam",
    "not_in_pfam": "have no regions in Pfam-A.regions and are not listed as dead",
}


# ---------------------------------------------------------------------------
# Downloads
# ---------------------------------------------------------------------------


def _open_url(url, retries=3, timeout=120):
    """Open a URL for streaming, retrying transient failures with 2**attempt backoff.

    Catches OSError rather than the URLError pair: a connect or read timeout can
    surface as a bare TimeoutError (socket.timeout is an alias for it) without
    being wrapped, and gzip.BadGzipFile on a truncated small download is another
    OSError subclass. HTTPError and URLError are both OSError subclasses too, so
    this is a widening, not a change of intent. Note it covers the *open* only --
    a truncation part-way through the 6.3 GB stream raises out of the caller's
    iteration, which is what FETCH_DOMAIN_META's 'error_retry' label is for.
    """
    req = urllib.request.Request(url, headers={"Accept-Encoding": "identity"})
    for attempt in range(retries):
        try:
            return urllib.request.urlopen(req, timeout=timeout, context=_ssl_context())
        except OSError as exc:
            if attempt < retries - 1:
                time.sleep(2**attempt)
                continue
            raise RuntimeError(f"Failed to fetch {url} after {retries} attempts: {exc}") from exc
    raise RuntimeError(f"_open_url called with retries={retries} <= 0")


def download_text(url, gzipped=False, retries=3):
    """Return a small download as text, decompressing gzip in memory."""
    with _open_url(url, retries=retries) as resp:
        raw = resp.read()
    if gzipped:
        raw = gzip.decompress(raw)
    return raw.decode("utf-8", "replace")


def download_to_file(url, path, retries=3):
    """Stream a download to disk (used only for the cached copy of the big FASTA)."""
    tmp = f"{path}.part"
    with _open_url(url, retries=retries) as resp, open(tmp, "wb") as fh:
        shutil.copyfileobj(resp, fh, length=1 << 20)
    os.replace(tmp, path)
    return path


def iter_gzip_lines(fileobj):
    """Yield raw bytes lines from a gzip stream, without landing it on disk.

    gzip.GzipFile streams from any read()-able object and handles concatenated
    gzip members natively, so this works identically on an HTTP response and on
    a local file. Lines keep their trailing newline; callers rstrip what they use.
    """
    gz = gzip.GzipFile(fileobj=fileobj, mode="rb")
    return io.BufferedReader(gz, buffer_size=1 << 20)


# ---------------------------------------------------------------------------
# Reference tables
# ---------------------------------------------------------------------------


def parse_release(text):
    """Return the Pfam release string, e.g. '38.2', from Pfam.version's contents."""
    for line in text.splitlines():
        if line.lower().startswith("pfam release"):
            return line.split(":", 1)[1].strip()
    raise RuntimeError(f"Could not find a release line in Pfam.version: {text.strip()[:200]!r}")


def parse_clans(text):
    """Return {family: clan}. A blank clan column means the family is its own singleton clan."""
    clans = {}
    for row in csv.reader(io.StringIO(text), delimiter="\t"):
        if not row or not row[0].strip():
            continue
        fam = row[0].strip()
        clan = row[1].strip() if len(row) > 1 else ""
        clans[fam] = clan if clan else f"CL_{fam}"
    return clans


def parse_uniprot_dat(stream, taxa=None):
    """Return {accession: (taxon_id, sequence)} from a Swiss-Prot flat-file stream.

    This is the whole protein universe: Pfam-A.regions gives coordinates and
    nothing else, so sequence *and* taxon come from here, and membership in this
    mapping is what "reviewed" means -- Swiss-Prot is the reviewed set, which is
    why no separate reviewed-accession download is needed any more.

    `taxa` optionally restricts the universe to those NCBI taxon ids while
    parsing. Under --tier human_only the eligible strata are human-only, so the
    caller passes {"9606"} and the universe is ~20k entries instead of ~575k --
    the difference between ~10 MB and several hundred MB held for the whole
    streaming pass.

    Every accession of an entry is indexed, not only the primary one: Pfam's
    pfamseq_acc follows UniProt's primary accession, but an accession demoted to
    secondary between the two releases would otherwise silently lose all of its
    regions.
    """
    universe = {}
    accessions, taxon, seq_parts, in_seq = [], "", [], False

    def flush():
        if accessions and taxon and (taxa is None or taxon in taxa):
            sequence = "".join(seq_parts)
            for acc in accessions:
                universe[acc] = (taxon, sequence)

    for raw in stream:
        line = raw.decode("ascii", "replace") if isinstance(raw, bytes) else raw
        if line.startswith("//"):
            flush()
            accessions, taxon, seq_parts, in_seq = [], "", [], False
        elif in_seq:
            seq_parts.append(line.strip().replace(" ", ""))
        elif line.startswith("AC   "):
            accessions.extend(a.strip() for a in line[5:].split(";") if a.strip())
        elif line.startswith("OX   "):
            match = DAT_TAXON_RE.search(line)
            taxon = match.group(1) if match else ""
        elif line.startswith("SQ   "):
            in_seq = True
    flush()
    return universe


def eligible_taxa(eligible):
    """The taxon filter implied by a set of eligible tier indices, or None for no filter.

    Tiers 0 and 1 are the human strata (see TIERS), so a run whose eligible set
    holds only those cannot use a non-human protein at all and the universe can be
    narrowed while it is being parsed. Any other eligible set admits every taxon.
    """
    return {HUMAN_TAXON} if set(eligible) <= {0, 1} else None


def parse_dead(text):
    """Return the set of killed Pfam families, so a missing family reads as dead, not as a typo."""
    return {match.group(1) for line in text.splitlines() if (match := DEAD_AC_RE.match(line))}


# ---------------------------------------------------------------------------
# Reservoir sampling
# ---------------------------------------------------------------------------


class Reservoir:
    """Algorithm-R reservoir over one (family, tier) stratum.

    Every reservoir derives from the run's single --seed (params.seed); the
    family and tier only namespace it, so which instances a family contributes
    depends on that family's own records and nothing else -- not on how many
    families precede it in the stream, nor on how many other strata were
    offered a record in between. One global RNG would instead make a family's
    sample depend on stream position, so adding one family to the samplesheet
    would resample every family after it.

    random.Random() seeds a str via SHA-512, not hash(), so this is stable
    across processes and independent of PYTHONHASHSEED.
    """

    __slots__ = ("k", "rng", "buf", "seen")

    def __init__(self, k, seed, family, tier):
        self.k = k
        self.rng = random.Random(f"{seed}:{family}:{tier}")
        self.buf = []
        self.seen = 0

    def offer(self):
        """Reserve a slot for the next record, or return None if it is not sampled.

        Deciding before the record is materialised lets the caller skip
        accumulating the sequence for the ~99.99 % of instances that are
        rejected, which is where the streaming pass spends its time.
        """
        self.seen += 1
        if len(self.buf) < self.k:
            self.buf.append(None)
            return len(self.buf) - 1
        j = self.rng.randrange(self.seen)
        return j if j < self.k else None


def tier_of(taxon_id, accession, reviewed):
    """Classify one instance into a disjoint tier index (0 = most preferred).

    Keyed on the numeric taxon straight out of the flat file's OX line, replacing
    the organism mnemonic Pfam-A.fasta's deflines used to carry. `reviewed` stays a
    parameter rather than becoming `True`: with a Swiss-Prot universe every
    accession is reviewed and tiers 1 and 3 never fill, but the ladder is what
    makes adding a TrEMBL source later a data change instead of a redesign.
    """
    is_human = taxon_id == HUMAN_TAXON
    is_reviewed = accession in reviewed
    if is_human and is_reviewed:
        return 0
    if is_human:
        return 1
    if is_reviewed:
        return 2
    return 3


# ---------------------------------------------------------------------------
# The streaming pass
# ---------------------------------------------------------------------------


def sample_instances(stream, wanted, pool_size, universe, seed, eligible=None):
    """One pass over Pfam-A.regions.tsv, returning ({family: [record]}, stats, seen_families).

    A record is a dict with instance_id, family, protein_id, taxon_id, start, end,
    tier and sequence. `eligible` is a set of tier indices (see TIER_ELIGIBILITY);
    records in any other tier are counted and dropped, so a family can appear in
    `seen_families` and still be absent from the returned mapping. That is what
    separates "Pfam has no such family" from "Pfam has it but not in a tier we
    accept" in the drop report.

    At most `pool_size // PARENT_CAP_DIVISOR` regions of any one parent protein are
    offered per (family, tier), so a family is spread over parents rather than
    filled from whichever protein happens to carry the most copies of the domain.

    `universe` is {accession: (taxon_id, sequence)} from parse_uniprot_dat. A region
    whose protein is absent from it is counted and skipped -- at full scale that is
    most of the file, because Pfam covers every reference proteome while the universe
    is Swiss-Prot. Membership in `universe` is also what `reviewed` means, so the
    tier lookup reads its keys.

    Unlike the Pfam-A.fasta pass this replaces, a record is complete at its own row:
    the sequence is sliced out of the parent rather than accumulated over following
    lines, so there is no pending-record state to carry across rows.

    Reservoirs for every requested family are held for the whole pass rather than
    flushed at each family boundary. Family blocks are contiguous in the real file
    (it is sorted by pfamA_acc), so flushing would bound live memory tighter -- but
    holding them makes the result correct even if the file ever stops being
    contiguous, rather than silently truncating a family to its first block.

    That bound is linear in pool_size = ddi_examples_target x
    ddi_examples_pool_factor: the current defaults (N = 5, factor = 5, M = 25) cost
    a few hundred MB of records against FETCH_DOMAIN_META's memory, on top of the
    universe itself. Re-check that headroom before raising either param -- the
    product is what matters, and the target scales it exactly as the factor does.
    """
    if eligible is None:
        eligible = TIER_ELIGIBILITY["any"]
    reservoirs = {}
    parent_cap = max(1, pool_size // PARENT_CAP_DIVISOR)
    parent_offers = defaultdict(Counter)
    stats = Counter()
    seen_families = set()
    last_family = None
    columns = None

    for raw in stream:
        line = raw.decode("utf-8", "replace") if isinstance(raw, bytes) else raw
        fields = line.rstrip("\r\n").split("\t")

        if columns is None:
            # Read the header by name. A file without one is fatal rather than
            # guessed at: a positional read of the wrong release would put the
            # envelope coordinates where the alignment ones belong and every
            # instance id would shift silently.
            columns = {name: i for i, name in enumerate(f.strip() for f in fields)}
            missing = [c for c in REGIONS_REQUIRED if c not in columns]
            if missing:
                raise RuntimeError(
                    f"Pfam-A.regions is missing column(s) {', '.join(missing)}; "
                    f"header was {fields}"
                )
            continue

        stats["records"] += 1
        if len(fields) <= max(columns[c] for c in REGIONS_REQUIRED):
            stats["malformed_rows"] += 1
            continue

        family = fields[columns[REGIONS_FAMILY_COL]].strip()
        if family not in wanted:
            continue

        if family != last_family:
            if family in seen_families:
                stats["noncontiguous_blocks"] += 1
            seen_families.add(family)
            last_family = family

        accession = fields[columns[REGIONS_PROTEIN_COL]].strip()
        entry = universe.get(accession)
        if entry is None:
            stats["not_in_universe"] += 1
            continue
        taxon_id, parent = entry

        tier = tier_of(taxon_id, accession, universe)
        stats[f"tier_{TIERS[tier]}_offered"] += 1
        if tier not in eligible:
            stats["ineligible_tier_records"] += 1
            continue

        key = (family, tier)
        # Bound one parent's share before the reservoir sees the row -- see
        # PARENT_CAP_DIVISOR. The cap is on rows *offered*, so the reservoir stays
        # unbiased across parents; within a parent it keeps that parent's first
        # `cap` regions, which is a positional bias deliberately traded for the
        # guarantee that `cap` parents' worth of regions cannot be crowded out.
        offered = parent_offers[key]
        if offered[accession] >= parent_cap:
            stats["parent_cap_skipped"] += 1
            continue
        offered[accession] += 1

        res = reservoirs.get(key)
        if res is None:
            res = reservoirs[key] = Reservoir(pool_size, seed, family, TIERS[tier])
        slot = res.offer()
        if slot is None:
            continue

        start_s, end_s = fields[columns[REGIONS_START_COL]].strip(), fields[columns[REGIONS_END_COL]].strip()
        try:
            start_i, end_i = int(start_s), int(end_s)
        except ValueError:
            stats["malformed_rows"] += 1
            res.buf[slot] = None
            continue
        # 1-based inclusive, which is what instance_id, instances.tsv's start/end
        # and domainsplit's start_pos/end_pos all mean.
        if start_i < 1 or end_i > len(parent) or end_i < start_i:
            stats["coords_out_of_range"] += 1
            res.buf[slot] = None
            continue
        sequence = parent[start_i - 1 : end_i]
        if not sequence:
            stats["empty_sequence"] += 1
            res.buf[slot] = None
            continue

        res.buf[slot] = {
            "instance_id": f"{family}_{accession}_{start_s}_{end_s}",
            "family": family,
            "protein_id": accession,
            "taxon_id": taxon_id,
            "start": start_s,
            "end": end_s,
            "tier": tier,
            "sequence": sequence,
        }

    if columns is None:
        raise RuntimeError("Pfam-A.regions stream was empty -- not even a header row")

    # Fill each family from tier 0 upward until pool_size. When a tier holds
    # more survivors than there is room for, sample rather than truncate --
    # taking the first few would bias towards whatever order the file happens
    # to be in within that tier.
    #
    # Parent proteins are preferred over extra regions of one protein: Pfam
    # routinely gives one protein several regions of the same family, and a
    # family filled from three proteins cannot supply an instance pair with
    # distinct parents, which is what the external test set and the negative
    # samplers both need. Within a parent the order is still the rng's.
    by_family = {}
    for family in sorted(seen_families):
        rng = random.Random(f"{seed}:{family}:pick")
        chosen, chosen_ids = [], set()
        # An ineligible tier simply has no reservoir -- the gate above never
        # created one -- so this loop needs no knowledge of `eligible`, and the
        # rng is consumed in the same order either way.
        for tier in range(N_TIERS):
            room = pool_size - len(chosen)
            if room <= 0:
                break
            res = reservoirs.get((family, tier))
            if res is None:
                continue
            buf = [r for r in res.buf if r and r["sequence"] and r["instance_id"] not in chosen_ids]
            if not buf:
                continue
            picked = buf if len(buf) <= room else _spread_over_parents(rng, buf, room)
            chosen.extend(picked)
            chosen_ids.update(r["instance_id"] for r in picked)
        if chosen:
            by_family[family] = sorted(chosen, key=lambda r: r["instance_id"])
            for rec in by_family[family]:
                stats[f"tier_{TIERS[rec['tier']]}_kept"] += 1

    return by_family, stats, seen_families


def _spread_over_parents(rng, records, room):
    """Pick `room` records, taking at most one per parent protein before repeating.

    Round-robin over shuffled per-parent queues. With one region per protein this
    is exactly the uniform sample it replaces; with several it maximises the number
    of distinct parents, which is the quantity SELECT_EXAMPLES and the external test
    set are actually short of.
    """
    by_parent = defaultdict(list)
    for rec in records:
        by_parent[rec["protein_id"]].append(rec)
    queues = [by_parent[p] for p in sorted(by_parent)]
    for queue in queues:
        rng.shuffle(queue)
    rng.shuffle(queues)

    picked = []
    while len(picked) < room:
        progressed = False
        for queue in queues:
            if not queue:
                continue
            picked.append(queue.pop())
            progressed = True
            if len(picked) == room:
                break
        if not progressed:
            break
    return picked


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------


def write_instances_tsv(path, by_family, clans, reviewed):
    """Write the one table MAKE_METIS, the splitters and SELECT_EXAMPLES all read.

    MAKE_METIS takes instance_id -> clan, the splitters invert to
    family -> {instance_id} and read family -> clan, SELECT_EXAMPLES takes
    protein_id. One table cannot contradict itself the way two would.
    """
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(INSTANCE_COLUMNS)
        for family in sorted(by_family):
            for rec in by_family[family]:
                writer.writerow(
                    [
                        rec["instance_id"],
                        family,
                        clans.get(family, f"CL_{family}"),
                        rec["protein_id"],
                        rec["start"],
                        rec["end"],
                        rec["taxon_id"],
                        "reviewed" if rec["protein_id"] in reviewed else "unreviewed",
                    ]
                )


def classify_dropped(families, by_family, seen_families, dead):
    """Return [(family, reason)] for every requested family that kept no instance.

    Three distinct causes, and conflating them would hide the one that is a
    configuration choice rather than a fact about Pfam: a family absent from the
    stream is dead or mistyped, while a family *present* in the stream with
    nothing in an eligible tier is --tier doing exactly what it was asked to.
    """
    dropped = []
    for family in sorted(families):
        if family in by_family:
            continue
        if family in seen_families:
            dropped.append((family, "no_eligible_instances"))
        elif family in dead:
            dropped.append((family, "dead"))
        else:
            dropped.append((family, "not_in_pfam"))
    return dropped


def write_dropped_families(path, dropped):
    """Write the drop report. Always written, header-only when nothing dropped.

    FETCH_DOMAIN_META declares it as a required output, so an unconditional file
    keeps the task from failing on the happy path.
    """
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(DROPPED_FAMILIES_COLUMNS)
        writer.writerows(dropped)


def read_dropped_families(path):
    """Read back what write_dropped_families wrote (cache path)."""
    with open(path) as fh:
        return [(row["family"], row["reason"]) for row in csv.DictReader(fh, delimiter="\t")]


def read_instances_tsv(path):
    """Read back what write_instances_tsv wrote, as {family: [record]} (cache path)."""
    by_family = defaultdict(list)
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by_family[row["family"]].append(row)
    return dict(by_family)


def build_species(by_family):
    """Return (ordered_ids, {id: taxon_id}) covering both instances and families.

    Instance rows serve QC's instance-level same_species; family rows serve the
    family-level taxon-pair term in sample_negatives_ilp.py. A family's taxon is
    the modal taxon of its sampled parents. One table, two consumers, no code
    change in either.
    """
    ids, species = [], {}
    for family in sorted(by_family):
        counts = Counter()
        for rec in by_family[family]:
            # Present on every record, whether it came from the streaming pass or
            # was read back out of the cache's instances.tsv.
            taxon = rec.get("taxon_id") or ""
            ids.append(rec["instance_id"])
            species[rec["instance_id"]] = taxon
            if taxon:
                counts[taxon] += 1
        ids.append(family)
        species[family] = counts.most_common(1)[0][0] if counts else ""
    return ids, species


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------


class Cache:
    """Release-keyed on-disk cache. A convenience, never a dependency.

    An absent, partial or stale cache falls through to the full stream and
    produces byte-identical output -- affordable precisely because a cold run
    is one gzip stream, i.e. minutes rather than the 242 h the per-family API
    alternative would have cost.
    """

    def __init__(self, root, release):
        self.dir = os.path.join(root, f"pfam-{release}") if root else None
        if self.dir:
            os.makedirs(self.dir, exist_ok=True)

    def path(self, name):
        return os.path.join(self.dir, name) if self.dir else None

    def text(self, name, url, gzipped=False):
        """Return a small reference table, reading from cache when present."""
        path = self.path(name)
        if path and os.path.exists(path):
            with open(path) as fh:
                return fh.read()
        text = download_text(url, gzipped=gzipped)
        if path:
            tmp = f"{path}.part"
            with open(tmp, "w") as fh:
                fh.write(text)
            os.replace(tmp, path)
        return text


def content_digest(items):
    """Stable digest of a set of ids, for the cache key. ~0.3 s over 570k accessions."""
    return hashlib.sha256("\n".join(sorted(items)).encode()).hexdigest()


def derived_key(release, families, pool_size, seed, universe_digest, tier):
    """Hash every input that can change the cached instance set or its columns.

    `tier` is in here for the obvious reason and one less obvious one: it decides
    which strata may be filled *and* how wide the protein universe is parsed
    (eligible_taxa), so a warm cache written under --tier any would otherwise serve
    a --tier human_only run non-human instances, which downstream (domainsplit's
    human-only ingest) is a hard failure several steps later rather than a wrong
    number here.

    `universe_digest` covers the Swiss-Prot flat file and belongs in here as much as
    the Pfam release does. It is now the *whole* join partner -- coordinates come
    from Pfam, but membership, taxon and sequence all come from the DAT -- so a new
    UniProt release, roughly every eight weeks, changes the sample and the sequences
    at the same --seed against the same Pfam release. Digesting its accession set
    also settles the Cache's pinning: the DAT is filed under pfam-{release}/, where
    the release string alone would let two machines with differently-aged copies
    disagree silently under one key.
    """
    payload = "|".join(
        [
            str(FORMAT_VERSION),
            release,
            str(pool_size),
            str(seed),
            tier,
            universe_digest,
            ",".join(sorted(families)),
        ]
    )
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _sample(items, limit=20):
    """Render a list for a log line, truncated -- a full-Pfam run can report
    thousands of families and the count is the actionable part, not the tail."""
    items = sorted(items)
    if len(items) <= limit:
        return ", ".join(items)
    return f"{', '.join(items[:limit])}, ... and {len(items) - limit} more"


def read_families(args):
    """Return the sorted set of requested Pfam families from --domains or --families."""
    if args.domains:
        rows = read_ppis(args.domains)
        return sorted({p for row in rows for p in (row["protein1"], row["protein2"]) if p})
    with open(args.families) as fh:
        return sorted({line.strip() for line in fh if line.strip()})


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--domains",
        help="interaction CSV whose protein1/protein2 columns hold Pfam family accessions",
    )
    source.add_argument("--families", help="plain text file, one Pfam family accession per line")
    parser.add_argument(
        "--out-prefix",
        default="",
        help="prefix for the outputs; default writes sequences.fasta, species.tsv, "
        "go_annotations.tsv, instances.tsv and dropped_families.tsv in the working directory",
    )
    parser.add_argument("--pool-size", type=int, default=5, help="instances sampled per family (M); default 5")
    parser.add_argument("--seed", type=int, default=42, help="sampling seed; default 42")
    parser.add_argument(
        "--tier",
        choices=sorted(TIER_ELIGIBILITY),
        default="any",
        help="which strata may be sampled: 'any' fills all four in preference order, "
        "'human_only' makes the two non-human strata ineligible, so a family with no "
        "human instance ends with zero instances; default any",
    )
    parser.add_argument("--clans", help="precomputed Pfam-A.clans.tsv(.gz), skipping that download")
    parser.add_argument("--pfam-regions", help="local Pfam-A.regions.tsv.gz, skipping the 4.7 GB download")
    parser.add_argument(
        "--uniprot-dat",
        help="local uniprot_sprot.dat.gz -- the protein universe: sequence, taxon and, by "
        "being Swiss-Prot, the reviewed set. Downloaded (~600 MB) when not given",
    )
    parser.add_argument("--interpro-cache", help="directory for cached downloads and sampled instances")
    parser.add_argument("--cache-regions", action="store_true", help="also cache Pfam-A.regions.tsv.gz (~4.7 GB)")
    parser.add_argument("--pfam-release", help="skip the Pfam.version lookup and use this release string")
    return parser.parse_args()


def main():
    args = parse_args()
    prefix = f"{args.out_prefix}." if args.out_prefix else ""
    fasta_out = f"{prefix}sequences.fasta"
    species_out = f"{prefix}species.tsv"
    go_out = f"{prefix}go_annotations.tsv"
    instances_out = f"{prefix}instances.tsv"
    dropped_out = f"{prefix}dropped_families.tsv"

    families = read_families(args)
    eligible = TIER_ELIGIBILITY[args.tier]
    print(
        f"Requested {len(families):,} Pfam families "
        f"(--tier {args.tier}: {', '.join(TIERS[t] for t in sorted(eligible))})",
        file=sys.stderr,
    )

    release = args.pfam_release or parse_release(download_text(PFAM_VERSION_URL, gzipped=True))
    print(f"Pfam release {release}", file=sys.stderr)
    cache = Cache(args.interpro_cache, release)

    if args.clans:
        opener = gzip.open if args.clans.endswith(".gz") else open
        with opener(args.clans, "rt") as fh:
            clans = parse_clans(fh.read())
    else:
        clans = parse_clans(cache.text("Pfam-A.clans.tsv", PFAM_CLANS_URL, gzipped=True))
    # The protein universe, narrowed to the taxa this tier can use while it is
    # parsed -- see eligible_taxa(). Parsed before the cache key is computed
    # because its accession set *is* the key's UniProt component.
    taxa_filter = eligible_taxa(eligible)
    dat_path = args.uniprot_dat or cache.path("uniprot_sprot.dat.gz")
    if not args.uniprot_dat:
        if not dat_path:
            raise RuntimeError(
                "--uniprot-dat is required when --interpro-cache is not set: the flat "
                "file is the protein universe and there is nowhere to cache it"
            )
        if not os.path.exists(dat_path):
            print(f"Caching {SPROT_URL} to {dat_path} (~600 MB)...", file=sys.stderr)
            download_to_file(SPROT_URL, dat_path)
    print(f"Reading the protein universe from {dat_path}...", file=sys.stderr)
    with open(dat_path, "rb") as fh:
        universe = parse_uniprot_dat(iter_gzip_lines(fh), taxa_filter)
    if not universe:
        raise RuntimeError(
            f"{dat_path} yielded no proteins"
            + (f" for taxon(s) {', '.join(sorted(taxa_filter))}" if taxa_filter else "")
            + " -- every region would be dropped as not_in_universe"
        )
    # Swiss-Prot *is* the reviewed set, so the universe doubles as it.
    reviewed = universe
    print(
        f"Reference tables: {len(clans):,} families, {len(universe):,} proteins"
        + (f" (taxon {', '.join(sorted(taxa_filter))})" if taxa_filter else " (all taxa)"),
        file=sys.stderr,
    )

    key = derived_key(
        release,
        families,
        args.pool_size,
        args.seed,
        content_digest(universe),
        args.tier,
    )
    cached_instances = cache.path(f"instances-{key}.tsv")
    cached_sequences = cache.path(f"sequences-{key}.fasta")
    cached_dropped = cache.path(f"dropped-{key}.tsv")
    stats = Counter()

    # All three or none: the drop report cannot be rebuilt from instances.tsv
    # (it names families that are *not* in it, and distinguishes a tier drop from
    # a dead accession), so a pair cached by an older version reads as a miss.
    if cached_instances and all(os.path.exists(p) for p in (cached_instances, cached_sequences, cached_dropped)):
        print(f"Reusing cached instances for key {key}", file=sys.stderr)
        by_family = read_instances_tsv(cached_instances)
        dropped = read_dropped_families(cached_dropped)
        shutil.copyfile(cached_instances, instances_out)
        shutil.copyfile(cached_sequences, fasta_out)
        shutil.copyfile(cached_dropped, dropped_out)
    else:
        wanted = set(families)
        regions_path = args.pfam_regions
        if not regions_path and args.cache_regions and cache.dir:
            regions_path = cache.path("Pfam-A.regions.tsv.gz")
            if not os.path.exists(regions_path):
                print(f"Caching {PFAM_REGIONS_URL} to {regions_path} (~4.7 GB)...", file=sys.stderr)
                download_to_file(PFAM_REGIONS_URL, regions_path)

        if regions_path:
            print(f"Streaming {regions_path}...", file=sys.stderr)
            with open(regions_path, "rb") as fh:
                by_family, stats, seen = sample_instances(
                    iter_gzip_lines(fh), wanted, args.pool_size, universe, args.seed, eligible
                )
        else:
            print(f"Streaming {PFAM_REGIONS_URL} (~4.7 GB gz)...", file=sys.stderr)
            with _open_url(PFAM_REGIONS_URL, timeout=600) as resp:
                by_family, stats, seen = sample_instances(
                    iter_gzip_lines(resp), wanted, args.pool_size, universe, args.seed, eligible
                )

        write_instances_tsv(instances_out, by_family, clans, reviewed)
        seqs = {rec["instance_id"]: rec["sequence"] for recs in by_family.values() for rec in recs}
        write_fasta(seqs, seqs.keys(), fasta_out)

        # The dead table is only downloaded when something actually dropped, so
        # the happy path still makes no extra request.
        absent = [f for f in families if f not in by_family and f not in seen]
        dead = parse_dead(cache.text("Pfam-A.dead", PFAM_DEAD_URL, gzipped=True)) if absent else set()
        dropped = classify_dropped(families, by_family, seen, dead)
        write_dropped_families(dropped_out, dropped)

        if cached_instances:
            shutil.copyfile(instances_out, cached_instances)
            shutil.copyfile(fasta_out, cached_sequences)
            shutil.copyfile(dropped_out, cached_dropped)

    ids, species = build_species(by_family)
    write_species_tsv(species_out, ids, species)
    # DDI mode drops the functional_relatedness attributes, not the --go_annotations
    # argument, so every consumer that requires the flag still gets a valid file.
    write_go_tsv(go_out, [], {})

    n_instances = sum(len(v) for v in by_family.values())
    print(
        f"Written {n_instances:,} instances across {len(by_family):,} families to {instances_out}",
        file=sys.stderr,
    )
    print(f"Written {n_instances:,} domain sequences to {fasta_out}", file=sys.stderr)
    print(f"Written {len(ids):,} species rows (instances + families) to {species_out}", file=sys.stderr)
    print(f"Written empty GO annotations to {go_out} (dropped in DDI mode)", file=sys.stderr)
    print(f"Written {len(dropped):,} dropped families to {dropped_out}", file=sys.stderr)

    if stats:
        print(f"Scanned {stats['records']:,} Pfam-A region rows", file=sys.stderr)
        for index, tier in enumerate(TIERS):
            offered, kept = stats[f"tier_{tier}_offered"], stats[f"tier_{tier}_kept"]
            gate = "" if index in eligible else "  [ineligible under --tier]"
            print(f"  tier {tier}: {offered:,} offered, {kept:,} kept{gate}", file=sys.stderr)
        if stats["parent_cap_skipped"]:
            print(
                f"  {stats['parent_cap_skipped']:,} region rows skipped: parent already at "
                f"its per-family cap of {max(1, args.pool_size // PARENT_CAP_DIVISOR)}",
                file=sys.stderr,
            )
        if stats["not_in_universe"]:
            print(
                f"  {stats['not_in_universe']:,} region rows skipped: parent protein outside "
                "the universe (expected -- Pfam covers every reference proteome)",
                file=sys.stderr,
            )
        for label, message in (
            ("malformed_rows", "region rows with missing or non-numeric fields"),
            ("coords_out_of_range", "regions whose coordinates fall outside the parent sequence"),
            ("empty_sequence", "sampled instances with an empty sequence"),
            ("noncontiguous_blocks", "families whose block was re-entered (file not contiguous)"),
        ):
            if stats[label]:
                print(f"Warning: {stats[label]:,} {message}", file=sys.stderr)

    # One warning per reason, so "--tier dropped 4,000 families" cannot hide
    # among dead accessions -- the first is a knob, the second is Pfam's history.
    by_reason = defaultdict(list)
    for family, reason in dropped:
        by_reason[reason].append(family)
    for reason, message in DROP_REASONS.items():
        hit = by_reason.get(reason)
        if hit:
            print(f"Warning: {len(hit)} requested families {message}: {_sample(hit)}", file=sys.stderr)


if __name__ == "__main__":
    main()
