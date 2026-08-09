#!/usr/bin/env python3
"""Fetch Pfam domain instances for all families in a DDI CSV (DDI mode's FETCH_DATA).

In DDI mode the interaction file's protein1/protein2 columns hold Pfam family
accessions rather than UniProt accessions. This script turns that family list
into the same four artefacts DATA_PREP produces in PPI mode -- sequences.fasta,
species.tsv, go_annotations.tsv -- plus instances.tsv, the single table that
maps every sampled domain instance back to its family, clan and parent protein.

Everything comes from four bulk downloads and one streaming pass, with no
per-family API requests: the per-family InterPro endpoint costs ~29 s, so the
full 30,134-family set would be ~242 sequential hours against a public service.

    Pfam-A.fasta.gz       every instance of every family: parent accession,
                          organism mnemonic, start-end, family, and the domain
                          sequence already cut (so we never cut one ourselves,
                          and sequence and coordinates come from one release)
    Pfam-A.clans.tsv.gz   family -> clan; a blank clan column means singleton
    speclist.txt          organism mnemonic -> numeric taxon id
    reviewed accessions   the Swiss-Prot accession set, for tier classification

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
import urllib.error
import urllib.request
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fetch_data import _ssl_context, write_go_tsv, write_species_tsv
from utils import read_ppis, write_fasta

PFAM_BASE = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
PFAM_FASTA_URL = f"{PFAM_BASE}/Pfam-A.fasta.gz"
PFAM_CLANS_URL = f"{PFAM_BASE}/Pfam-A.clans.tsv.gz"
PFAM_DEAD_URL = f"{PFAM_BASE}/Pfam-A.dead.gz"
PFAM_VERSION_URL = f"{PFAM_BASE}/Pfam.version.gz"
SPECLIST_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/speclist.txt"
REVIEWED_URL = "https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true&fields=accession&format=list"

# Bumped whenever the sampling or output format changes, so a cache written by
# an older version of this script is ignored rather than silently reused.
FORMAT_VERSION = 1

# >A0A067SRH6_GALM3/383-505 A0A067SRH6.1 PF26733.1;03009_C;
#  ^entry (acc_mnemonic)    ^accession   ^family
HEADER_RE = re.compile(
    r"^(?P<entry>\S+)/(?P<start>\d+)-(?P<end>\d+)\s+(?P<acc>[A-Za-z0-9]+)(?:\.\d+)?\s+(?P<family>PF\d+)"
)

# HUMAN E    9606: N=Homo sapiens
SPECLIST_RE = re.compile(r"^(?P<mnemonic>\S+)\s+\S+\s+(?P<taxon>\d+):")

DEAD_AC_RE = re.compile(r"^#=GF\s+AC\s+(PF\d+)")

# Disjoint strata, in fill order: human-reviewed -> human -> reviewed -> any.
# Disjoint rather than nested so an instance lands in exactly one reservoir and
# the cascade needs no deduplication at the family boundary.
TIERS = ("human_reviewed", "human", "reviewed", "any")
N_TIERS = len(TIERS)


# ---------------------------------------------------------------------------
# Downloads
# ---------------------------------------------------------------------------


def _open_url(url, retries=3, timeout=120):
    """Open a URL for streaming, retrying transient failures with 2**attempt backoff."""
    req = urllib.request.Request(url, headers={"Accept-Encoding": "identity"})
    for attempt in range(retries):
        try:
            return urllib.request.urlopen(req, timeout=timeout, context=_ssl_context())
        except (urllib.error.HTTPError, urllib.error.URLError) as exc:
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


def parse_speclist(text):
    """Return {organism_mnemonic: taxon_id} from UniProt's speclist.txt."""
    taxa = {}
    for line in text.splitlines():
        match = SPECLIST_RE.match(line)
        if match:
            taxa.setdefault(match["mnemonic"], match["taxon"])
    return taxa


def parse_reviewed(text):
    """Return the set of reviewed (Swiss-Prot) accessions."""
    return {line.strip() for line in text.splitlines() if line.strip()}


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


def tier_of(mnemonic, accession, reviewed):
    """Classify one instance into a disjoint tier index (0 = most preferred)."""
    is_human = mnemonic == "HUMAN"
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


def sample_instances(stream, wanted, pool_size, reviewed, seed):
    """One pass over Pfam-A.fasta, returning ({family: [record]}, stats).

    A record is a dict with instance_id, family, protein_id, mnemonic, start,
    end, tier and sequence.

    Reservoirs for every requested family are held for the whole pass rather
    than flushed at each family boundary. Family blocks are contiguous in
    practice (141 blocks for 141 families in an 84 MB sample), so flushing
    would bound live memory tighter -- but holding them costs at most
    |families| x 4 x pool_size records (~150 MB at full Pfam with M = 5) and
    makes the result correct even if the file ever stops being contiguous,
    rather than silently truncating a family to its first block.
    """
    reservoirs = {}
    stats = Counter()
    seen_families = set()
    last_family = None

    pending = None  # (reservoir, slot, record) awaiting its sequence lines
    seq_parts = []

    def finish():
        if pending is None:
            return
        res, slot, rec = pending
        rec["sequence"] = b"".join(seq_parts).decode("ascii", "replace")
        if not rec["sequence"]:
            stats["empty_sequence"] += 1
        res.buf[slot] = rec

    for raw in stream:
        if raw.startswith(b">"):
            finish()
            pending, seq_parts = None, []
            stats["records"] += 1

            match = HEADER_RE.match(raw[1:].decode("ascii", "replace").rstrip())
            if not match:
                stats["unparsed_headers"] += 1
                continue
            family = match["family"]
            if family not in wanted:
                continue

            if family != last_family:
                if family in seen_families:
                    stats["noncontiguous_blocks"] += 1
                seen_families.add(family)
                last_family = family

            entry = match["entry"]
            mnemonic = entry.rsplit("_", 1)[1] if "_" in entry else ""
            accession = match["acc"]
            tier = tier_of(mnemonic, accession, reviewed)
            stats[f"tier_{TIERS[tier]}_offered"] += 1

            key = (family, tier)
            res = reservoirs.get(key)
            if res is None:
                res = reservoirs[key] = Reservoir(pool_size, seed, family, TIERS[tier])
            slot = res.offer()
            if slot is None:
                continue

            start, end = match["start"], match["end"]
            pending = (
                res,
                slot,
                {
                    "instance_id": f"{family}_{accession}_{start}_{end}",
                    "family": family,
                    "protein_id": accession,
                    "mnemonic": mnemonic,
                    "start": start,
                    "end": end,
                    "tier": tier,
                },
            )
        elif pending is not None:
            seq_parts.append(raw.rstrip())
    finish()

    # Fill each family from tier 0 upward until pool_size. When a tier holds
    # more survivors than there is room for, sample rather than truncate --
    # taking the first few would bias towards whatever order the file happens
    # to be in within that tier.
    by_family = {}
    for family in sorted(seen_families):
        rng = random.Random(f"{seed}:{family}:pick")
        chosen, chosen_ids = [], set()
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
            picked = buf if len(buf) <= room else rng.sample(buf, room)
            chosen.extend(picked)
            chosen_ids.update(r["instance_id"] for r in picked)
        if chosen:
            by_family[family] = sorted(chosen, key=lambda r: r["instance_id"])
            for rec in by_family[family]:
                stats[f"tier_{TIERS[rec['tier']]}_kept"] += 1

    return by_family, stats


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------

INSTANCE_COLUMNS = ["instance_id", "family", "clan", "protein_id", "start", "end", "taxon_id", "source_db"]


def write_instances_tsv(path, by_family, clans, taxa, reviewed):
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
                        taxa.get(rec["mnemonic"], ""),
                        "reviewed" if rec["protein_id"] in reviewed else "unreviewed",
                    ]
                )


def read_instances_tsv(path):
    """Read back what write_instances_tsv wrote, as {family: [record]} (cache path)."""
    by_family = defaultdict(list)
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by_family[row["family"]].append(row)
    return dict(by_family)


def build_species(by_family, taxa):
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
            taxon = rec.get("taxon_id") if "taxon_id" in rec else taxa.get(rec["mnemonic"], "")
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


def derived_key(release, families, pool_size, seed):
    """Hash every input that can change the sampled instance set."""
    payload = "|".join([str(FORMAT_VERSION), release, str(pool_size), str(seed), ",".join(sorted(families))])
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
        help="prefix for the four outputs; default writes sequences.fasta, species.tsv, "
        "go_annotations.tsv and instances.tsv in the working directory",
    )
    parser.add_argument("--pool-size", type=int, default=5, help="instances sampled per family (M); default 5")
    parser.add_argument("--seed", type=int, default=42, help="sampling seed; default 42")
    parser.add_argument("--clans", help="precomputed Pfam-A.clans.tsv(.gz), skipping that download")
    parser.add_argument("--pfam-fasta", help="local Pfam-A.fasta.gz, skipping the 6.3 GB download")
    parser.add_argument("--interpro-cache", help="directory for cached downloads and sampled instances")
    parser.add_argument("--cache-fasta", action="store_true", help="also cache Pfam-A.fasta.gz (~6.3 GB)")
    parser.add_argument("--pfam-release", help="skip the Pfam.version lookup and use this release string")
    return parser.parse_args()


def main():
    args = parse_args()
    prefix = f"{args.out_prefix}." if args.out_prefix else ""
    fasta_out = f"{prefix}sequences.fasta"
    species_out = f"{prefix}species.tsv"
    go_out = f"{prefix}go_annotations.tsv"
    instances_out = f"{prefix}instances.tsv"

    families = read_families(args)
    print(f"Requested {len(families):,} Pfam families", file=sys.stderr)

    release = args.pfam_release or parse_release(download_text(PFAM_VERSION_URL, gzipped=True))
    print(f"Pfam release {release}", file=sys.stderr)
    cache = Cache(args.interpro_cache, release)

    if args.clans:
        opener = gzip.open if args.clans.endswith(".gz") else open
        with opener(args.clans, "rt") as fh:
            clans = parse_clans(fh.read())
    else:
        clans = parse_clans(cache.text("Pfam-A.clans.tsv", PFAM_CLANS_URL, gzipped=True))
    taxa = parse_speclist(cache.text("speclist.txt", SPECLIST_URL))
    reviewed = parse_reviewed(cache.text("reviewed.list", REVIEWED_URL))
    print(
        f"Reference tables: {len(clans):,} families, {len(taxa):,} organism mnemonics, "
        f"{len(reviewed):,} reviewed accessions",
        file=sys.stderr,
    )

    key = derived_key(release, families, args.pool_size, args.seed)
    cached_instances = cache.path(f"instances-{key}.tsv")
    cached_sequences = cache.path(f"sequences-{key}.fasta")
    stats = Counter()

    if cached_instances and os.path.exists(cached_instances) and os.path.exists(cached_sequences):
        print(f"Reusing cached instances for key {key}", file=sys.stderr)
        by_family = read_instances_tsv(cached_instances)
        shutil.copyfile(cached_instances, instances_out)
        shutil.copyfile(cached_sequences, fasta_out)
    else:
        wanted = set(families)
        fasta_path = args.pfam_fasta
        if not fasta_path and args.cache_fasta and cache.dir:
            fasta_path = cache.path("Pfam-A.fasta.gz")
            if not os.path.exists(fasta_path):
                print(f"Caching {PFAM_FASTA_URL} to {fasta_path} (~6.3 GB)...", file=sys.stderr)
                download_to_file(PFAM_FASTA_URL, fasta_path)

        if fasta_path:
            print(f"Streaming {fasta_path}...", file=sys.stderr)
            with open(fasta_path, "rb") as fh:
                by_family, stats = sample_instances(iter_gzip_lines(fh), wanted, args.pool_size, reviewed, args.seed)
        else:
            print(f"Streaming {PFAM_FASTA_URL} (~6.3 GB gz)...", file=sys.stderr)
            with _open_url(PFAM_FASTA_URL, timeout=600) as resp:
                by_family, stats = sample_instances(iter_gzip_lines(resp), wanted, args.pool_size, reviewed, args.seed)

        write_instances_tsv(instances_out, by_family, clans, taxa, reviewed)
        seqs = {rec["instance_id"]: rec["sequence"] for recs in by_family.values() for rec in recs}
        write_fasta(seqs, seqs.keys(), fasta_out)
        if cached_instances:
            shutil.copyfile(instances_out, cached_instances)
            shutil.copyfile(fasta_out, cached_sequences)

    ids, species = build_species(by_family, taxa)
    write_species_tsv(species_out, ids, species)
    # R9 drops the functional_relatedness attributes, not the --go_annotations
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

    if stats:
        print(f"Scanned {stats['records']:,} Pfam-A records", file=sys.stderr)
        for tier in TIERS:
            offered, kept = stats[f"tier_{tier}_offered"], stats[f"tier_{tier}_kept"]
            print(f"  tier {tier}: {offered:,} offered, {kept:,} kept", file=sys.stderr)
        for label, message in (
            ("unparsed_headers", "headers not matching the expected format"),
            ("empty_sequence", "sampled instances with an empty sequence"),
            ("noncontiguous_blocks", "families whose block was re-entered (file not contiguous)"),
        ):
            if stats[label]:
                print(f"Warning: {stats[label]:,} {message}", file=sys.stderr)

    missing = [f for f in families if f not in by_family]
    if missing:
        dead = parse_dead(cache.text("Pfam-A.dead", PFAM_DEAD_URL, gzipped=True))
        killed = sorted(f for f in missing if f in dead)
        unknown = sorted(f for f in missing if f not in dead)
        if killed:
            print(f"Warning: {len(killed)} requested families are dead in Pfam: {_sample(killed)}", file=sys.stderr)
        if unknown:
            print(
                f"Warning: {len(unknown)} requested families have no instances in Pfam-A.fasta "
                f"and are not listed as dead: {_sample(unknown)}",
                file=sys.stderr,
            )


if __name__ == "__main__":
    main()
