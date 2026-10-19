"""
Primer Specificity Module (BLAST-based in-silico PCR)
-----------------------------------------------------
Checks primer pair specificity by BLASTing primers against the genome
and identifying all amplicons that would be produced.

Primer FASTA IDs must use "|" as delimiter:
    >SSR{ssr_id}|L   (left primer)
    >SSR{ssr_id}|R   (right primer)

This is produced automatically by primers_to_blast_fasta() in primer_design.py.
"""

import os
import subprocess
import sys
import tempfile
import hashlib
from dataclasses import dataclass
from typing import List, Dict, Any, Optional


# ---------------------------------------------------------------------------
# Bundled BLAST binary resolution
# ---------------------------------------------------------------------------

def get_blast_bin_dir() -> str:
    """
    Return the directory containing bundled BLAST binaries.

    Priority order:
      1. PyInstaller frozen bundle  → tries multiple known locations
      2. Development (source)       → <project_root>/blast/<platform>/
      3. System PATH fallback       → empty string (let OS resolve)
    """
    if sys.platform == "win32":
        subfolder = "windows"
    elif sys.platform == "darwin":
        subfolder = "macos"
    else:
        subfolder = "linux"

    if getattr(sys, "frozen", False):
        # PyInstaller 6.x puts bundled files in _internal/ next to the .exe
        # sys._MEIPASS points to _internal/ but let's check both locations
        candidates = [
            # Standard _MEIPASS location
            os.path.join(sys._MEIPASS, "blast", subfolder),
            # Relative to the exe itself (PyInstaller 6.x onedir)
            os.path.join(os.path.dirname(sys.executable), "_internal", "blast", subfolder),
            # Directly next to exe (older PyInstaller behaviour)
            os.path.join(os.path.dirname(sys.executable), "blast", subfolder),
        ]
    else:
        # Running from source — go up to project root from core/
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        candidates = [
            os.path.join(project_root, "blast", subfolder),
        ]

    for path in candidates:
        if os.path.isdir(path):
            return path

    # Fall back to system PATH
    return ""


# ---------------------------------------------------------------------------
# Parameter dataclasses
# ---------------------------------------------------------------------------

@dataclass
class BlastParams:
    # word_size, task, dust, soft_masking, penalty/reward/gapopen are
    # standard settings for short-sequence/primer BLAST screening.
    # evalue is set to BLAST's own native default (10) — a stricter
    # value trades away real sensitivity: it can suppress genuine
    # off-target hits from ever being reported at all, letting primers
    # pass specificity screening that a less restrictive search would
    # correctly have failed. Runs can take considerably longer on large,
    # repeat-dense genomes as a result — accepted deliberately in favour
    # of correctness. BLAST's own internal memory during the search
    # still scales with how permissive this value is, so a very large
    # genome remains something to watch for memory pressure.
    task: str = "blastn-short"
    word_size: int = 7
    evalue: float = 10.0
    reward: int = 1
    penalty: int = -3
    gapopen: int = 5
    gapextend: int = 2
    dust: str = "no"
    soft_masking: str = "false"


@dataclass
class SpecificityParams:
    # Defaults below follow standard primer-specificity screening
    # stringency.
    min_product: int = 100
    max_product: int = 300

    # A target is ignored only once it has 6+ total mismatches to a
    # primer — i.e. up to 5 are tolerated. min_identity is set low
    # enough (70%) that it doesn't become the accidental binding
    # constraint ahead of max_mismatches: 5 mismatches on an 18bp primer
    # (our min_align_len) is 72% identity, so a stricter min_identity
    # here would silently veto exactly the hits max_mismatches is meant
    # to still evaluate. The 3'-end check below is what does the real
    # precision work in this scheme — these whole-primer thresholds are
    # deliberately a coarse net, not the primary filter.
    max_mismatches: int = 5
    min_align_len: int = 18
    min_identity: float = 70.0

    allow_offtarget_amplicons: bool = True
    max_offtarget_amplicon_size: int = 1500
    require_unique_hits: bool = True
    # require_unique_hits: when True, adds an ADDITIONAL gate requiring each
    # primer to map to at most max_total_hits genomic locations (checked at
    # the same mismatch/identity thresholds above), on top of the amplicon
    # check. This does not replace the amplicon check — it's a stricter
    # layer on top of it.
    max_total_hits: int = 1

    require_3prime_match: bool = True
    three_prime_window: int = 5
    three_prime_min_mismatches_to_ignore: int = 2
    # require_3prime_match: whole-primer identity/mismatch counts alone
    # can't tell you WHERE a mismatch sits — and that's what actually
    # determines whether a hit is a real PCR competitor, since a mismatch
    # at or near the 3' end is what stops polymerase extension, while a
    # mismatch in the middle or 5' end of the primer often doesn't matter.
    #
    # By default, a hit is only considered safe to ignore once it has AT LEAST
    # three_prime_min_mismatches_to_ignore (2) mismatches within the last
    # three_prime_window (5) bases. A hit with 0 or 1 mismatch in that
    # terminal window is still treated as a genuine competitor — this is
    # NOT the same as requiring an exact match; a single terminal
    # mismatch alone does not clear a hit.


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _run_cmd(cmd: List[str], cwd: Optional[str] = None) -> str:
    """Run a shell command, raise RuntimeError on failure, return stdout."""
    kwargs = {}
    if sys.platform == "win32":
        kwargs["creationflags"] = subprocess.CREATE_NO_WINDOW

    result = subprocess.run(
        cmd,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        **kwargs,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed: {' '.join(cmd)}\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )
    return result.stdout


def _write_primers_fasta(primer_results: List[Dict[str, Any]], path: str) -> None:
    """
    Write primers to FASTA using "|" delimiter for safe ID parsing.
    Format: >SSR{ssr_id}_P{pair_rank}|L  and  >SSR{ssr_id}_P{pair_rank}|R

    pair_rank is included because a single SSR locus can have multiple
    candidate primer pairs (rank 0, 1, 2, ...) with genuinely different
    sequences. Without pair_rank in the query name, every rank's primers
    would share the identical query ID, causing BLAST to report hits
    against one merged, meaningless identity rather than each pair's own
    actual specificity.
    """
    with open(path, "w") as fh:
        for i, p in enumerate(primer_results):
            ssr_id    = p.get("ssr_id", f"pair{i}")
            pair_rank = p.get("pair_rank", 0)
            fh.write(f">SSR{ssr_id}_P{pair_rank}|L\n{p['left_primer']}\n")
            fh.write(f">SSR{ssr_id}_P{pair_rank}|R\n{p['right_primer']}\n")


def _parse_blast_tab(path: str, specificity: "SpecificityParams" = None) -> List[Dict[str, Any]]:
    """
    Parse BLAST tabular output (outfmt 6 + sstrand, qseq, sseq).

    If `specificity` is given, hits are filtered against its cheap
    numeric criteria (min_align_len, max_mismatches, min_identity) DURING
    parsing, before qseq/sseq are even extracted — on a repeat-dense
    genome, BLAST's own permissive search settings can produce hundreds
    of thousands of very weak candidate hits per primer, and every one
    of those that fails this same numeric check downstream in _good()
    (the first, cheap half of it) never has its qseq/sseq actually used
    for anything — those two alignment strings are the largest part of
    each hit's memory footprint, so holding all of them for hits that
    were always going to be discarded is a real, avoidable memory cost
    for large output volumes. This is a pure performance change: every
    caller that uses the parsed hits already only sees hits that pass
    this exact same filter, so nothing about GRACE's visible output
    (including the raw-hits export) changes as a result.
    """
    hits = []
    if not os.path.exists(path):
        return hits
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) != 15:
                continue  # skip malformed lines

            (
                qseqid, sseqid, pident, length,
                mismatch, gapopen, qstart, qend,
                sstart, send, evalue, bitscore, sstrand,
                qseq, sseq,
            ) = cols

            pident_v, length_v, mismatch_v = float(pident), int(length), int(mismatch)

            if specificity is not None:
                if (length_v   < specificity.min_align_len
                        or mismatch_v > specificity.max_mismatches
                        or pident_v   < specificity.min_identity):
                    continue  # would fail _good()'s base_ok check anyway — skip entirely

            hits.append({
                "qseqid": qseqid,
                "sseqid": sseqid,
                "pident": pident_v,
                "length": length_v,
                "mismatch": mismatch_v,
                "gapopen": int(gapopen),
                "qstart": int(qstart),
                "qend": int(qend),
                "sstart": int(sstart),
                "send": int(send),
                "evalue": float(evalue),
                "bitscore": float(bitscore),
                "sstrand": sstrand.strip(),
                "qseq": qseq.strip(),
                "sseq": sseq.strip(),
            })
    return hits


def _pair_hits(
    hits: List[Dict[str, Any]],
    specificity: SpecificityParams,
    query_lengths: Dict[str, int] = None,
    hit_writer=None,
) -> Dict[tuple, Dict[str, Any]]:
    """
    Pair left and right primer hits into amplicons and assess specificity.
    Returns a dict keyed by (ssr_id, pair_rank) — NOT ssr_id alone, since a
    single locus can have multiple candidate primer pairs (different
    ranks) with genuinely different specificity outcomes.

    query_lengths: {query_name: full primer length}, needed to check
        whether a hit's alignment actually reaches the primer's true 3'
        end (see _three_prime_matches). If not supplied, the 3'-end
        check is skipped (fails open) rather than incorrectly rejecting
        every hit.

    hit_writer: optional callable(ssr_id, pair_rank, side, hit_dict),
        called once for every surviving (post-_good()) hit. This is how
        raw hit-level detail is written straight to disk incrementally,
        batch by batch, rather than being accumulated into this
        function's own returned dict — for a large primer set against a
        repeat-dense genome, holding every surviving hit's full record
        (including complete alignment strings) in memory for the whole
        run, only to serialise it all at the very end, was a real,
        avoidable memory cost. The returned dict never carries the raw
        hit records themselves, regardless of whether hit_writer is
        given — only their counts (n_left_hits/n_right_hits), which is
        all any downstream consumer (the results table, panel design,
        exports) has ever actually used.
    """
    query_lengths = query_lengths or {}

    by_query: Dict[str, List[Dict]] = {}
    for h in hits:
        by_query.setdefault(h["qseqid"], []).append(h)

    ssr_ids = set()
    for q in by_query.keys():
        if "|" in q:
            ssr_ids.add(q.rsplit("|", 1)[0])

    results: Dict[str, Dict[str, Any]] = {}

    def _three_prime_mismatch_count(h, primer_len, window):
        """
        Count mismatches within the primer's actual 3' terminal `window`
        bases, based on the literal BLAST alignment (qseq/sseq) — not
        just whether they match exactly. Returns None if the alignment
        doesn't reach the primer's true 3' end at all (in which case the
        3' end isn't confirmed to match, so it's treated the same as
        before: not a confirmed competitor).
        """
        if h["qend"] != primer_len:
            return None  # alignment doesn't even reach the primer's true 3' end
        qseq, sseq = h.get("qseq", ""), h.get("sseq", "")
        if len(qseq) < window or len(sseq) < window:
            return None
        tail_q, tail_s = qseq[-window:].upper(), sseq[-window:].upper()
        return sum(1 for a, b in zip(tail_q, tail_s) if a != b or a == "-" or b == "-")

    def _good(h, sp):
        base_ok = (
            h["length"]   >= sp.min_align_len
            and h["mismatch"] <= sp.max_mismatches
            and h["pident"]   >= sp.min_identity
        )
        if not base_ok:
            return False
        if sp.require_3prime_match:
            primer_len = query_lengths.get(h["qseqid"])
            if primer_len is None:
                return True  # can't verify — fail open rather than reject everything
            mismatches = _three_prime_mismatch_count(h, primer_len, sp.three_prime_window)
            if mismatches is None:
                return False  # 3' end not confirmed to match — not a confirmed competitor
            # A hit is only safe to
            # ignore once it has AT LEAST three_prime_min_mismatches_to_ignore
            # mismatches in the terminal window. Fewer than that (including
            # a single terminal mismatch) still counts as a real competitor.
            return mismatches < sp.three_prime_min_mismatches_to_ignore
        return True

    def _unique_locations(hit_list):
        return {(h["sseqid"], h["sstrand"], h["sstart"]) for h in hit_list}

    for ssr_key in ssr_ids:
        # Parse (ssr_id, pair_rank) up front — needed both for the final
        # result key and, now, for hit_writer calls made right below.
        raw = ssr_key[3:] if ssr_key.startswith("SSR") else ssr_key
        if "_P" in raw:
            ssr_id_str, pair_rank_str = raw.rsplit("_P", 1)
        else:
            ssr_id_str, pair_rank_str = raw, "0"  # defensive fallback
        try:
            numeric_ssr_id = int(ssr_id_str)
        except ValueError:
            numeric_ssr_id = ssr_id_str
        try:
            pair_rank = int(pair_rank_str)
        except ValueError:
            pair_rank = pair_rank_str

        left_hits_all  = by_query.get(f"{ssr_key}|L", [])
        right_hits_all = by_query.get(f"{ssr_key}|R", [])

        left_hits  = [h for h in left_hits_all  if _good(h, specificity)]
        right_hits = [h for h in right_hits_all if _good(h, specificity)]

        if hit_writer is not None:
            for h in left_hits:
                hit_writer(numeric_ssr_id, pair_rank, "Forward", h)
            for h in right_hits:
                hit_writer(numeric_ssr_id, pair_rank, "Reverse", h)

        left_locations  = _unique_locations(left_hits)
        right_locations = _unique_locations(right_hits)
        n_left_hits     = len(left_locations)
        n_right_hits    = len(right_locations)

        MAX_HITS_PER_SIDE = 200
        left_search  = left_hits[:MAX_HITS_PER_SIDE]
        right_search = right_hits[:MAX_HITS_PER_SIDE]

        # Build EVERY possible amplicon genome-wide first — do NOT discard
        # anything by size yet. Classifying into "intended" (in-range) vs
        # "off-target" (out-of-range) has to happen AFTER collection, or
        # off-target amplicons can never be seen/evaluated at all (which
        # was the previous bug: allow_offtarget_amplicons/
        # max_offtarget_amplicon_size silently did nothing, since anything
        # they were meant to evaluate had already been discarded).
        all_amplicons = []
        for lh in left_search:
            for rh in right_search:
                if lh["sseqid"] != rh["sseqid"]:
                    continue
                if lh["sstrand"] == rh["sstrand"]:
                    continue
                coords = [lh["sstart"], lh["send"], rh["sstart"], rh["send"]]
                product = abs(max(coords) - min(coords)) + 1
                all_amplicons.append({
                    "contig": lh["sseqid"],
                    "product_size": product,
                    "left_hit": lh,
                    "right_hit": rh,
                })

        amplicons = []           # intended-size (in-range) amplicons
        offtarget_amplicons = [] # out-of-range amplicons
        for a in all_amplicons:
            if specificity.min_product <= a["product_size"] <= specificity.max_product:
                amplicons.append(a)
            else:
                offtarget_amplicons.append(a)

        n_amplicons = len(amplicons)

        # Gate 1: exactly one intended-size amplicon must form
        amplicon_ok = (n_amplicons == 1)

        # Gate 2: off-target products — either none are tolerated at all,
        # or any that exist must be LARGER than max_offtarget_amplicon_size
        # to be implausible competitors (a 3kb stray product essentially
        # never co-amplifies alongside a 150bp intended one under standard
        # PCR conditions). Products at or below that size are plausible
        # competitors and fail. Products shorter than min_product always
        # fail — short products amplify at least as efficiently as the
        # intended one — even if max_offtarget_amplicon_size is set lower.
        if specificity.allow_offtarget_amplicons:
            offtarget_ok = all(
                a["product_size"] > specificity.max_offtarget_amplicon_size
                and a["product_size"] >= specificity.min_product
                for a in offtarget_amplicons
            )
        else:
            offtarget_ok = (len(offtarget_amplicons) == 0)

        # Gate 3 (optional): each primer must map to at most max_total_hits
        # genomic locations, independent of whether an amplicon forms.
        if specificity.require_unique_hits:
            hits_ok = (
                n_left_hits  <= specificity.max_total_hits
                and n_right_hits <= specificity.max_total_hits
            )
        else:
            hits_ok = True

        status = "PASS" if (amplicon_ok and offtarget_ok and hits_ok) else "FAIL"

        # Amplicons/offtarget_amplicons keep their summary fields (contig,
        # product_size) but drop the nested full left_hit/right_hit
        # records — confirmed unused anywhere downstream, and each one
        # carries a complete alignment string pair that adds up for a
        # large primer set for no actual benefit.
        clean_amplicons = [
            {"contig": a["contig"], "product_size": a["product_size"]} for a in amplicons
        ]
        clean_offtarget_amplicons = [
            {"contig": a["contig"], "product_size": a["product_size"]} for a in offtarget_amplicons
        ]

        results[(numeric_ssr_id, pair_rank)] = {
            "ssr_id":              numeric_ssr_id,
            "pair_rank":           pair_rank,
            "specificity_status":  status,
            "n_left_hits":         n_left_hits,
            "n_right_hits":        n_right_hits,
            "n_amplicons":         n_amplicons,
            "amplicons":           clean_amplicons,
            "offtarget_amplicons": clean_offtarget_amplicons,
        }

    return results


# ---------------------------------------------------------------------------
# Performance: thread recommendation and BLAST database caching
# ---------------------------------------------------------------------------

def recommend_blast_threads(genome_size_bytes: int, cpu_count: int = None, available_memory_bytes=None) -> int:
    """
    Recommend a safe -num_threads value for blastn, based on available
    memory relative to genome size.

    Unlike SSR detection's own multiprocessing (which duplicates genome
    data into each separate worker process, especially costly on Windows
    with no copy-on-write process sharing), blastn's built-in threading
    shares ONE loaded database across all its threads — so the marginal
    memory cost of each additional thread is considerably lower than an
    additional multiprocessing worker would be. This function is
    correspondingly less conservative than recommend_worker_count() in
    core/ssr_detection.py, while still protecting against the kind of
    out-of-memory crash already observed on tight-memory hardware when
    BLAST is run against a large genome.

    This is a starting heuristic based on BLAST's general threading
    behaviour, not a precisely-benchmarked formula the way the SSR
    scanner's memory scaling is — if real-world experience suggests these
    thresholds should shift, they're a reasonable, isolated place to tune.
    """
    if cpu_count is None:
        cpu_count = os.cpu_count() or 1
    if available_memory_bytes is None:
        from core.ssr_detection import get_available_memory_bytes
        available_memory_bytes = get_available_memory_bytes()

    if available_memory_bytes is None:
        return min(2, cpu_count)

    ratio = available_memory_bytes / max(genome_size_bytes, 1)
    if ratio < 0.5:
        return 1
    elif ratio < 2.0:
        return min(2, cpu_count)
    else:
        return cpu_count


def _genome_content_key(genome: dict) -> str:
    """
    A fast, content-based identifier for a genome dict — based on sorted
    (contig name, sequence length, content sample) triples. Hashing every
    base of a multi-gigabase genome on every BLAST run would be its own
    significant cost, so this samples a bounded portion of each contig
    (its start, middle, and end) rather than the full sequence — cheap
    regardless of genome size, while still actually distinguishing two
    genomes that happen to share the same contig names and lengths but
    genuinely different sequence content (structural metadata alone,
    tested directly, is not sufficient: two such genomes would otherwise
    be silently treated as identical and incorrectly share a cached BLAST
    database, giving wrong results rather than just missing the cache).

    This is what lets the reconstructed FASTA (see
    _get_or_write_genome_fasta) and the BLAST database built from it
    (_get_or_build_blast_db) be correctly reused across repeat runs on
    the SAME loaded genome, rather than only matching by coincidence of
    an identical file path.
    """
    SAMPLE_LEN = 64
    parts = []
    for name, seq in sorted(genome.items()):
        n = len(seq)
        start = seq[:SAMPLE_LEN]
        end = seq[-SAMPLE_LEN:] if n > SAMPLE_LEN else ""
        mid = seq[n // 2: n // 2 + SAMPLE_LEN] if n > 2 * SAMPLE_LEN else ""
        parts.append(f"{name}:{n}:{start}:{mid}:{end}")
    raw = "|".join(parts)
    return hashlib.sha1(raw.encode("utf-8")).hexdigest()[:16]


def _get_or_write_genome_fasta(genome: dict) -> str:
    """
    Write `genome` to a FASTA file at a STABLE path derived from its
    content key, reusing an already-written file for the same genome
    instead of writing a fresh copy (to a fresh, differently-named temp
    path) on every single call. This matters specifically because the
    original file a user loaded might not be FASTA at all (GBFF, FASTQ
    are both supported load formats) — reconstructing a genuine FASTA is
    necessary regardless, but doing it to an unpredictable path each time
    silently defeated _get_or_build_blast_db()'s own path-based caching,
    since that cache can only recognise a repeat run if it's actually
    given the same path as before.
    """
    from core.fasta_loader import write_fasta
    cache_root = os.path.join(tempfile.gettempdir(), "grace_genome_fasta_cache")
    os.makedirs(cache_root, exist_ok=True)
    key = _genome_content_key(genome)
    path = os.path.join(cache_root, f"{key}.fasta")

    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path

    write_fasta(genome, path)
    return path


def _blast_db_cache_key(genome_fasta: str) -> str:
    """
    A stable identifier for a genome file's BLAST database cache entry,
    based on path + size + modification time — cheap to compute even for
    a multi-gigabase file (unlike a full checksum, which would need to
    read the entire file), while still correctly invalidating the cache
    if the genome file is replaced or edited.
    """
    st = os.stat(genome_fasta)
    raw = f"{os.path.abspath(genome_fasta)}|{st.st_size}|{st.st_mtime}"
    return hashlib.sha1(raw.encode("utf-8")).hexdigest()[:16]


def _get_or_build_blast_db(genome_fasta: str, makeblastdb: str) -> str:
    """
    Build the BLAST database for `genome_fasta` once and cache it in a
    persistent location keyed by _blast_db_cache_key(), reusing it on
    later calls against the same, unchanged genome file instead of
    rebuilding (and re-copying the whole genome) on every single
    specificity-screening run. Returns the db_prefix path to pass to
    blastn's -db argument.
    """
    cache_root = os.path.join(tempfile.gettempdir(), "grace_blast_db_cache")
    os.makedirs(cache_root, exist_ok=True)
    key = _blast_db_cache_key(genome_fasta)
    db_dir = os.path.join(cache_root, key)
    db_prefix = os.path.join(db_dir, "genome_db")

    marker = db_prefix + ".nin"
    marker_alt = db_prefix + ".nal"
    if os.path.exists(marker) or os.path.exists(marker_alt):
        return db_prefix

    os.makedirs(db_dir, exist_ok=True)
    _run_cmd([
        makeblastdb,
        "-in", genome_fasta,
        "-dbtype", "nucl",
        "-out", db_prefix,
    ])
    return db_prefix


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def check_specificity_blast(
    genome_fasta: str,
    primer_results: List[Dict[str, Any]],
    blast_params: Optional[BlastParams] = None,
    specificity_params: Optional[SpecificityParams] = None,
    blast_bin_dir: Optional[str] = None,
    progress_callback=None,
    batch_size: int = 25,
    raw_hits_path: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """
    Run BLAST-based specificity check on primer pairs.

    Args:
        genome_fasta: path to genome FASTA file (used to build BLAST db)
        primer_results: list of dicts with keys: ssr_id, left_primer, right_primer
        blast_params: BlastParams dataclass (defaults to BlastParams())
        specificity_params: SpecificityParams dataclass (defaults to SpecificityParams())
        blast_bin_dir: optional override path to BLAST binaries directory.
                       If None, automatically resolves bundled or system BLAST.
        progress_callback: optional callable(primers_done, primers_total),
            called after each batch completes. BLAST itself has no
            built-in progress mechanism for a single invocation covering
            every query — this is only possible at all because the
            database is now cached (see _get_or_build_blast_db()), which
            makes running blastn once per batch cheap: each call just
            loads the existing index rather than rebuilding it.
        batch_size: number of primer PAIRS (not individual query
            sequences — each pair's L/R queries always stay together
            within one batch) processed per blastn invocation.
        raw_hits_path: path to write every surviving BLAST hit's full
            detail (including alignment strings) to, as a tab-separated
            file, written incrementally as each batch completes rather
            than accumulated in memory for the whole run. If None, a
            fresh path in the system temp directory is generated
            automatically. Unlike the BLAST database/genome FASTA
            caches, this file is written fresh every call — it holds
            this run's actual result detail, not a reusable index.

    Returns:
        List of primer result dicts with added specificity keys. Hit-
        level detail is NOT included here at all (see raw_hits_path) —
        only per-locus summary fields (specificity_status, hit counts,
        amplicon sizes), which is all any downstream consumer has ever
        actually used. Holding every surviving hit's full record for
        every locus in memory for the whole run, only to serialise it
        all at the very end, was a real and avoidable cost for a large
        primer set against a repeat-dense genome.
    """
    if blast_params is None:
        blast_params = BlastParams()
    if specificity_params is None:
        specificity_params = SpecificityParams()

    if raw_hits_path is None:
        raw_hits_fd, raw_hits_path = tempfile.mkstemp(suffix=".tsv", prefix="grace_raw_blast_hits_")
        os.close(raw_hits_fd)

    exe_suffix = ".exe" if sys.platform == "win32" else ""

    resolved_bin_dir = blast_bin_dir or get_blast_bin_dir()

    if resolved_bin_dir:
        blastn      = os.path.join(resolved_bin_dir.strip(), f"blastn{exe_suffix}")
        makeblastdb = os.path.join(resolved_bin_dir.strip(), f"makeblastdb{exe_suffix}")
    else:
        blastn      = f"blastn{exe_suffix}"
        makeblastdb = f"makeblastdb{exe_suffix}"

    genome_fasta = genome_fasta.strip()

    # Built once per genome file and reused across runs — see
    # _get_or_build_blast_db()'s docstring. This also reads directly from
    # genome_fasta's own location rather than copying it first, which for
    # a multi-gigabase genome was previously a substantial, unnecessary
    # chunk of wall-clock time before BLAST had even started.
    db_prefix = _get_or_build_blast_db(genome_fasta, makeblastdb)

    n_threads = recommend_blast_threads(os.path.getsize(genome_fasta))

    query_lengths: Dict[str, int] = {}
    for p in primer_results:
        ssr_id    = p.get("ssr_id")
        pair_rank = p.get("pair_rank", 0)
        query_lengths[f"SSR{ssr_id}_P{pair_rank}|L"] = len(p["left_primer"])
        query_lengths[f"SSR{ssr_id}_P{pair_rank}|R"] = len(p["right_primer"])

    total = len(primer_results)
    paired: Dict[tuple, Dict[str, Any]] = {}

    RAW_HITS_HEADER = (
        "ssr_id\tpair_rank\tprimer_side\tqseqid\tsseqid\tpident\tlength\t"
        "mismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tsstrand\n"
    )
    raw_hits_fh = open(raw_hits_path, "w", encoding="utf-8")
    raw_hits_fh.write(RAW_HITS_HEADER)

    def _write_hit_row(ssr_id, pair_rank, side, h):
        raw_hits_fh.write(
            f"{ssr_id}\t{pair_rank}\t{side}\t{h['qseqid']}\t{h['sseqid']}\t{h['pident']}\t"
            f"{h['length']}\t{h['mismatch']}\t{h['gapopen']}\t{h['qstart']}\t{h['qend']}\t"
            f"{h['sstart']}\t{h['send']}\t{h['evalue']}\t{h['bitscore']}\t{h['sstrand']}\n"
        )

    try:
        for batch_start in range(0, total, batch_size):
            batch = primer_results[batch_start:batch_start + batch_size]

            with tempfile.TemporaryDirectory() as tmpdir:
                primers_fa = os.path.join(tmpdir, "primers.fasta")
                out_tab    = os.path.join(tmpdir, "blast_out.tsv")

                _write_primers_fasta(batch, primers_fa)

                _run_cmd([
                    blastn,
                    "-task", blast_params.task,
                    "-query", primers_fa,
                    "-db", db_prefix,
                    "-out", out_tab,
                    "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qseq sseq",
                    "-word_size", str(blast_params.word_size),
                    "-evalue", str(blast_params.evalue),
                    "-reward", str(blast_params.reward),
                    "-penalty", str(blast_params.penalty),
                    "-gapopen", str(blast_params.gapopen),
                    "-gapextend", str(blast_params.gapextend),
                    "-dust", blast_params.dust,
                    "-soft_masking", blast_params.soft_masking,
                    "-num_threads", str(n_threads),
                ])

                hits_raw = _parse_blast_tab(out_tab, specificity_params)

            # Same naming scheme as _write_primers_fasta — needed so the
            # 3'-end check can tell whether a hit's alignment actually
            # reaches each primer's true terminal base, not just how long
            # the alignment is.
            batch_paired = _pair_hits(hits_raw, specificity_params, query_lengths, hit_writer=_write_hit_row)
            paired.update(batch_paired)  # keys are (ssr_id, pair_rank) — unique per locus/rank, safe to merge

            if progress_callback is not None:
                progress_callback(min(batch_start + len(batch), total), total)
    finally:
        raw_hits_fh.close()

    out: List[Dict[str, Any]] = []
    for p in primer_results:
        ssr_id    = p.get("ssr_id")
        pair_rank = p.get("pair_rank", 0)
        spec = paired.get((ssr_id, pair_rank))
        if spec is None:
            merged = dict(p)
            merged["specificity_status"] = "UNKNOWN"
            merged["amplicons"] = []
            merged["offtarget_amplicons"] = []
            merged["n_left_hits"] = 0
            merged["n_right_hits"] = 0
        else:
            merged = dict(p)
            merged.update({
                "specificity_status": spec["specificity_status"],
                "amplicons": spec["amplicons"],
                "offtarget_amplicons": spec["offtarget_amplicons"],
                "n_left_hits": spec["n_left_hits"],
                "n_right_hits": spec["n_right_hits"],
            })
        out.append(merged)

    return out