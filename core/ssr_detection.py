"""
SSR Detection Module — High Performance
----------------------------------------
Detects SSRs (microsatellites) in a genome dictionary:
    genome_dict = {contig_name: sequence}

Performance strategy
--------------------
1. Multiprocessing  — contigs are distributed across all available CPU cores
                      using a process pool. Each worker scans its own batch
                      independently with no shared state.

2. Compiled regex   — one pre-compiled pattern per motif length using a
                      backreference: ([ACGTN]{k})(?:\\1){min_rep-1,}
                      re.finditer() runs in C, far faster than Python loops.

3. Pre-computed     — all possible canonical motifs for each k are computed
   canonical table    upfront (4^k possibilities) and stored in a plain dict.
                      No lru_cache overhead on the hot path.

4. Early skip       — contigs shorter than k * min_rep are skipped immediately.

5. Single uppercase — seq.upper() called once per contig, not per pattern.

6. Chunk-based      — for very large contigs the sequence is scanned in
   scanning           overlapping chunks to keep memory pressure manageable,
                      though for most genomes a full-contig scan is fine.
"""

import re
import os
import multiprocessing as mp
from functools import lru_cache
from itertools import product as iproduct

try:
    import numpy as _np
except ImportError:
    _np = None


# ---------------------------------------------------------
# CANDIDATE WINDOW FILTERING (fast coarse pass)
# ---------------------------------------------------------
#
# Regex backreference matching is correct but does a full backtracking
# scan across the ENTIRE sequence for every motif length. Most of a real
# genome has no repeat potential at all for a given k, so this finds
# candidate windows via a much cheaper numpy vectorized check first, and
# only runs the exact (unchanged, trusted) regex within those windows.
#
# IMPORTANT: this is a coarse FILTER only, never a replacement for the
# regex match itself — the final match boundaries always come from the
# same regex pattern used before, just applied to a smaller amount of
# text. An earlier prototype tried to determine match boundaries directly
# from the numpy periodicity check and was WRONG: a single coincidental
# base match at the boundary between two adjacent-but-different repeats
# (e.g. "ACACAC" immediately followed by "GCGCGC") can bridge them into
# one distorted, incorrectly-bounded region. Using numpy only to find
# candidate regions — with generous padding — and always deferring to the
# real regex for the actual match extraction avoids this entirely.

_WINDOW_PADDING = 10  # bp of slack on each side of a candidate region


def _find_candidate_windows(arr, k, threshold, padding=_WINDOW_PADDING):
    """
    Return a sorted, non-overlapping list of (start, end) windows where a
    period-k repeat of at least `threshold` bp could plausibly exist,
    found via a fast vectorized periodicity check (arr[i] == arr[i+k]).

    This is deliberately a superset of the true match regions (padded and
    then merged) — it exists only to skip full-sequence regex scanning
    over stretches with zero repeat potential, never to determine exact
    boundaries itself.
    """
    n = len(arr)
    if n <= k:
        return []
    match = (arr[:-k] == arr[k:])
    padded = _np.concatenate(([False], match, [False]))
    diff = _np.diff(padded.astype(_np.int8))
    run_starts = _np.where(diff == 1)[0]
    run_ends   = _np.where(diff == -1)[0]
    run_lengths = run_ends - run_starts
    qualifying = run_lengths >= threshold
    spans = []
    for s, e in zip(run_starts[qualifying], run_ends[qualifying]):
        spans.append((max(0, int(s) - padding), min(n, int(e) + k + padding)))
    if not spans:
        return []
    spans.sort()
    merged = [spans[0]]
    for s, e in spans[1:]:
        if s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def _finditer_windowed(pattern, seq, arr, k, threshold):
    """
    Yield regex matches for `pattern` (which matches unit-length k with
    the given repeat threshold), but only search within candidate windows
    identified by the fast numpy periodicity filter — not the whole
    sequence. Falls back to a full-sequence scan if numpy isn't available
    (e.g. not installed), preserving correctness either way.

    Match objects yielded have their positions already offset back to the
    original sequence's coordinates, and `.group(1)` still returns the
    matched unit exactly as the caller expects.
    """
    if _np is None:
        yield from pattern.finditer(seq)
        return

    windows = _find_candidate_windows(arr, k, threshold)
    for w_start, w_end in windows:
        window_seq = seq[w_start:w_end]
        for m in pattern.finditer(window_seq):
            yield _OffsetMatch(m, w_start)


class _OffsetMatch:
    """Thin wrapper so a match found within a window reports absolute
    (whole-sequence) positions via .start()/.end(), while .group() still
    delegates to the real match object unchanged."""
    __slots__ = ("_m", "_offset")

    def __init__(self, m, offset):
        self._m = m
        self._offset = offset

    def start(self, *args):
        return self._m.start(*args) + self._offset

    def end(self, *args):
        return self._m.end(*args) + self._offset

    def group(self, *args):
        return self._m.group(*args)


# ---------------------------------------------------------
# CANONICAL MOTIF
# ---------------------------------------------------------

def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def is_compound_motif(motif: str) -> bool:
    """
    Return True if motif is just a shorter unit repeated.
    e.g. CTCT -> True (it's CT*2), ATAT -> True (AT*2),
         CTAG -> False (no shorter period divides it evenly).
    Checks all divisors of len(motif) smaller than len(motif).
    """
    n = len(motif)
    for period in range(1, n):
        if n % period == 0:
            unit = motif[:period]
            if unit * (n // period) == motif:
                return True
    return False


@lru_cache(maxsize=131072)
def canonical_motif(motif: str, level: int = 4) -> str:
    """Compute the canonical form of a motif (cached)."""
    if level == 0:
        return motif
    motif = motif.upper()
    rotations = [motif[i:] + motif[:i] for i in range(len(motif))]
    best = min(rotations)
    if level >= 2:
        rc = reverse_complement(best)
        rc_rots = [rc[i:] + rc[:i] for i in range(len(rc))]
        best = min(best, min(rc_rots))
    return best


def _build_canonical_table(k: int, level: int) -> dict:
    """
    Pre-compute canonical forms for all possible k-mers over ACGT.
    Returns a dict {motif_str: canonical_str}.
    Much faster than calling canonical_motif() on every match at runtime.
    """
    bases = "ACGT"
    table = {}
    for combo in iproduct(bases, repeat=k):
        motif = "".join(combo)
        table[motif] = canonical_motif(motif, level)
    return table


# ---------------------------------------------------------
# PATTERN BUILDER
# ---------------------------------------------------------

def _build_pattern(k: int, min_rep: int) -> re.Pattern:
    """Compile a regex that matches any k-mer repeated >= min_rep times."""
    return re.compile(
        rf"([ACGT]{{{k}}})(?:\1){{{min_rep - 1},}}",
        re.IGNORECASE,
    )


# ---------------------------------------------------------
# SINGLE-CONTIG WORKER  (must be top-level for pickling)
# ---------------------------------------------------------

def _scan_contig(args):
    """
    Scan a single contig for SSRs.
    Args is a tuple for multiprocessing compatibility.
    Patterns and canonical tables are pre-built and passed in to avoid
    rebuilding them for every contig.
    Returns a list of SSR dicts (without ssr_id).
    """
    (
        contig, seq,
        motif_lengths, min_repeats,
        search_reverse,
        std_level,
        patterns, canonical_tables,   # pre-built, passed in
    ) = args

    seq = seq.upper()
    n   = len(seq)
    results = []

    seq_arr = _np.frombuffer(seq.encode("ascii", errors="replace"), dtype=_np.uint8) if _np is not None else None

    for k in motif_lengths:
        min_rep = min_repeats.get(k, 3)
        min_len = k * min_rep

        if n < min_len:
            continue

        pattern       = patterns[k]
        canonical_tbl = canonical_tables[k]
        threshold     = k * (min_rep - 1)

        # Forward strand
        for m in _finditer_windowed(pattern, seq, seq_arr, k, threshold):
            motif_str = m.group(1).upper()

            if k >= 2 and len(set(motif_str)) == 1:
                continue

            # Skip compound motifs — e.g. CTCT is just CT repeated, not a
            # true tetranucleotide. is_compound_motif catches all such cases.
            if is_compound_motif(motif_str):
                continue

            full_len     = m.end() - m.start()
            repeat_count = full_len // k
            canon        = canonical_tbl.get(motif_str) or canonical_motif(motif_str, std_level)

            results.append({
                "contig":          contig,
                "start":           m.start() + 1,  # 1-based,
                "end":             m.end(),         # 1-based inclusive
                "motif":           motif_str,
                "canonical_motif": canon,
                "repeat_count":    repeat_count,
            })

        # Reverse complement strand (optional)
        if search_reverse:
            rc_seq = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            rc_arr = _np.frombuffer(rc_seq.encode("ascii", errors="replace"), dtype=_np.uint8) if _np is not None else None
            for m in _finditer_windowed(pattern, rc_seq, rc_arr, k, threshold):
                motif_str = m.group(1).upper()

                if k >= 2 and len(set(motif_str)) == 1:
                    continue

                if is_compound_motif(motif_str):
                    continue

                full_len     = m.end() - m.start()
                repeat_count = full_len // k
                canon        = canonical_tbl.get(motif_str) or canonical_motif(motif_str, std_level)

                fwd_start = n - m.end()
                fwd_end   = n - m.start() - 1

                results.append({
                    "contig":          contig,
                    "start":           fwd_start + 1,  # 1-based
                    "end":             fwd_end + 1,     # 1-based inclusive
                    "motif":           motif_str,
                    "canonical_motif": canon,
                    "repeat_count":    repeat_count,
                })

    return results


# ---------------------------------------------------------
# MEMORY-AWARE WORKER COUNT
# ---------------------------------------------------------

def get_available_memory_bytes():
    """
    Return currently available system RAM in bytes, or None if it can't be
    determined (psutil not installed). Cross-platform via psutil.
    """
    try:
        import psutil
        return psutil.virtual_memory().available
    except Exception:
        return None


def recommend_worker_count(genome_size_bytes: int, cpu_count: int = None, available_memory_bytes=None) -> int:
    """
    Recommend a safe number of worker processes for SSR detection, scaling
    down automatically on memory-constrained machines without sacrificing
    parallelism on well-resourced ones.

    Rationale: `available_memory_bytes` should be queried AFTER the genome
    is already loaded into memory (as it always is by the time SSR
    detection runs) — so it already reflects the genome's own memory cost.
    What's left is the "surplus" room available for the extra pickled
    copies that multiprocessing IPC will create (unavoidable on Windows,
    which has no fork-based copy-on-write sharing). The safe worker count
    depends on how that surplus compares to the genome's own size:

        - Tight/moderate surplus (less than ~75% of the genome's own size)
          -> fully sequential (1 worker, zero duplication — the existing
             sequential code path already avoids multiprocessing entirely)
        - Generous surplus (75%-300% of genome size)
          -> a small, conservative worker count (2)
        - Very generous surplus (>300% of genome size)
          -> up to 4 workers — capped even here, rather than using every
             available core unconditionally, since a real crash occurred
             with a previous, less conservative version of this function
             (uncapped at the top band). More workers always means more
             duplicated memory regardless of how much surplus the
             available-memory estimate suggests; that estimate is itself
             imprecise, and being wrong in the risky direction costs a
             crash, not just a slower run.

    If available memory can't be determined (psutil not installed), falls
    back to a simpler size-only heuristic, itself tightened for the same
    reason: cap workers conservatively for genomes above 750Mb (previously
    1.5Gb), and even below that cap at 4 workers rather than the full CPU
    count, since this fallback path has zero visibility into actual
    memory pressure at all.

    Args:
        genome_size_bytes: total genome size (sum of all sequence lengths)
        cpu_count: number of CPU cores available (defaults to os.cpu_count())
        available_memory_bytes: override for testing; auto-detected otherwise

    Returns:
        Recommended n_workers value.
    """
    if cpu_count is None:
        cpu_count = max(1, os.cpu_count() or 1)

    if available_memory_bytes is None:
        available_memory_bytes = get_available_memory_bytes()

    if available_memory_bytes is None:
        # Can't detect available memory — fall back to a size-only heuristic,
        # tightened further since this path has zero visibility into
        # actual memory pressure.
        if genome_size_bytes > 750_000_000:
            return min(cpu_count, 2)
        return min(cpu_count, 4)

    if genome_size_bytes <= 0:
        return min(cpu_count, 4)

    # available_memory_bytes is queried AFTER the genome is already loaded,
    # so it already reflects the genome's own memory cost — don't subtract
    # it again here. Compare available "surplus" room directly against the
    # genome's size, since that's the scale of what a worker's chunk of
    # duplicated data could plausibly cost.
    ratio = available_memory_bytes / genome_size_bytes

    if ratio < 0.75:
        return 1
    elif ratio < 3.0:
        return min(cpu_count, 2)
    else:
        return min(cpu_count, 4)


# ---------------------------------------------------------
# PUBLIC API
# ---------------------------------------------------------

def find_ssrs(
    genome_dict: dict,
    motif_lengths: tuple = (2, 3, 4, 5, 6),
    min_repeats: dict = None,
    search_reverse: bool = False,
    motif_standardisation_level: int = 4,
    progress_callback=None,
    n_workers: int = None,
    min_contig_len: int = 300,
) -> list:
    """
    Detect SSRs in the genome using multiprocessing and compiled regex.

    Args:
        genome_dict: dict of {contig_name: sequence}
        motif_lengths: tuple of motif sizes to search, 1 (mononucleotide)
                       through 6 (hexanucleotide), e.g. (2, 3, 4)
        min_repeats: dict of {motif_length: min_repeat_count}; defaults to 3 for all
        search_reverse: also scan reverse complement strand
        motif_standardisation_level: 0-4, controls canonical motif computation
        progress_callback: optional callable(done, total) for UI progress updates
        n_workers: number of worker processes (default: all available CPU cores)
        min_contig_len: skip contigs shorter than this (default 300bp).
                        Contigs this short cannot yield usable SSR markers
                        and skipping them dramatically reduces overhead for
                        fragmented assemblies.

    Returns:
        List of dicts with fields:
            ssr_id, contig, start, end, motif, canonical_motif, repeat_count

    Note: a motif length of 1 (mononucleotide/homopolymer, e.g. AAAAAA) is
    searched only if explicitly included in motif_lengths. Runs of a single
    repeated base are always excluded from being counted as a motif length
    of 2 or higher (e.g. AAAA never counts as a fake dinucleotide "AA"
    repeat) — this is unconditional, not a separate toggle.
    """
    if min_repeats is None:
        min_repeats = {k: 3 for k in motif_lengths}

    if n_workers is None:
        n_workers = max(1, os.cpu_count() or 1)

    contigs       = list(genome_dict.items())
    total_contigs = len(contigs)

    # Filter out contigs that are too short to yield usable SSR markers
    # and build argument tuples for the remaining contigs
    skipped = sum(1 for _, seq in contigs if len(seq) < min_contig_len)
    contigs = [(c, s) for c, s in contigs if len(s) >= min_contig_len]
    total_contigs = len(contigs)

    if skipped:
        import sys
        print(
            f"[ssr_detection] Skipped {skipped:,} contigs shorter than "
            f"{min_contig_len}bp. Scanning {total_contigs:,} contigs.",
            file=sys.stderr,
        )

    if total_contigs == 0:
        return []

    # Pre-build patterns and canonical tables once — not once per contig
    patterns         = {k: _build_pattern(k, min_repeats.get(k, 3)) for k in motif_lengths}
    canonical_tables = {k: _build_canonical_table(k, motif_standardisation_level) for k in motif_lengths}

    task_args = [
        (
            contig, seq,
            motif_lengths, min_repeats,
            search_reverse,
            motif_standardisation_level,
            patterns, canonical_tables,
        )
        for contig, seq in contigs
    ]

    ssrs = []

    # Use multiprocessing for large genomes, single-process for small ones
    # (process pool overhead isn't worth it for < ~100 contigs)
    # Only fire the progress callback every 1% of contigs to avoid
    # hammering Streamlit with up to 250,000 re-renders on large genomes
    update_interval = max(1, total_contigs // 100)

    if n_workers > 1 and total_contigs >= 10:
        # Chunksize must scale DOWN to total_contigs/n_workers so work
        # actually spreads across every worker. The old max(50, ...) floor
        # broke this badly for chromosome-scale genomes (few, huge
        # contigs): with e.g. 30 contigs and 8 workers, it forced
        # chunksize=50 — larger than the whole contig list — meaning ALL
        # contigs collapsed into a single chunk sent to just ONE worker.
        # That means the entire genome gets pickled and IPC-transmitted to
        # one process instead of distributed across cores, multiplying
        # peak memory well beyond what's needed and giving zero
        # parallelism benefit for exactly the large-genome case where it
        # matters most. The upper cap (5000) still protects against
        # excessive per-task overhead for genomes with hundreds of
        # thousands of tiny scaffolds.
        chunksize = max(1, min(total_contigs // n_workers, 5000))

        # Use spawn context — prevents deadlock on re-run in PyInstaller frozen apps
        ctx = mp.get_context("spawn")
        with ctx.Pool(processes=n_workers) as pool:
            completed = 0
            for contig_results in pool.imap(
                _scan_contig, task_args, chunksize=chunksize
            ):
                ssrs.extend(contig_results)
                completed += 1
                if progress_callback and (
                    completed % update_interval == 0 or completed == total_contigs
                ):
                    progress_callback(completed, total_contigs)
    else:
        for i, args in enumerate(task_args):
            ssrs.extend(_scan_contig(args))
            if progress_callback and (
                (i + 1) % update_interval == 0 or i + 1 == total_contigs
            ):
                progress_callback(i + 1, total_contigs)

    # Sort by contig then position, then assign sequential IDs
    ssrs.sort(key=lambda x: (x["contig"], x["start"]))
    for i, s in enumerate(ssrs, start=1):
        s["ssr_id"] = i

    return ssrs


# ---------------------------------------------------------
# COMPOUND SSR MERGING (post-processing)
# ---------------------------------------------------------

def merge_compound_ssrs(ssrs: list, genome: dict, max_gap: int = 10) -> list:
    """
    Merge adjacent perfect SSRs of DIFFERENT canonical motifs into a single
    compound-locus record, e.g. (CA)8 immediately followed by (GA)6 becomes
    one compound SSR spanning both.

    This does NOT handle interrupted/imperfect repeats — i.e. a single
    motif's array broken up by a short mismatch (e.g. (CA)8-TT-(CA)5) is
    left untouched, since that requires fuzzy-matching detection rather
    than a merge of independently-detected perfect repeats. It IS however
    aware of components that were already extended by
    extend_imperfect_ssrs() — see below.

    Args:
        ssrs: list of perfect SSR dicts as returned by find_ssrs(), already
              sorted by (contig, start) — pass find_ssrs()'s own output
              directly, before any other filtering.
        genome: {contig: sequence} dict — needed to read the literal bases
                sitting in the gap between two components being merged.
        max_gap: maximum bp gap between one SSR's end and the next SSR's
                 start for them to be considered part of the same compound
                 locus. A gap of 0 means strictly back-to-back.

    Returns:
        A new list. Merged entries are flagged with "is_compound": True and
        gain a "components" key (the original SSR dicts that were combined,
        in genomic order). Their "motif" field becomes a fully-rendered
        display string showing each component as (UNIT)count (or its own
        imperfect structure string, if that component was itself extended
        by extend_imperfect_ssrs()), with the LITERAL lowercase gap
        sequence spliced in between whenever there was one — e.g.
        "(CA)8-ggg-(GA)6" for a 3bp gap, or "(CA)8" immediately followed by
        "(GA)6" with no separator at all for a 0bp (back-to-back) merge.
        "canonical_motif" stays a plain "/"-joined string (e.g. "AC/AG")
        for internal grouping purposes — it does not carry gap detail.
        "repeat_count" becomes the sum of component repeat counts
        (informational total, not a single meaningful count, since the
        components have different periods).

        Non-merged SSRs pass through unchanged, including their "ssr_id" —
        note that IDs are NOT reassigned by this function, since it may be
        called before or after other filters; call sites are responsible
        for reassigning "ssr_id" sequentially once all processing is done.
    """
    if not ssrs:
        return ssrs

    merged = []
    i = 0
    n = len(ssrs)
    while i < n:
        run = [ssrs[i]]
        j = i + 1
        while j < n:
            nxt = ssrs[j]
            if nxt["contig"] != run[0]["contig"]:
                break
            gap = nxt["start"] - run[-1]["end"] - 1
            if gap < 0 or gap > max_gap:
                # gap < 0: overlapping detections — don't merge these,
                # leave them as separate perfect entries rather than
                # guess at how to combine overlapping spans.
                break
            if nxt["canonical_motif"] == run[-1]["canonical_motif"]:
                # Same motif repeated nearby is an interrupted/imperfect
                # repeat scenario, not a compound one — stop the run here
                # rather than silently absorbing it.
                break
            run.append(nxt)
            j += 1

        if len(run) == 1:
            merged.append(run[0])
            i += 1
            continue

        first, last = run[0], run[-1]
        seq = genome[first["contig"]].upper()

        parts = []
        for idx, comp in enumerate(run):
            parts.append(comp.get("imperfect_structure") or f"({comp['motif']}){comp['repeat_count']}")
            if idx < len(run) - 1:
                gap_start0 = comp["end"]                    # 0-based, inclusive
                gap_end0   = run[idx + 1]["start"] - 1      # 0-based, exclusive
                gap_seq = seq[gap_start0:gap_end0]
                if gap_seq:
                    parts.append(gap_seq.lower())

        merged.append({
            "contig":          first["contig"],
            "start":           first["start"],
            "end":             last["end"],
            "motif":           "-".join(parts),
            "canonical_motif": "/".join(c["canonical_motif"] for c in run),
            "repeat_count":    sum(c["repeat_count"] for c in run),
            "is_compound":     True,
            "components":      run,
        })
        i = j

    return merged


# ---------------------------------------------------------
# IMPERFECT / INTERRUPTED SSR DETECTION (seed-and-extend)
# ---------------------------------------------------------

def _extend_forward(seq, unit, k, end0, max_interruption, min_extension_repeats):
    """
    Try to absorb ONE interruption immediately after `end0` (0-based,
    exclusive), provided `unit` resumes for at least `min_extension_repeats`
    full units right after a gap of at most `max_interruption` bp.

    Tries gap sizes smallest-first and returns on the first one that works
    (prefers the shortest interruption that still satisfies the resumption
    requirement).

    Returns (new_end0, absorbed_repeat_units, gap_seq_or_None) on success,
    else None. gap_seq is the literal interrupting bases (empty string
    when gap=0, i.e. no interruption was actually needed).

    Note: two different single-call str.find()-based restructurings of
    this loop were benchmarked and rejected — for the short patterns
    involved (2-6bp units), CPython's C-level short-string comparison is
    already close to optimal, and the restructured versions were measured
    to be marginally SLOWER in practice due to added string-building
    overhead, not faster. The genuine lever for this stage's performance
    at scale is parallelising extend_imperfect_ssrs() across candidates
    (see below), not micro-optimising this per-candidate function.
    """
    n = len(seq)
    for gap in range(0, max_interruption + 1):
        probe = end0 + gap
        count = 0
        pos = probe
        while pos + k <= n and seq[pos:pos + k] == unit:
            count += 1
            pos += k
        if count >= min_extension_repeats:
            gap_seq = seq[end0:probe] if gap > 0 else ""
            return pos, count, gap_seq
    return None


def _extend_backward(seq, unit, k, start0, max_interruption, min_extension_repeats):
    """Mirror of _extend_forward, growing `start0` (0-based, inclusive) backward."""
    for gap in range(0, max_interruption + 1):
        probe_end = start0 - gap  # exclusive
        count = 0
        pos = probe_end - k
        while pos >= 0 and seq[pos:pos + k] == unit:
            count += 1
            pos -= k
        if count >= min_extension_repeats:
            new_start0 = pos + k
            gap_seq = seq[probe_end:start0] if gap > 0 else ""
            return new_start0, count, gap_seq
    return None


def _render_structure(structure: list) -> str:
    """
    Render a locus's segment structure into a human-readable string like
    "(CA)8-tt-(CA)5". Repeat segments show as UPPERCASE (UNIT)count; gap
    segments show as the literal lowercase interrupting bases.
    """
    parts = []
    for seg in structure:
        if seg["type"] == "repeat":
            parts.append(f"({seg['unit']}){seg['count']}")
        else:
            parts.append(seg["seq"].lower())
    return "-".join(parts)


def _extend_contig_group(args):
    """
    Run the imperfect-extension pass for a single contig's perfect-SSR
    candidates. This is the unit of parallelisation: different contigs
    are fully independent, but candidates WITHIN one contig have a genuine
    sequential dependency (whether entry i gets consumed depends on
    whether an earlier entry already engulfed it), so per-contig is the
    correct and safe granularity — exactly mirroring how find_ssrs()
    itself parallelises across contigs.
    """
    contig, seq, group, max_interruption, min_extension_repeats = args
    seq = seq.upper()
    group_sorted = sorted(group, key=lambda x: x["start"])
    n_items = len(group_sorted)
    consumed = [False] * n_items
    computed = [None] * n_items

    for i, ssr in enumerate(group_sorted):
        if consumed[i]:
            continue
        k = len(ssr["motif"])
        unit = ssr["motif"].upper()
        start0 = ssr["start"] - 1
        end0   = ssr["end"]

        structure = [{"type": "repeat", "unit": unit, "count": ssr["repeat_count"]}]
        interruptions = []
        total_extra_units = 0

        while True:
            ext = _extend_forward(seq, unit, k, end0, max_interruption, min_extension_repeats)
            if ext is None:
                break
            new_end0, absorbed, gap_seq = ext
            if gap_seq:
                interruptions.append({"position": end0, "length": len(gap_seq), "sequence": gap_seq})
                structure.append({"type": "gap", "seq": gap_seq})
            structure.append({"type": "repeat", "unit": unit, "count": absorbed})
            end0 = new_end0
            total_extra_units += absorbed

        while True:
            ext = _extend_backward(seq, unit, k, start0, max_interruption, min_extension_repeats)
            if ext is None:
                break
            new_start0, absorbed, gap_seq = ext
            prepend = [{"type": "repeat", "unit": unit, "count": absorbed}]
            if gap_seq:
                interruptions.insert(0, {"position": new_start0 + absorbed * k, "length": len(gap_seq), "sequence": gap_seq})
                prepend.append({"type": "gap", "seq": gap_seq})
            structure = prepend + structure
            start0 = new_start0
            total_extra_units += absorbed

        if total_extra_units == 0:
            computed[i] = ssr
            continue

        new_ssr = dict(ssr)
        new_ssr["start"]               = start0 + 1
        new_ssr["end"]                 = end0
        new_ssr["repeat_count"]        = ssr["repeat_count"] + total_extra_units
        new_ssr["is_imperfect"]        = True
        new_ssr["interruptions"]       = interruptions
        new_ssr["imperfect_structure"] = _render_structure(structure)
        computed[i] = new_ssr

        for j in range(n_items):
            if j == i or consumed[j]:
                continue
            other = group_sorted[j]
            o_start0 = other["start"] - 1
            o_end0   = other["end"]
            if start0 <= o_start0 and o_end0 <= end0:
                consumed[j] = True

    return [computed[i] for i in range(n_items) if not consumed[i]]


def extend_imperfect_ssrs(
    ssrs: list,
    genome: dict,
    max_interruption: int = 3,
    min_extension_repeats: int = 2,
) -> list:
    """
    Extend perfect SSR seeds through short interruptions to detect
    imperfect/interrupted repeats, e.g. (CA)8-TT-(CA)5.

    Regex can only find perfect repeats, so this runs as a separate
    seed-and-extend pass: for each perfect SSR, tries growing its boundary
    outward in both directions, tolerating a short non-matching gap (at
    most `max_interruption` bp) provided the EXACT SAME repeat unit (the
    literal matched motif, not just an equivalent canonical form) resumes
    for at least `min_extension_repeats` more full units immediately after
    the gap. Multiple interruptions can be absorbed one at a time, as far
    outward as this rule keeps succeeding.

    Args:
        ssrs: perfect SSR list from find_ssrs(), BEFORE compound merging —
              extension changes start/end boundaries, and compound merging
              should operate on the final (possibly extended) boundaries.
        genome: {contig: sequence} dict — needed to inspect bases beyond
                each seed's current boundary.
        max_interruption: maximum bp gap tolerated per interruption.
        min_extension_repeats: minimum additional full repeat units that
                must immediately follow an interruption for it to be
                absorbed into the locus.

    Returns:
        A new list. Extended entries gain:
            is_imperfect:        True
            interruptions:       list of {"position", "length", "sequence"}
                                 dicts for each absorbed interruption, in
                                 genomic order (position is 0-based, the
                                 index where the gap starts)
            imperfect_structure: human-readable string of the whole locus,
                                 e.g. "(CA)8-tt-(CA)5" — repeat segments
                                 shown as (UNIT)count, interruptions shown
                                 as the literal lowercase interrupting bases
            repeat_count:        updated to the total perfect units across
                                 the whole extended array (original +
                                 absorbed)
        Non-extended SSRs pass through completely unchanged. Any other
        perfect SSR entries whose span ends up fully inside an extended
        locus are dropped (they're now part of the imperfect array, not
        a separate locus).

        Known limitation: this is a heuristic, not an alignment — in rare
        cases where a different repeat happens to sit immediately next to
        one being extended, and its sequence coincidentally matches the
        literal unit being extended after some gap, it could be absorbed
        as if it were an interruption.
    """
    if not ssrs:
        return ssrs

    from collections import defaultdict
    by_contig = defaultdict(list)
    for ssr in ssrs:
        by_contig[ssr["contig"]].append(ssr)

    result = []
    for contig, group in by_contig.items():
        task = (contig, genome[contig], group, max_interruption, min_extension_repeats)
        result.extend(_extend_contig_group(task))

    return result