"""
Primer Design Module (backend-safe)
-----------------------------------

Specificity is NOT integrated here by design.
Primers can be exported in BLAST-ready FASTA format.

Note on FASTA IDs for BLAST compatibility:
    Uses "|" as delimiter: SSR{id}|L and SSR{id}|R
    This avoids ambiguity when contig names or ssr_ids contain underscores.

Performance note:
    Multiprocessing is intentionally disabled. Spawning child processes from
    inside a PyQt6 QThread on Windows causes deadlocks due to Qt's event loop
    initialisation in the spawned process. primer3 is a C extension and runs
    fast enough in a single QThread for all practical SSR dataset sizes.

amplicon mode:
    Genotype-by-sequencing mode uses shorter flanks and smaller product sizes
    to ensure amplicons fit within Illumina read lengths.
"""

import re
import primer3

# Extra bp of template beyond the bare minimum needed to reach the
# requested max product size. Without this, Primer3 would have zero
# room to actually choose where to place primers (it would be forced
# to use almost the entire template edge-to-edge every time), making
# it far harder to also satisfy Tm/GC/complementarity constraints.
PRODUCT_SIZE_SAFETY_MARGIN = 50


# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------

# RECOMMENDED_PRESET - Used for Capillary Electrophoresis mode (default)
# Standard parameters for SSR genotyping with fragment analysis
RECOMMENDED_PRESET = {
    # Core primer constraints
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MAX_SIZE": 23,       # Primer3's own documented default
    "PRIMER_MIN_TM": 57.0,
    "PRIMER_OPT_TM": 58.0,
    "PRIMER_MAX_TM": 60.0,
    "PRIMER_MIN_GC": 30.0,
    "PRIMER_OPT_GC": 40.0,
    "PRIMER_MAX_GC": 60.0,

    # Secondary structure — hardcoded to Primer3's own documented
    # defaults, not exposed in the UI (advanced settings are no longer
    # user-configurable; see primer_panel.py)
    "PRIMER_MAX_POLY_X": 4,
    "PRIMER_GC_CLAMP": 0,
    "PRIMER_MAX_SELF_ANY": 8.0,
    "PRIMER_MAX_SELF_END": 3.0,
    "PRIMER_PAIR_MAX_COMPL_ANY": 8.0,
    "PRIMER_PAIR_MAX_COMPL_END": 3.0,
    "PRIMER_MAX_HAIRPIN_TH": 24.0,
    
    # Pair constraints
    "PRIMER_PAIR_MAX_DIFF_TM": 5.0,
    
    # Thermodynamic parameters — salt/dNTP hardcoded to Primer3's own
    # documented defaults; Tm formula and salt corrections deliberately
    # kept at SantaLucia (1998, value 1), Primer3's own RECOMMENDED
    # setting (see primer_panel.py tooltips) rather than its bare
    # backward-compatibility default of 0.
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_SALT_DIVALENT": 1.5,
    "PRIMER_DNTP_CONC": 0.6,
    "PRIMER_DNA_CONC": 50.0,
    "PRIMER_SALT_CORRECTIONS": 1,
    "PRIMER_TM_FORMULA": 1,
}

# amplicon_PRESET - Used for Amplicon Sequencing mode
# Optimized for Illumina amplicon sequencing
# Stricter Tm range and shorter products to fit 300bp paired-end reads
amplicon_PRESET = {
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MAX_SIZE": 23,       # Primer3's own documented default
    "PRIMER_MIN_TM": 58.0,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MAX_TM": 62.0,
    "PRIMER_MIN_GC": 35.0,
    "PRIMER_MAX_GC": 65.0,
    # Advanced settings hardcoded to the same Primer3 defaults as
    # capillary mode — the two modes now only differ in the basic,
    # user-adjustable settings above (Tm range, GC range, product size).
    "PRIMER_MAX_POLY_X": 4,
    "PRIMER_GC_CLAMP": 0,
    "PRIMER_MAX_SELF_ANY": 8.0,
    "PRIMER_MAX_SELF_END": 3.0,
    "PRIMER_PAIR_MAX_COMPL_ANY": 8.0,
    "PRIMER_PAIR_MAX_COMPL_END": 3.0,
    "PRIMER_MAX_HAIRPIN_TH": 24.0,
    "PRIMER_PAIR_MAX_DIFF_TM": 1.0,
    # Thermodynamic parameters
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_SALT_DIVALENT": 1.5,
    "PRIMER_DNTP_CONC": 0.6,
    "PRIMER_DNA_CONC": 50.0,
    "PRIMER_SALT_CORRECTIONS": 1,
    "PRIMER_TM_FORMULA": 1,
}

PRESETS = {
    "recommended": RECOMMENDED_PRESET,
    "amplicon":    amplicon_PRESET,
}

INT_KEYS = {
    "PRIMER_MIN_SIZE", "PRIMER_OPT_SIZE", "PRIMER_MAX_SIZE", "PRIMER_MAX_POLY_X",
    "PRIMER_GC_CLAMP", "PRIMER_SALT_CORRECTIONS", "PRIMER_TM_FORMULA",
}
FLOAT_KEYS = {
    "PRIMER_MIN_TM", "PRIMER_OPT_TM", "PRIMER_MAX_TM",
    "PRIMER_MIN_GC", "PRIMER_MAX_GC",
    "PRIMER_MAX_SELF_ANY", "PRIMER_MAX_SELF_END",
    "PRIMER_PAIR_MAX_COMPL_ANY", "PRIMER_PAIR_MAX_COMPL_END",
    "PRIMER_MAX_HAIRPIN_TH", "PRIMER_PAIR_MAX_DIFF_TM",
    "PRIMER_SALT_MONOVALENT", "PRIMER_SALT_DIVALENT", "PRIMER_DNTP_CONC",
    "PRIMER_DNA_CONC",
}



# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def normalize_types(opts):
    fixed = {}
    for key, value in opts.items():
        if key in INT_KEYS:
            fixed[key] = int(value)
        elif key in FLOAT_KEYS:
            fixed[key] = float(value)
        else:
            fixed[key] = value
    return fixed


def calc_gc(seq):
    seq = seq.upper()
    gc  = seq.count("G") + seq.count("C")
    return round((gc / len(seq)) * 100, 2)


_REPETITIVE_PRIMER_RE = {
    period: re.compile(rf"([ACGT]{{{period}}})\1{{3,}}")  # unit repeated 4+ times consecutively
    for period in range(2, 7)
}


def is_repetitive_primer(seq):
    """
    True if `seq` contains a run of a short repeated unit (period 2-6bp)
    appearing at least 4 times consecutively — e.g. a primer built
    substantially out of "CTCTCTCTCT" (a dinucleotide repeat).

    Primer3's own PRIMER_MAX_POLY_X setting only limits single-base
    homopolymer runs (AAAAA) — it has no concept of an ALTERNATING
    pattern like this at all, confirmed directly against Primer3's own
    documentation ("Max Homopolymer Length — maximum allowable length of
    a MONOnucleotide repeat"). A primer built mostly from a short
    repeated motif is exactly the kind of sequence expected to sit on,
    or very near, an SSR locus elsewhere in the genome — a real
    off-target risk that isn't reliably caught by BLAST screening either
    (BLAST's DUST low-complexity filter is deliberately disabled for
    this application, for good reason — see BlastParams — but that
    means BLAST's own statistics for a low-complexity query are less
    reliable than for a typical primer, not more). This check catches
    the problem at its source: such a primer should never be considered
    a candidate in the first place, independent of how BLAST would or
    wouldn't have scored it. period=1 (plain homopolymers) is
    deliberately excluded here since PRIMER_MAX_POLY_X already covers it.
    """
    seq = seq.upper()
    return any(pattern.search(seq) for pattern in _REPETITIVE_PRIMER_RE.values())


def calc_tm(seq):
    return round(primer3.calc_tm(seq), 2)


# Nearest-neighbor thermodynamic parameters (kcal/mol)
# These are NEGATIVE values - more negative = stronger binding
_NN_DG = {
    "AA": -1.0,  "AT": -0.88, "TA": -0.58, "CA": -1.45,
    "GT": -1.44, "CT": -1.28, "GA": -1.30, "CG": -2.17,
    "GC": -2.24, "GG": -1.84, "AC": -1.44, "AG": -1.28,
    "TC": -1.30, "TG": -1.45, "TT": -1.0,  "CC": -1.84,
}


def calc_3prime_dg(seq):
    """
    Calculate 3' end stability using nearest-neighbor thermodynamics.
    
    Evaluates the last 5 nucleotides (4 dinucleotide pairs) of the primer
    and sums their free energy (ΔG) values.
    
    Returns:
        ΔG in kcal/mol (NEGATIVE values)
        More negative = more stable binding
        Less negative = weaker binding
    
    Optimal range: -9 to -5 kcal/mol
        - Too stable (< -9 kcal/mol): Risk of non-specific priming
        - Too weak (> -5 kcal/mol): Poor primer extension by polymerase
    
    Note: Uses simplified nearest-neighbor parameters. For more accurate
    calculations accounting for salt concentration, use primer3.calc_end_stability()
    """
    if len(seq) < 5:
        return 0.0
    
    pentamer = seq[-5:].upper()
    # Sum the ΔG values (keep them NEGATIVE - do NOT use abs())
    total = sum(_NN_DG.get(pentamer[i:i+2], -1.0) for i in range(4))
    return round(total, 2)


# ---------------------------------------------------------------------------
# Low-complexity flank filter
# ---------------------------------------------------------------------------

# Pre-compiled patterns for fast repeat detection

# ---------------------------------------------------------------------------
# Core design function
# ---------------------------------------------------------------------------

def _design_one(ssr, template, left_bound, target_start, target_len,
                product_size_range, opts, num_pairs,
                amplicon_mode=False):
    """Design primers for a single SSR. Returns (success_list, failed_list)."""

    primer3_input = {
        "SEQUENCE_ID":               f"{ssr['contig']}|{ssr['ssr_id']}",
        "SEQUENCE_TEMPLATE":         template,
        "SEQUENCE_TARGET":           [target_start, target_len],
        "PRIMER_PRODUCT_SIZE_RANGE": [list(product_size_range)],
        # Keep stricter masking: exclude the SSR itself from primer binding
        # This prevents primers landing within the repeat 
        "SEQUENCE_EXCLUDED_REGION":  [[target_start, target_len]],
    }

    try:
        result = primer3.design_primers(primer3_input, opts)
    except Exception as e:
        # A genuinely degenerate template (e.g. an SSR sitting right at a
        # contig edge with almost no flanking sequence) could in principle
        # cause primer3's C bindings to raise rather than gracefully return
        # zero pairs. Catch it here so one bad locus can't crash the whole
        # batch — everyone else's SSRs still get designed, and this one
        # shows up as a clear failure rather than an unhandled crash.
        return [], [{**ssr, "fail_reason": "design_error", "design_error": str(e)}]

    success   = []
    found_any = False

    # Get salt concentrations for accurate 3' stability calculation
    mv_conc = opts.get("PRIMER_SALT_MONOVALENT", 50.0)
    dv_conc = opts.get("PRIMER_SALT_DIVALENT", 0.0)

    for rank in range(num_pairs):
        left_key = f"PRIMER_LEFT_{rank}_SEQUENCE"
        if left_key not in result:
            break
        found_any = True
        left_seq  = result[left_key]
        right_seq = result[f"PRIMER_RIGHT_{rank}_SEQUENCE"]

        # Reject primers that are themselves substantially a short
        # repeated unit (e.g. "CTCTCTCTCT") — see is_repetitive_primer()
        # for why Primer3's own settings can't catch this. Skipping here
        # (rather than filtering the whole `success` list afterward)
        # means this rank simply doesn't contribute a candidate, exactly
        # like any other Primer3-side rejection would.
        if is_repetitive_primer(left_seq) or is_repetitive_primer(right_seq):
            continue

        # Primer3's own reported (position, length) tuples, in TEMPLATE-
        # relative coordinates. Per Primer3's documented convention: for
        # a left primer, "position" is its LEFTMOST base; for a right
        # primer, "position" is its RIGHTMOST base (since it's on the
        # reverse strand) — these are NOT symmetric, and treating them
        # the same way would silently misplace every right primer's
        # computed genomic location.
        left_pos, left_len   = result[f"PRIMER_LEFT_{rank}"]
        right_pos, right_len = result[f"PRIMER_RIGHT_{rank}"]

        # Convert to absolute genomic coordinates, 0-based half-open
        # (i.e. genome[contig][start:end] gives exactly that primer's
        # sequence) — left_bound is the genomic position where this
        # locus's template itself begins.
        left_primer_genomic_start  = left_bound + left_pos
        left_primer_genomic_end    = left_primer_genomic_start + left_len
        right_primer_genomic_end   = left_bound + right_pos + 1
        right_primer_genomic_start = right_primer_genomic_end - right_len

        rec = {
            "ssr_id":          ssr["ssr_id"],
            "pair_rank":       rank,
            "contig":          ssr["contig"],
            "start":           ssr["start"],
            "end":             ssr["end"],
            "motif":           ssr.get("imperfect_structure") or ssr["motif"],
            "canonical_motif": ssr.get("canonical_motif"),
            "repeat_count":    ssr["repeat_count"],
            "left_primer":     left_seq,
            "right_primer":    right_seq,
            "product_size":    result[f"PRIMER_PAIR_{rank}_PRODUCT_SIZE"],
            "left_gc":         calc_gc(left_seq),
            "right_gc":        calc_gc(right_seq),
            "left_tm":         round(result[f"PRIMER_LEFT_{rank}_TM"], 2),
            "right_tm":        round(result[f"PRIMER_RIGHT_{rank}_TM"], 2),
            # Use standard thermodynamic ΔG calculation (negative values, -6 to -9 kcal/mol)
            "left_3end_dg":    calc_3prime_dg(left_seq),
            "right_3end_dg":   calc_3prime_dg(right_seq),
            "amplicon_mode":   amplicon_mode,
            # Genomic coordinates, 0-based half-open — needed to later
            # extract and highlight the actual amplicon (primers + SSR +
            # any sequence between them) directly from the genome for
            # the flanking-sequence export, rather than re-deriving these
            # positions from primer3's output a second time somewhere else.
            "left_primer_genomic_start":  left_primer_genomic_start,
            "left_primer_genomic_end":    left_primer_genomic_end,
            "right_primer_genomic_start": right_primer_genomic_start,
            "right_primer_genomic_end":   right_primer_genomic_end,
        }

        # Store genomic feature if available
        if "genomic_feature" in ssr:
            rec["genomic_feature"] = ssr["genomic_feature"]

        success.append(rec)

    if success:
        return success, []
    if found_any:
        # Primer3 did return pairs, but every one was rejected above as
        # repetitive — previously this locus was counted neither as
        # designed nor as failed, and silently vanished from the summary.
        return [], [{**ssr, "fail_reason": "repetitive_only"}]
    return [], [{**ssr, "fail_reason": "no_pairs",
                 "fail_detail": summarise_primer3_explain(result)}]


def _parse_explain(text):
    """
    Parse one primer3 *_EXPLAIN string, e.g.
    "considered 1083, GC content failed 1000, high tm 83, ok 0",
    into (considered, ok, {reason: count}).
    """
    considered, ok, reasons = 0, 0, {}
    for part in (text or "").split(","):
        part = part.strip()
        name, _, num = part.rpartition(" ")
        try:
            n = int(num)
        except ValueError:
            continue
        if name == "considered":
            considered = n
        elif name == "ok":
            ok = n
        elif n > 0:
            reasons[name.replace(" failed", "")] = n
    return considered, ok, reasons


def summarise_primer3_explain(result):
    """
    Turn primer3's PRIMER_{LEFT,RIGHT,PAIR}_EXPLAIN strings into a short,
    human-readable reason why no pair was found, naming the side that ran
    out of candidates and its most common rejection reasons.
    Returns "" if primer3 gave no explanation.
    """
    for key, label in (("PRIMER_LEFT_EXPLAIN",  "No acceptable forward primer"),
                       ("PRIMER_RIGHT_EXPLAIN", "No acceptable reverse primer"),
                       ("PRIMER_PAIR_EXPLAIN",  "No compatible primer pair")):
        considered, ok, reasons = _parse_explain(result.get(key))
        if ok == 0 and (considered > 0 or key == "PRIMER_PAIR_EXPLAIN"):
            if considered == 0:
                return f"{label} (no candidates in range)"
            top = sorted(reasons.items(), key=lambda kv: -kv[1])[:2]
            detail = ", ".join(f"{r} {n:,}" for r, n in top)
            return f"{label} — of {considered:,} considered: {detail}" if detail else label
    return ""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def design_primers_for_all_ssrs(
    genome,
    ssr_list,
    flank=100,
    product_size_range=(100, 350),
    preset="recommended",
    primer_opts=None,
    num_pairs=1,
    progress_callback=None,
    amplicon_mode=False,
):
    """
    Design primers for a list of SSRs using primer3.

    Runs entirely in the calling thread — no multiprocessing. Spawning child
    processes from inside a PyQt6 QThread on Windows causes deadlocks, and
    primer3 (a C extension) is fast enough in a single thread for all
    practical SSR dataset sizes.

    Args:
        genome: dict of {contig_name: sequence}
        ssr_list: list of SSR dicts
        flank: minimum bp of flanking sequence on each side of the SSR.
               This is a FLOOR, not a fixed value — if product_size_range's
               max requires more flanking sequence than this to be
               physically achievable, the flank actually used is grown
               automatically (see PRODUCT_SIZE_SAFETY_MARGIN). It will
               never be shrunk below what you pass in.
        product_size_range: (min, max) tuple for amplicon size
        preset: one of "strict", "recommended", "relaxed", "amplicon"
        primer_opts: dict of primer3 parameter overrides
        num_pairs: number of primer pairs to return per SSR (1–5)
        progress_callback: optional callable(done, total)
        amplicon_mode: if True, applies amplicon-specific design logic and adds
                  SEQUENCE_EXCLUDED_REGION to keep primers out of the repeat

    Returns:
        dict with keys:
            "success":      list of primer result dicts
            "failed":       list of SSR dicts for which no primers were found
                            (primer3 attempted every SSR). Each entry has
                            "fail_reason": "no_pairs" (primer3 found no valid
                            pair; primer3's explanation summarised under
                            "fail_detail"), "repetitive_only" (every pair
                            returned was rejected as repetitive), or
                            "design_error" (message under "design_error")
            "size_limited": list of SSR dicts where the contig was too short
                            (or the SSR too close to a contig edge) to reach
                            the requested max product size — primers were
                            still designed, just within whatever smaller
                            template was actually available
    """
    opts = PRESETS[preset].copy()

    if primer_opts:
        opts.update(primer_opts)

    opts = normalize_types(opts)

    num_pairs = max(1, min(int(num_pairs), 5))
    opts["PRIMER_NUM_RETURN"] = num_pairs

    total     = len(ssr_list)
    success   = []
    failed    = []
    size_limited = []  # designed, but contig too short to reach requested max product size

    # The bare minimum flank (per side) needed so that a product of
    # product_size_range[1] bp is physically achievable at all, given
    # this SSR's own length. Previously `flank` alone (a fixed 100bp by
    # default) silently determined the template size regardless of what
    # max product size the user actually asked for — e.g. with flank=100
    # the template is only ~2*100 + SSR_length bp, capping every product
    # at roughly 200-240bp even if product_max was set to 350 or higher.
    # `flank` is now treated as a floor: the effective flank used is
    # whichever is larger.
    requested_max = product_size_range[1]

    update_interval = max(1, total // 100)

    for i, ssr in enumerate(ssr_list):
        contig      = ssr["contig"]
        start       = ssr["start"]
        end         = ssr["end"]
        full_seq    = genome[contig]
        target_len  = end - (start - 1)

        min_flank_for_size = max(0, (requested_max - target_len + 1) // 2) + PRODUCT_SIZE_SAFETY_MARGIN
        effective_flank = max(flank, min_flank_for_size)

        left_bound   = max(0, (start - 1) - effective_flank)
        right_bound  = min(len(full_seq), end + effective_flank)
        template     = full_seq[left_bound:right_bound]
        target_start = (start - 1) - left_bound

        # Flag (but still attempt) SSRs where the contig is too short /
        # the SSR sits too close to a contig edge to actually reach the
        # requested max product size — Primer3 will just do its best
        # within whatever template is available.
        achievable_max = len(template)
        if achievable_max < requested_max:
            size_limited.append(ssr)

        # Every SSR gets a genuine primer3 attempt, regardless of flank
        # complexity or contig-edge proximity — primer3's own internal
        # thermodynamic scoring determines whether a valid pair exists,
        # and downstream BLAST specificity checking is what actually
        # catches primers landing in repetitive/non-unique regions. A
        # separate pre-filter here was rejecting loci before primer3 (or
        # BLAST) ever got a chance to make that call.
        s, f = _design_one(
            ssr, template, left_bound,
            target_start, target_len,
            product_size_range, opts, num_pairs,
            amplicon_mode=amplicon_mode,
        )
        success.extend(s)
        failed.extend(f)

        if progress_callback and (
            (i + 1) % update_interval == 0 or i + 1 == total
        ):
            progress_callback(i + 1, total)

    return {
        "success": success,
        "failed": failed,
        "size_limited": size_limited,
    }


def _extract_flanking_regions(genome: dict, pass_results: list, flank_bp: int = 10):
    """
    Shared extraction logic for both the HTML and FASTA flanking-sequence
    exports — the coordinate math (converting genomic primer/SSR
    positions into a tagged, contiguous-run breakdown of the exported
    region) lives here ONCE, so the two export formats can never quietly
    drift apart or develop different bugs in the same underlying logic.

    For every record in `pass_results`, extracts the genomic region
    spanning from `flank_bp` bases before the leftmost primer to
    `flank_bp` bases after the rightmost primer, and tags each base as
    "primer", "ssr", or plain flanking sequence.

    Returns (regions, skipped_count), where each region is a dict:
        {rec, contig, export_start, export_end, seq, spans}
    spans is a list of (tag_or_None, chunk_string) tuples in order,
    covering the whole exported sequence.

    Records missing the genomic-coordinate fields (added specifically
    for this export; older session data saved before this feature
    existed won't have them) are skipped and counted rather than
    causing the whole export to fail.
    """
    regions = []
    skipped = 0

    for rec in pass_results:
        required = ("contig", "left_primer_genomic_start", "left_primer_genomic_end",
                    "right_primer_genomic_start", "right_primer_genomic_end", "start", "end")
        if any(k not in rec or rec[k] is None for k in required):
            skipped += 1
            continue

        contig = rec["contig"]
        if contig not in genome:
            skipped += 1
            continue
        contig_seq = genome[contig]

        amplicon_start = rec["left_primer_genomic_start"]
        amplicon_end   = rec["right_primer_genomic_end"]  # 0-based, exclusive

        export_start = max(0, amplicon_start - flank_bp)
        export_end   = min(len(contig_seq), amplicon_end + flank_bp)

        # Per-base tag, built once and sliced into contiguous runs below
        # — robust even if regions were ever adjacent or (in principle)
        # overlapping, since primer tags are applied first and the SSR
        # tag applied after simply takes priority for its own span,
        # without needing special-case merge logic.
        tags = [None] * (export_end - export_start)

        def _tag(g_start, g_end, label):
            s = max(g_start, export_start) - export_start
            e = min(g_end, export_end) - export_start
            for i in range(max(0, s), max(0, e)):
                tags[i] = label

        _tag(rec["left_primer_genomic_start"], rec["left_primer_genomic_end"], "primer")
        _tag(rec["right_primer_genomic_start"], rec["right_primer_genomic_end"], "primer")
        # SSR coordinates (rec["start"]/["end"]) are 1-based inclusive,
        # matching the convention used throughout the rest of GRACE —
        # convert to the same 0-based-exclusive convention used above.
        _tag(rec["start"] - 1, rec["end"], "ssr")

        seq = str(contig_seq[export_start:export_end]).upper()

        spans = []
        run_start = 0
        for i in range(1, len(tags) + 1):
            if i == len(tags) or tags[i] != tags[run_start]:
                spans.append((tags[run_start], seq[run_start:i]))
                run_start = i
        if not spans and seq:
            spans = [(None, seq)]

        regions.append({
            "rec": rec, "contig": contig,
            "export_start": export_start, "export_end": export_end,
            "seq": seq, "spans": spans,
        })

    return regions, skipped


def build_flanking_sequence_html(genome: dict, pass_results: list, flank_bp: int = 10) -> str:
    """
    Build a self-contained HTML document showing, for every PASS-status
    primer record in `pass_results`, the genomic region spanning from
    `flank_bp` bases before the leftmost primer to `flank_bp` bases
    after the rightmost primer — with the forward primer, reverse
    primer, and the SSR itself highlighted in distinct colours for easy
    visual inspection, and plain flanking sequence shown around them.

    Returns the complete HTML document as a string.
    """
    PRIMER_COLOR = "#a8d5ff"
    SSR_COLOR    = "#ffd699"

    def esc(s):
        return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    regions, skipped = _extract_flanking_regions(genome, pass_results, flank_bp)

    blocks = []
    for region in regions:
        rec = region["rec"]
        html_seq = ""
        for tag, chunk in region["spans"]:
            if tag == "primer":
                html_seq += f'<span class="primer">{esc(chunk)}</span>'
            elif tag == "ssr":
                html_seq += f'<span class="ssr">{esc(chunk)}</span>'
            else:
                html_seq += esc(chunk)

        pair_rank = rec.get("pair_rank", 0)
        blocks.append(f"""
        <div class="locus">
            <div class="locus-header">
                SSR {esc(rec['ssr_id'])} &middot; pair {esc(pair_rank)} &middot;
                motif {esc(rec.get('canonical_motif') or rec.get('motif', ''))} &middot;
                product {esc(rec.get('product_size', '?'))}bp &middot;
                {esc(region['contig'])}:{region['export_start'] + 1}-{region['export_end']}
            </div>
            <div class="seq">{html_seq}</div>
        </div>""")

    n_included = len(blocks)
    skipped_note = (
        f'<p class="note">{skipped} record(s) were skipped — missing genomic '
        f'coordinate data (likely from a session saved before this export existed).</p>'
        if skipped else ""
    )

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>GRACE — Flanking Sequences</title>
<style>
    body {{ font-family: -apple-system, Segoe UI, Arial, sans-serif; background: #fafbfc; color: #1f2328; margin: 24px; }}
    h1 {{ font-size: 20px; }}
    .legend span {{ display: inline-block; padding: 2px 8px; border-radius: 3px; margin-right: 12px; font-size: 13px; }}
    .legend .primer {{ background: {PRIMER_COLOR}; }}
    .legend .ssr {{ background: {SSR_COLOR}; }}
    .note {{ color: #b45309; font-size: 13px; }}
    .locus {{ margin: 18px 0; padding: 12px; background: #ffffff; border: 1px solid #d0d7de; border-radius: 6px; }}
    .locus-header {{ font-size: 12px; color: #59636e; margin-bottom: 8px; }}
    .seq {{ font-family: Consolas, "Courier New", monospace; font-size: 14px; line-height: 1.6; word-break: break-all; }}
    .seq .primer {{ background: {PRIMER_COLOR}; }}
    .seq .ssr {{ background: {SSR_COLOR}; }}
</style>
</head>
<body>
<h1>GRACE &mdash; Flanking Sequences</h1>
<p class="legend">
    <span class="primer">Primer (forward &amp; reverse)</span>
    <span class="ssr">SSR repeat region</span>
    Plain text: flanking sequence (&plusmn;{flank_bp}bp beyond the primers)
</p>
{skipped_note}
<p>{n_included} locus/loci included.</p>
{''.join(blocks)}
</body>
</html>
"""


def build_flanking_sequence_fasta(genome: dict, pass_results: list, flank_bp: int = 10) -> str:
    """
    Build a plain-text, FASTA-formatted export of the same flanking
    regions as build_flanking_sequence_html() — for users who want the
    raw sequence for further use (pasting into another tool, manual
    BLAST, etc.) rather than a visual inspection view.

    Since plain FASTA has no way to show colour, the position of each
    primer and the SSR within the sequence is instead encoded directly
    in the header line as 1-based, inclusive positions relative to the
    start of that entry's own sequence — e.g. a header noting
    "left_primer=1-20 ssr=35-50 right_primer=65-84" tells you exactly
    where each feature sits without needing to count bases by eye.

    Returns the complete FASTA text as a string.
    """
    regions, skipped = _extract_flanking_regions(genome, pass_results, flank_bp)

    lines = []
    for region in regions:
        rec = region["rec"]
        # Convert each span into a 1-based, inclusive, sequence-relative
        # position range for the header.
        pos = 1
        feature_ranges = {"primer": [], "ssr": []}
        for tag, chunk in region["spans"]:
            if tag in feature_ranges:
                feature_ranges[tag].append((pos, pos + len(chunk) - 1))
            pos += len(chunk)

        primer_desc = " ".join(f"primer={s}-{e}" for s, e in feature_ranges["primer"])
        ssr_desc = " ".join(f"ssr={s}-{e}" for s, e in feature_ranges["ssr"])

        header = (
            f">SSR{rec['ssr_id']}_pair{rec.get('pair_rank', 0)} "
            f"motif={rec.get('canonical_motif') or rec.get('motif', '')} "
            f"product={rec.get('product_size', '?')}bp "
            f"location={region['contig']}:{region['export_start'] + 1}-{region['export_end']} "
            f"{primer_desc} {ssr_desc}"
        )
        lines.append(header)
        lines.append(region["seq"])

    if skipped:
        lines.append(
            f"; {skipped} record(s) were skipped — missing genomic coordinate data "
            f"(likely from a session saved before this export existed)."
        )

    return "\n".join(lines) + "\n" if lines else ""



def primers_to_blast_fasta(primer_results):
    """
    Convert primer results into a BLAST-ready FASTA string.
    Uses bare primer sequences for BLAST.
    """
    lines = []
    for rec in primer_results:
        ssr_id = rec["ssr_id"]
        lines.append(f">SSR{ssr_id}|L")
        lines.append(rec["left_primer"])
        lines.append(f">SSR{ssr_id}|R")
        lines.append(rec["right_primer"])
    return "\n".join(lines) + "\n" if lines else ""