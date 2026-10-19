"""
report_data.py
---------------
Pure-Python computation of every statistic and chart dataset used in the
PDF summary report. Deliberately has NO Qt/PyQt6 dependency at all — it
takes an AppState and returns plain dicts/lists of numbers, so it can be
unit-tested directly and reused identically by both the on-screen SSR
Summary charts and the PDF report generator, guaranteeing the two never
disagree with each other.

Each compute_*_summary() function returns None if that stage's data isn't
available in state — the report generator uses this to skip sections
cleanly, since GRACE is explicitly meant to produce a useful report at
whatever stage a user has actually reached, not only after a full run
through the entire pipeline.
"""

from collections import Counter


def compute_genome_summary(state) -> dict:
    if not state.genome:
        return None
    genome = state.genome
    lengths = {contig: len(seq) for contig, seq in genome.items()}
    total_bp = sum(lengths.values())
    n_contigs = len(lengths)
    largest = max(lengths.values()) if lengths else 0
    top_contigs = sorted(lengths.items(), key=lambda x: x[1], reverse=True)[:15]
    return {
        "filename": state.genome_filename,
        "n_contigs": n_contigs,
        "total_bp": total_bp,
        "largest_contig_bp": largest,
        "top_contigs_bar": top_contigs,  # list of (name, length_bp)
    }


def compute_ssr_summary(state) -> dict:
    if not state.ssrs:
        return None
    ssrs = state.ssrs
    genome = state.genome
    has_gff = state.gff_features is not None

    total_ssrs = len(ssrs)
    genome_size = sum(len(s) for s in genome.values()) if genome else 0
    density = (total_ssrs / (genome_size / 1_000_000)) if genome_size > 0 else 0

    motifs = [ssr.get("canonical_motif", ssr.get("motif", "")) for ssr in ssrs]
    motif_counts = Counter(motifs)
    unique_motifs = len(motif_counts)

    repeat_counts = [ssr.get("repeat_count", 0) for ssr in ssrs]
    mean_repeats = sum(repeat_counts) / total_ssrs if total_ssrs else 0

    motif_types = {
        "Mononucleotide": 0, "Dinucleotide": 0, "Trinucleotide": 0,
        "Tetranucleotide": 0, "Pentanucleotide": 0, "Hexanucleotide": 0,
        "Compound": 0, "Imperfect": 0,
    }
    for ssr in ssrs:
        if ssr.get("is_compound"):
            motif_types["Compound"] += 1
            continue
        if ssr.get("is_imperfect"):
            motif_types["Imperfect"] += 1
            continue
        length = len(ssr.get("canonical_motif", ssr.get("motif", "")))
        if length == 1:
            motif_types["Mononucleotide"] += 1
        elif length == 2:
            motif_types["Dinucleotide"] += 1
        elif length == 3:
            motif_types["Trinucleotide"] += 1
        elif length == 4:
            motif_types["Tetranucleotide"] += 1
        elif length == 5:
            motif_types["Pentanucleotide"] += 1
        elif length >= 6:
            motif_types["Hexanucleotide"] += 1
    motif_type_pie = {k: v for k, v in motif_types.items() if v > 0}

    feature_pie = None
    if has_gff:
        features = [ssr.get("genomic_feature", "intergenic") for ssr in ssrs]
        feature_counts = Counter(features)
        feature_pie = {
            "Exon": feature_counts.get("exon", 0),
            "Intron": feature_counts.get("intron", 0),
            "Intergenic": feature_counts.get("intergenic", 0),
        }
        feature_pie = {k: v for k, v in feature_pie.items() if v > 0}

    top_motifs_bar = motif_counts.most_common(15)

    counts = {}
    for ssr in ssrs:
        display = ssr.get("chrom_name", ssr.get("contig", "?"))
        counts[display] = counts.get(display, 0) + 1
    chrom_bar = sorted(counts.items(), key=lambda x: x[1], reverse=True)[:15]

    return {
        "total_ssrs": total_ssrs,
        "density_per_mb": density,
        "unique_motifs": unique_motifs,
        "mean_repeats": mean_repeats,
        "motif_type_pie": motif_type_pie,
        "feature_pie": feature_pie,
        "top_motifs_bar": top_motifs_bar,       # list of (motif, count)
        "chrom_bar": chrom_bar,                 # list of (name, count)
    }


def _bin_values(values, n_bins=10):
    """
    Simple equal-width histogram binning. Returns a list of
    (bin_label, count) tuples suitable for the shared bar-chart drawer.
    Returns [] for an empty input rather than raising.
    """
    if not values:
        return []
    lo, hi = min(values), max(values)
    if lo == hi:
        return [(f"{lo:g}", len(values))]
    width = (hi - lo) / n_bins
    bins = [0] * n_bins
    for v in values:
        idx = int((v - lo) / width)
        if idx >= n_bins:
            idx = n_bins - 1
        bins[idx] += 1

    # Pick just enough decimal places that the bin width survives rounding
    # — with whole-number rounding, a narrow-range dataset (e.g. Tm values
    # spanning only ~3°C across 10 bins) makes every bin edge round to the
    # same handful of values, so adjacent bins become visually identical.
    decimals = 0
    while round(width, decimals) == 0 and decimals < 4:
        decimals += 1
    fmt = f"{{:.{decimals}f}}"
    labels = [f"{fmt.format(lo + i*width)}-{fmt.format(lo + (i+1)*width)}" for i in range(n_bins)]
    return list(zip(labels, bins))


def compute_primer_summary(state) -> dict:
    if not state.primer_results:
        return None
    primers = state.primer_results
    n_designed = len(primers)
    n_ssrs_with_primers = len({p.get("ssr_id") for p in primers})

    product_sizes = [p["product_size"] for p in primers if p.get("product_size") is not None]
    tms = []
    for p in primers:
        lt, rt = p.get("left_tm"), p.get("right_tm")
        if lt is not None:
            tms.append(lt)
        if rt is not None:
            tms.append(rt)

    return {
        "n_designed": n_designed,
        "n_unique_ssrs": n_ssrs_with_primers,
        "product_size_hist": _bin_values(product_sizes, n_bins=10),
        "tm_hist": _bin_values(tms, n_bins=10),
    }


def compute_specificity_summary(state) -> dict:
    if not state.specificity_results:
        return None
    results = state.specificity_results
    n_pass = sum(1 for r in results if r.get("specificity_status") == "PASS")
    n_fail = sum(1 for r in results if r.get("specificity_status") == "FAIL")
    n_unknown = sum(1 for r in results if r.get("specificity_status") == "UNKNOWN")
    total = max(len(results), 1)
    pass_rate = 100 * n_pass / total
    n_unique_pass = len({r.get("ssr_id") for r in results if r.get("specificity_status") == "PASS"})

    status_pie = {"PASS": n_pass, "FAIL": n_fail, "UNKNOWN": n_unknown}
    status_pie = {k: v for k, v in status_pie.items() if v > 0}

    return {
        "n_pass": n_pass,
        "n_fail": n_fail,
        "n_unknown": n_unknown,
        "pass_rate": pass_rate,
        "n_unique_pass": n_unique_pass,
        "status_pie": status_pie,
    }


def compute_panel_summary(state) -> dict:
    result = state.panel_design_result
    if result is None:
        return None

    dye_bar = list(result.per_dye_counts.items())
    pool_counts = {}
    for m in result.markers:
        pool_counts[m.multiplex_pool] = pool_counts.get(m.multiplex_pool, 0) + 1
    pool_bar = sorted(
        [(f"Pool {p}", c) for p, c in pool_counts.items()],
        key=lambda x: int(x[0].split()[1])
    )

    markers_table = [
        {
            "pool": m.multiplex_pool,
            "dye": m.dye,
            "ssr_id": m.ssr_id,
            "motif": m.motif,
            "product_size": m.product_size,
            "forward_primer": m.record.get("left_primer", ""),
            "reverse_primer": m.record.get("right_primer", ""),
        }
        for m in sorted(result.markers, key=lambda m: (m.multiplex_pool, m.dye, m.product_size))
    ]

    return {
        "achieved": result.achieved_count,
        "requested": result.requested_count,
        "n_pools": len(result.pools),
        "dye_bar": dye_bar,     # list of (dye, count)
        "pool_bar": pool_bar,   # list of ("Pool N", count)
        "unplaced_reason": result.unplaced_reason,
        "markers_table": markers_table,
    }


def compute_full_report_data(state) -> dict:
    """
    Compute every section's data in one call. Sections with no available
    data are simply absent (None) — the PDF generator skips those,
    which is exactly what lets the report be generated usefully at any
    stage of the pipeline, not only after a full run through all of it.
    """
    return {
        "genome":      compute_genome_summary(state),
        "ssrs":        compute_ssr_summary(state),
        "primers":     compute_primer_summary(state),
        "specificity": compute_specificity_summary(state),
        "panel":       compute_panel_summary(state),
    }