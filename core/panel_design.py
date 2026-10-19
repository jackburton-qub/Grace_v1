"""
panel_design.py
----------------
Builds an optimised genotyping panel from PASS-specificity primer pairs:
groups loci into Ta-compatible multiplex PCR pools, assigns each pool a
balanced dye spread, and enforces same-dye product-size separation WITHIN
each multiplex pool to avoid fragment overlap during capillary electrophoresis.

This is a greedy, constructive optimiser — not an exhaustive/proven-optimal
solver. For realistic candidate counts (tens to low hundreds of PASS loci)
an exhaustive search over dye/rank/grouping assignments is computationally
infeasible and unnecessary; a well-designed greedy approach gives a strong,
wet-lab-practical result.

Key design decisions:

- Ta-compatible clustering happens FIRST, at the level of individual
  candidate (locus, rank) pairs — not after selection. A single locus can
  have multiple ranks with genuinely different Ta values, so it isn't
  tied to one cluster; whichever rank actually gets used depends on which
  multiplex it ends up filling. Clusters are processed largest-available-
  locus-count first, so the panel favours fewer, larger multiplexes over
  needlessly fragmenting into many small ones — reaching a requested
  count of e.g. 20 markers might mean one multiplex of 20, or two of
  10+10, or 9+11, entirely driven by what the actual Ta distribution
  among PASS candidates supports.

- Each multiplex pool gets its own internally-balanced dye spread (not
  just the panel as a whole) — a real multiplex PCR reaction should
  itself use its available dyes sensibly, not end up e.g. all-FAM by
  chance while another reaction has none.

- Product-size separation on a shared dye is enforced WITHIN each
  multiplex pool independently, not across the whole panel — a marker in
  one multiplex reaction does not constrain product sizes used in a
  different multiplex, even on the same dye. Pools are built candidate-
  by-candidate using first-fit: each locus tries every already-open pool
  for its Ta-cluster before a new one is opened, which is what lets
  max_per_multiplex-capped pools fill up evenly in practice without a
  separate rebalancing step.

- Motif-type preference (e.g. "prioritise tri/tetranucleotides") is a
  SOFT priority, not a hard filter: preferred-motif loci are always
  considered first when filling a multiplex, but if there aren't enough
  of them to reach the requested count, other motif types fill the
  remainder rather than leaving the panel short. Compound loci (joined
  canonical_motif strings like "AC/AG") don't represent one specific
  motif length, so they're treated as their own category and only
  prioritised when no specific preference is set.

- Primer-primer interactions are checked as each locus is placed: its
  two primers are scored against every primer already in the candidate
  pool with primer3's heterodimer and 3'-end stability calculations. A
  pool where any 3'-anchored interaction or any dimer is more negative
  than its threshold rejects the placement, and the locus moves on to
  the next pool (or a fresh one) exactly as it would for a size clash.
  Interactions that come within NEAR_THRESHOLD_MARGIN of a threshold in
  the final pools are reported, as are the rejected placements.
"""

from dataclasses import dataclass, field
from typing import Optional

from primer3.thermoanalysis import ThermoAnalysis


DYE_ORDER = ["FAM", "HEX", "NED", "PET"]

DEFAULT_END_DG_THRESHOLD = -6.0     # kcal/mol, 3'-anchored interactions
DEFAULT_DIMER_DG_THRESHOLD = -9.0   # kcal/mol, any heterodimer
NEAR_THRESHOLD_MARGIN = 1.0         # kcal/mol above a threshold still worth reporting

# Same salt/dNTP/DNA concentrations primer design uses (see the presets in
# core/primer_design.py), so interaction energies are on the same footing
# as the per-primer checks Primer3 already applied.
_THERMO_ARGS = {"mv_conc": 50.0, "dv_conc": 1.5, "dntp_conc": 0.6, "dna_conc": 50.0}

INTERACTION_END = "3' end"
INTERACTION_DIMER = "dimer"


@dataclass
class PanelCandidate:
    """One PASS-specificity primer pair, ready for panel assignment."""
    ssr_id: object
    pair_rank: int
    contig: str
    motif: str
    canonical_motif: str
    product_size: int
    ta: float               # representative annealing temperature for this pair
    record: dict             # the original primer/specificity row, kept for output


@dataclass
class PanelMarker:
    """A locus assigned a place in the final panel."""
    ssr_id: object
    pair_rank: int
    contig: str
    motif: str
    product_size: int
    ta: float
    dye: str
    multiplex_pool: int
    record: dict


@dataclass
class PrimerInteraction:
    """
    One primer-primer interaction worth reporting.

    status is "REJECTED" when it stopped locus A being placed in that
    pool (A was then placed elsewhere), "FLAGGED" when it is beyond its
    threshold inside a final pool, or "NEAR" when it is within
    NEAR_THRESHOLD_MARGIN of the threshold inside a final pool.
    """
    status: str
    kind: str               # INTERACTION_END or INTERACTION_DIMER
    dg: float               # kcal/mol
    threshold: float        # kcal/mol
    multiplex_pool: int
    ssr_id_a: object
    pair_rank_a: int
    primer_a: str           # "F" or "R"
    seq_a: str
    ssr_id_b: object
    pair_rank_b: int
    primer_b: str
    seq_b: str


@dataclass
class PanelResult:
    markers: list                  # list of PanelMarker, the final panel
    pools: list                    # list of list[PanelMarker], grouped by multiplex
    requested_count: int
    achieved_count: int
    dyes_used: list
    per_dye_counts: dict
    unplaced_reason: Optional[str] = None  # explanation if achieved < requested
    interactions: list = field(default_factory=list)  # list of PrimerInteraction
    end_dg_threshold: Optional[float] = None
    dimer_dg_threshold: Optional[float] = None


def _sort_key(ssr_id):
    """Stable sort key that works whether ssr_id is int, str, or mixed
    across a dataset."""
    return (str(type(ssr_id)), ssr_id)


def _spread_order(items: list) -> list:
    """
    Reorder items (already in genome order) so that any prefix of the
    result is spread evenly across the original list: first item, then
    the midpoint, then the quarter points, then the eighths, and so on
    (a base-2 van der Corput sequence over the indices). The greedy
    placement takes loci in this order, so stopping at target_count
    still samples the whole genome rather than just its first contigs.
    """
    n = len(items)
    if n <= 2:
        return list(items)
    order, seen = [], set()
    k = 0
    while len(order) < n:
        frac, denom, x = 0.0, 1.0, k
        while x:
            denom *= 2
            frac += (x & 1) / denom
            x >>= 1
        idx = int(frac * n)
        if idx not in seen:
            seen.add(idx)
            order.append(items[idx])
        k += 1
    return order


def _build_candidates(specificity_results: list) -> dict:
    """
    Group PASS rows from specificity_results into {ssr_id: [PanelCandidate, ...]}.
    Ta is approximated as the mean of left_tm/right_tm — primer3 reports
    each primer's own Tm, not a single pair "annealing temperature", so
    this mean is used as a practical single representative value per pair
    for multiplex-compatibility clustering.
    """
    by_locus = {}
    for row in specificity_results:
        if row.get("specificity_status") != "PASS":
            continue
        ssr_id = row.get("ssr_id")
        left_tm = row.get("left_tm")
        right_tm = row.get("right_tm")
        if left_tm is None or right_tm is None:
            continue
        ta = (left_tm + right_tm) / 2.0
        cand = PanelCandidate(
            ssr_id=ssr_id,
            pair_rank=row.get("pair_rank", 0),
            contig=row.get("contig"),
            motif=row.get("motif", ""),
            canonical_motif=row.get("canonical_motif") or "",
            product_size=row.get("product_size"),
            ta=ta,
            record=row,
        )
        by_locus.setdefault(ssr_id, []).append(cand)
    return by_locus


def _fits_dye(product_size: int, occupied_sizes: list, min_gap: int) -> bool:
    """True if product_size is at least min_gap bp from every size already
    assigned to this dye."""
    return all(abs(product_size - s) >= min_gap for s in occupied_sizes)


def _cluster_by_ta(flat_candidates: list, ta_tolerance: float) -> list:
    """
    Sort candidates by Ta and greedily group into clusters where every
    member is within ta_tolerance of the cluster's minimum Ta. A simple,
    effective 1-D interval clustering.
    """
    ordered = sorted(flat_candidates, key=lambda c: c.ta)
    clusters = []
    current = [ordered[0]]
    cluster_min_ta = ordered[0].ta
    for c in ordered[1:]:
        if c.ta - cluster_min_ta <= ta_tolerance:
            current.append(c)
        else:
            clusters.append(current)
            current = [c]
            cluster_min_ta = c.ta
    clusters.append(current)
    return clusters


def primer_interaction_energies(seq_a: str, seq_b: str, thermo: Optional[ThermoAnalysis] = None) -> tuple:
    """
    Return (dimer_dg, end_dg) in kcal/mol for two primers. dimer_dg is
    primer3's heterodimer ΔG; end_dg is the more negative of the two
    3'-end stabilities (seq_a's 3' end on seq_b, and seq_b's on seq_a).
    Pass a ThermoAnalysis to reuse across many calls - it gives the same
    numbers as primer3.calc_heterodimer / calc_end_stability, faster.
    """
    thermo = thermo or ThermoAnalysis(**_THERMO_ARGS)
    dimer = thermo.calc_heterodimer(seq_a, seq_b).dg / 1000.0
    end = min(
        thermo.calc_end_stability(seq_a, seq_b).dg,
        thermo.calc_end_stability(seq_b, seq_a).dg,
    ) / 1000.0
    return dimer, end


def _pair_primers(record: dict) -> list:
    """[(label, sequence)] for the forward/reverse primers that are present."""
    out = []
    for label, key in (("F", "left_primer"), ("R", "right_primer")):
        seq = (record.get(key) or "").strip().upper()
        if seq:
            out.append((label, seq))
    return out


def _motif_matches_preference(canonical_motif: str, preferred_motif_lengths) -> bool:
    """
    True if this locus's motif length is among the user's preferred
    lengths. Compound loci (joined "AC/AG"-style canonical_motif) don't
    represent one specific length and never match an explicit preference.
    """
    if not preferred_motif_lengths:
        return True  # no preference set — everything matches equally
    if "/" in canonical_motif:
        return False
    return len(canonical_motif) in preferred_motif_lengths


def design_panel(
    specificity_results: list,
    dyes: list,
    target_count: int,
    min_dye_gap: int = 75,
    ta_tolerance: float = 1.0,
    enforce_contig_spread: bool = True,
    preferred_motif_lengths: Optional[list] = None,
    max_per_multiplex: Optional[int] = None,
    end_dg_threshold: float = DEFAULT_END_DG_THRESHOLD,
    dimer_dg_threshold: float = DEFAULT_DIMER_DG_THRESHOLD,
) -> PanelResult:
    """
    Build an optimised multiplex genotyping panel.

    Args:
        specificity_results: full list of primer/specificity rows (as
            produced by check_specificity_blast) — only PASS rows are used.
            May contain multiple rows (ranks) per ssr_id.
        dyes: list of dye names to use, e.g. ["FAM", "HEX", "NED", "PET"].
        target_count: total number of markers wanted in the final panel,
            across ALL multiplex pools combined — the system decides how
            many pools are needed and how the count splits between them,
            driven by the actual Ta distribution among PASS candidates.
        min_dye_gap: minimum bp separation required between two markers
            assigned the SAME dye, enforced WITHIN each individual
            multiplex pool only — a marker in one multiplex does not
            constrain product sizes used in a different multiplex, even
            on the same dye. If your workflow pools multiple multiplex
            PCRs together before one capillary electrophoresis run, keep
            min_dye_gap generous enough yourself to account for that.
        ta_tolerance: maximum Ta range (°C) allowed within one multiplex
            pool. Default 1.0°C.
        enforce_contig_spread: if True, GRACE tries to select at most one
            marker per contig across the WHOLE panel. This is tried first
            as a strict preference; if it can't fully reach target_count
            using only fresh contigs, a second pass allows reusing a
            contig to fill the remainder rather than leaving the panel
            short.
        preferred_motif_lengths: optional list of motif lengths (e.g. [3,4]
            for tri/tetranucleotides) to prioritise. This is a soft
            priority — if not enough preferred-type loci exist, other
            motif types fill the remainder rather than leaving the panel
            short of the requested count.
        max_per_multiplex: optional cap on how many markers a single
            multiplex pool may contain. A Ta-compatible cluster larger
            than this cap is split across multiple pools rather than
            forced into one oversized multiplex reaction.
        end_dg_threshold: a locus is not placed in a pool if any 3'-end
            stability between one of its primers and a primer already in
            that pool is more negative than this (kcal/mol).
        dimer_dg_threshold: likewise for any heterodimer ΔG (kcal/mol).

    Returns:
        PanelResult with the final marker list, multiplex pool groupings,
        and a clear count of achieved vs requested markers.
    """
    if not dyes:
        raise ValueError("At least one dye must be selected")

    by_locus = _build_candidates(specificity_results)
    if not by_locus:
        return PanelResult(
            markers=[], pools=[], requested_count=target_count, achieved_count=0,
            dyes_used=list(dyes), per_dye_counts={d: 0 for d in dyes},
            unplaced_reason="No PASS-specificity primer pairs are available to build a panel from.",
            end_dg_threshold=end_dg_threshold, dimer_dg_threshold=dimer_dg_threshold,
        )

    flat = [c for cands in by_locus.values() for c in cands]
    ta_clusters = _cluster_by_ta(flat, ta_tolerance)

    # Process clusters with the most distinct candidate loci first, so the
    # panel favours fewer, larger multiplexes over unnecessary fragmentation.
    ta_clusters.sort(key=lambda cluster: -len(set(c.ssr_id for c in cluster)))

    per_dye_target = {
        d: target_count // len(dyes) + (1 if i < target_count % len(dyes) else 0)
        for i, d in enumerate(dyes)
    }
    per_dye_count = {d: 0 for d in dyes}       # GLOBAL count, used only to bias dye choice toward panel-wide balance
    used_contigs = set()
    used_loci = set()
    all_selected = []

    # One list of open pools PER CLUSTER — a pool's per-dye occupied sizes
    # are local to that pool, since size separation now applies only
    # within a single multiplex, not across the whole panel. Different
    # multiplex PCRs are no longer assumed to be pooled together before
    # capillary electrophoresis; if your workflow does pool them, keep
    # min_dye_gap generous enough to allow for that yourself.
    cluster_pools_by_idx = [[] for _ in ta_clusters]
    deferred = []  # (ssr_id, cluster_idx) pairs skipped only due to the contig-spread constraint

    thermo = ThermoAnalysis(**_THERMO_ARGS)
    energy_cache = {}
    rejections = []  # (PrimerInteraction, pool dict) - pool numbers are only known at the end

    def _energies(seq_a, seq_b):
        key = (seq_a, seq_b) if seq_a <= seq_b else (seq_b, seq_a)
        if key not in energy_cache:
            energy_cache[key] = primer_interaction_energies(*key, thermo=thermo)
        return energy_cache[key]

    def _interactions(record_a, record_b, margin):
        """Every (kind, dg, threshold, (label_a, seq_a), (label_b, seq_b))
        between the two pairs' primers that is within `margin` of, or
        beyond, its threshold."""
        found = []
        for pa in _pair_primers(record_a):
            for pb in _pair_primers(record_b):
                dimer, end = _energies(pa[1], pb[1])
                if end < end_dg_threshold + margin:
                    found.append((INTERACTION_END, end, end_dg_threshold, pa, pb))
                if dimer < dimer_dg_threshold + margin:
                    found.append((INTERACTION_DIMER, dimer, dimer_dg_threshold, pa, pb))
        return found

    def _worst_violation(cand, pool):
        """The most severe over-threshold interaction between this
        candidate's primers and the pool's, as (overshoot, PrimerInteraction),
        or None if the candidate is compatible with the whole pool."""
        worst = None
        for m in pool["markers"]:
            for kind, dg, thr, (la, sa), (lb, sb) in _interactions(cand.record, m.record, 0.0):
                if worst is None or dg - thr < worst[0]:
                    worst = (dg - thr, PrimerInteraction(
                        status="REJECTED", kind=kind, dg=dg, threshold=thr, multiplex_pool=-1,
                        ssr_id_a=cand.ssr_id, pair_rank_a=cand.pair_rank, primer_a=la, seq_a=sa,
                        ssr_id_b=m.ssr_id, pair_rank_b=m.pair_rank, primer_b=lb, seq_b=sb,
                    ))
        return worst

    def _new_pool(target_size):
        return {"markers": [], "dye_occupied": {d: [] for d in dyes}, "dye_count": {d: 0 for d in dyes},
                "target_size": target_size}

    def _try_place(ssr_id, candidates_for_locus, contig, pools_for_cluster, size_target_queue):
        """
        Try every existing open pool for this cluster first (first-fit);
        only open a brand-new pool if none of them can take it — either
        because they're all already at their target size, or because
        every dye in every open pool genuinely clashes in size with this
        locus's candidates. Returns True if placed.

        size_target_queue: a shared list acting as a queue of pre-computed
        EVEN pool sizes for this cluster (e.g. [8,8,8,9,9] for 42 markers
        capped at 10) — popped from as each new pool is opened, so pools
        fill up evenly rather than one being greedily filled to the raw
        cap before the next is even started. Falls back to
        max_per_multiplex itself once the precomputed queue is exhausted
        (can happen if size clashes forced more pools than estimated).
        """
        for pool in pools_for_cluster:
            if pool["target_size"] and len(pool["markers"]) >= pool["target_size"]:
                continue
            compatible = []
            worst = None
            for c in candidates_for_locus:
                violation = _worst_violation(c, pool)
                if violation is None:
                    compatible.append(c)
                elif worst is None or violation[0] < worst[0]:
                    worst = violation
            if not compatible:
                rejections.append((worst[1], pool))
                continue
            under_target = [d for d in dyes if per_dye_count[d] < per_dye_target[d]]
            eligible_dyes = under_target if under_target else list(dyes)
            eligible_dyes.sort(key=lambda d: (pool["dye_count"][d], per_dye_count[d]))
            for dye in eligible_dyes:
                occupied = pool["dye_occupied"][dye]
                fitting = [c for c in compatible if _fits_dye(c.product_size, occupied, min_dye_gap)]
                if not fitting:
                    continue
                fitting.sort(key=lambda c: c.pair_rank)
                chosen = fitting[0]
                marker = PanelMarker(
                    ssr_id=chosen.ssr_id, pair_rank=chosen.pair_rank, contig=chosen.contig,
                    motif=chosen.motif, product_size=chosen.product_size, ta=chosen.ta,
                    dye=dye, multiplex_pool=-1, record=chosen.record,
                )
                pool["markers"].append(marker)
                pool["dye_occupied"][dye].append(chosen.product_size)
                pool["dye_count"][dye] += 1
                per_dye_count[dye] += 1
                all_selected.append(marker)
                used_contigs.add(contig)
                used_loci.add(ssr_id)
                return True

        # No existing pool could take it — open a fresh one. An empty
        # pool has no occupied sizes yet, so this always succeeds.
        if size_target_queue:
            target_size = size_target_queue.pop(0)
        else:
            target_size = max_per_multiplex  # overflow beyond the estimate — fall back to the raw cap
        new_pool = _new_pool(target_size)
        under_target = [d for d in dyes if per_dye_count[d] < per_dye_target[d]]
        eligible_dyes = under_target if under_target else list(dyes)
        eligible_dyes.sort(key=lambda d: per_dye_count[d])
        dye = eligible_dyes[0]
        fitting = sorted(candidates_for_locus, key=lambda c: c.pair_rank)
        chosen = fitting[0]
        marker = PanelMarker(
            ssr_id=chosen.ssr_id, pair_rank=chosen.pair_rank, contig=chosen.contig,
            motif=chosen.motif, product_size=chosen.product_size, ta=chosen.ta,
            dye=dye, multiplex_pool=-1, record=chosen.record,
        )
        new_pool["markers"].append(marker)
        new_pool["dye_occupied"][dye].append(chosen.product_size)
        new_pool["dye_count"][dye] += 1
        per_dye_count[dye] += 1
        all_selected.append(marker)
        used_contigs.add(contig)
        used_loci.add(ssr_id)
        pools_for_cluster.append(new_pool)
        return True

    for cluster_idx, cluster in enumerate(ta_clusters):
        if len(all_selected) >= target_count:
            continue

        # ssr_ids are assigned in genome order (contig, start), so sorting
        # by them gives genome order; _spread_order then interleaves that
        # so placement samples the whole genome. Preferred-motif loci still
        # come first, each group spread independently.
        genome_ordered = sorted(set(c.ssr_id for c in cluster if c.ssr_id not in used_loci), key=_sort_key)
        preferred = [lid for lid in genome_ordered
                     if _motif_matches_preference(by_locus[lid][0].canonical_motif, preferred_motif_lengths)]
        preferred_set = set(preferred)
        others = [lid for lid in genome_ordered if lid not in preferred_set]
        locus_ids_in_cluster = _spread_order(preferred) + _spread_order(others)

        # Estimate how many candidates this cluster could plausibly
        # contribute (bounded by what's actually needed to reach the
        # target), and pre-compute EVEN pool-size targets from that
        # estimate — e.g. an estimated 42 with a cap of 10 becomes
        # [8,8,8,9,9], not filling one pool to 10 before starting the
        # next. If real size clashes mean more pools end up genuinely
        # needed than estimated, _try_place() falls back to the raw
        # max_per_multiplex cap for any pools beyond this queue.
        estimate = min(len(locus_ids_in_cluster), target_count - len(all_selected))
        if max_per_multiplex and estimate > 0:
            n_pools_estimate = -(-estimate // max_per_multiplex)  # ceil division
            base, remainder = divmod(estimate, n_pools_estimate)
            size_target_queue = [base + (1 if i < remainder else 0) for i in range(n_pools_estimate)]
        else:
            size_target_queue = []  # no cap — pools are effectively unbounded (target_size=None below)

        for ssr_id in locus_ids_in_cluster:
            if len(all_selected) >= target_count:
                break

            contig = by_locus[ssr_id][0].contig
            if enforce_contig_spread and contig in used_contigs:
                # Don't exclude outright — defer it. If nothing else is
                # available once every cluster has had its first pass,
                # a second pass below will allow same-contig loci rather
                # than leaving the panel short of its target.
                deferred.append((ssr_id, cluster_idx))
                continue

            candidates_for_locus = [c for c in cluster if c.ssr_id == ssr_id]
            _try_place(ssr_id, candidates_for_locus, contig, cluster_pools_by_idx[cluster_idx], size_target_queue)

    # Second pass: only matters when contig-spread is enabled and the
    # target wasn't reached using fresh contigs alone. Revisit deferred
    # same-contig candidates, now allowing them, to fill the remaining
    # slots — spread across contigs as its first choice, but fall back to
    # reusing one rather than leaving the panel short.
    if enforce_contig_spread and len(all_selected) < target_count:
        for ssr_id, cluster_idx in deferred:
            if len(all_selected) >= target_count:
                break
            if ssr_id in used_loci:
                continue
            cluster = ta_clusters[cluster_idx]
            candidates_for_locus = [c for c in cluster if c.ssr_id == ssr_id]
            contig = by_locus[ssr_id][0].contig
            _try_place(ssr_id, candidates_for_locus, contig, cluster_pools_by_idx[cluster_idx], [])

    # Flatten every cluster's pools into the final pool list, numbering
    # them in the order they were actually created.
    pools = []
    for pools_for_cluster in cluster_pools_by_idx:
        for pool in pools_for_cluster:
            if pool["markers"]:
                pools.append(pool["markers"])

    for pool_num, pool in enumerate(pools, start=1):
        for m in pool:
            m.multiplex_pool = pool_num

    # Interactions to report: placements that were rejected, then any
    # pair within each final pool that comes close to (or, defensively,
    # beyond) a threshold.
    interactions = []
    for inter, pool in rejections:
        inter.multiplex_pool = pool["markers"][0].multiplex_pool
        interactions.append(inter)
    for pool_num, pool in enumerate(pools, start=1):
        for i, ma in enumerate(pool):
            for mb in pool[i + 1:]:
                for kind, dg, thr, (la, sa), (lb, sb) in _interactions(ma.record, mb.record, NEAR_THRESHOLD_MARGIN):
                    interactions.append(PrimerInteraction(
                        status="FLAGGED" if dg < thr else "NEAR", kind=kind, dg=dg, threshold=thr,
                        multiplex_pool=pool_num,
                        ssr_id_a=ma.ssr_id, pair_rank_a=ma.pair_rank, primer_a=la, seq_a=sa,
                        ssr_id_b=mb.ssr_id, pair_rank_b=mb.pair_rank, primer_b=lb, seq_b=sb,
                    ))

    achieved = len(all_selected)
    unplaced_reason = None
    if achieved < target_count:
        n_candidate_loci = len(by_locus)
        unplaced_reason = (
            f"Requested {target_count} markers but only {achieved} could be placed. "
            f"{n_candidate_loci} loci had at least one PASS primer pair available. "
            f"This is usually caused by the size-separation requirement (min_dye_gap="
            f"{min_dye_gap}bp, enforced within each multiplex) running out of "
            f"non-overlapping size slots on one or more dyes, and/or the Ta tolerance "
            f"({ta_tolerance}\u00b0C) limiting how many loci can share a multiplex."
        )

    return PanelResult(
        markers=all_selected,
        pools=pools,
        requested_count=target_count,
        achieved_count=achieved,
        dyes_used=list(dyes),
        per_dye_counts=dict(per_dye_count),
        unplaced_reason=unplaced_reason,
        interactions=interactions,
        end_dg_threshold=end_dg_threshold,
        dimer_dg_threshold=dimer_dg_threshold,
    )