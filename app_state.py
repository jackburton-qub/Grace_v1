"""
app_state.py
------------
Single dataclass holding all application data.
Lives for the entire session — no serialisation, no stale keys.
All panels read from and write to this object directly.

Session save/resume: see save_session()/load_session() at the bottom.
Genome sequence data is deliberately NOT embedded in a saved session file
(a multi-Gb genome would make the session file itself multi-Gb and slow
to save/load) — only the original file path is stored, and the genome is
re-loaded from disk when a session is restored. gff_features (a GFFIndex
object, not a plain data structure) is handled the same way: rebuilt from
gff_path on restore rather than serialised directly.
"""

import json
import os
from dataclasses import dataclass, asdict
from typing import Optional


@dataclass
class AppState:
    # ── Genome ────────────────────────────────────────────
    genome:          Optional[dict]  = None   # {contig_name: sequence}
    genome_path:     Optional[str]   = None   # path to FASTA file on disk
    genome_filename: Optional[str]   = None

    # ── GFF Annotation (optional) ─────────────────────────
    gff_path:        Optional[str]    = None   # path to GFF/GTF file on disk
    gff_filename:    Optional[str]    = None
    gff_features:    Optional[object] = None   # parsed GFFIndex (built lazily)
    chrom_names:     Optional[dict]   = None   # contig accession → chromosome display name

    # ── SSR Detection ─────────────────────────────────────
    ssrs:               Optional[list]  = None  # list of SSR dicts
    selected_ssrs:      Optional[list]  = None  # user-selected subset
    ssr_detection_time: Optional[float] = None

    # ── Primer Design ─────────────────────────────────────
    primer_results:          Optional[list] = None  # all designed primers
    filtered_primer_results: Optional[list] = None  # after quality filters
    # Run summary for the Primer Design page: {"n_attempted", "elapsed",
    # "mode", "issues": [one dict per locus that failed or was size-limited]}
    primer_design_summary:   Optional[dict] = None
    product_min:             int  = 100
    product_max:             int  = 250

    # Existing primer pairs imported from a file (core/primer_import.py),
    # screened alongside any designed pairs. Independent of the genome and
    # SSRs, so loading a new genome or re-running design keeps them.
    imported_primers:        Optional[list] = None

    # Workflow mode: "capillary", "amplicon", or "gbs_re"
    workflow_mode:           str  = "capillary"

    # ── Specificity ───────────────────────────────────────
    specificity_results: Optional[list] = None
    blast_raw_hits_path: Optional[str] = None

    # ── Panel Design ───────────────────────────────────────
    panel_design_result: Optional[object] = None

    # ── Run version stamps ────────────────────────────────
    # Panels cache the version they last rendered (_last_rendered_version)
    # and skip re-rendering when the current value matches — this is how
    # navigating between panels without new data avoids needless rebuilds.
    # IMPORTANT: these must only ever INCREMENT, never be reset to 0 or
    # any previously-seen value. A reset can make a future run's version
    # number coincide with a value some panel already has cached from an
    # earlier cycle, causing that panel to wrongly conclude nothing
    # changed and silently skip a legitimate re-render (only catching up
    # on a SUBSEQUENT run, once the number reaches something genuinely
    # new) — this was a real, confirmed bug. Clearing state still counts
    # as a real change worth a bump, not just a fresh run producing data.
    ssr_version:     int = 0
    primer_version:  int = 0
    blast_version:   int = 0

    # ---------------------------------------------------------
    # Clear helpers
    # ---------------------------------------------------------

    def clear_downstream_of_genome(self):
        """Call when a new genome is loaded."""
        self.ssrs               = None
        self.selected_ssrs      = None
        self.ssr_detection_time = None
        self.ssr_version        += 1  # bump, don't reset — see module docstring
        self.clear_downstream_of_ssrs()

    def clear_gff(self):
        """Call when GFF is cleared or a new one loaded."""
        self.gff_path     = None
        self.gff_filename = None
        self.gff_features = None
        self.chrom_names  = None

    def clear_downstream_of_ssrs(self):
        """Call when new SSRs are detected."""
        self.primer_results          = None
        self.filtered_primer_results = None
        self.primer_design_summary   = None
        self.primer_version          += 1  # bump, don't reset — see module docstring
        self.clear_downstream_of_primers()

    def clear_downstream_of_primers(self):
        """Call when new primers are designed."""
        self.specificity_results = None
        self.blast_raw_hits_path = None
        self.blast_version       += 1  # bump, don't reset — see module docstring
        self.panel_design_result = None

    # ---------------------------------------------------------
    # Properties
    # ---------------------------------------------------------

    @property
    def has_genome(self) -> bool:
        return self.genome is not None

    @property
    def has_gff(self) -> bool:
        return self.gff_path is not None

    @property
    def has_ssrs(self) -> bool:
        return self.ssrs is not None and len(self.ssrs) > 0

    @property
    def has_primers(self) -> bool:
        return self.primer_results is not None and len(self.primer_results) > 0

    @property
    def has_specificity(self) -> bool:
        return (self.specificity_results is not None and
                len(self.specificity_results) > 0)

    @property
    def has_panel_design(self) -> bool:
        return self.panel_design_result is not None

    @property
    def has_imported_primers(self) -> bool:
        return bool(self.imported_primers)

    @property
    def blast_primers(self) -> list:
        """Designed primers (filtered if available, otherwise full set),
        followed by any imported pairs."""
        if self.filtered_primer_results:
            designed = self.filtered_primer_results
        else:
            designed = self.primer_results or []
        return list(designed) + list(self.imported_primers or [])

    def get_display_name(self, contig: str) -> str:
        """Return chromosome name if available, otherwise contig."""
        if self.chrom_names and contig in self.chrom_names:
            return self.chrom_names[contig]
        if self.gff_features and hasattr(self.gff_features, 'get_display_name'):
            return self.gff_features.get_display_name(contig)
        return contig

    # ---------------------------------------------------------
    # Session save/resume
    # ---------------------------------------------------------

    # Fields excluded from the saved session file — either too large to
    # embed (genome), rebuilt on restore rather than serialised directly
    # (gff_features, a GFFIndex object graph, not plain data), or
    # transient run bookkeeping that doesn't need to survive a save
    # (version stamps reset naturally as each stage re-runs anyway).
    _SESSION_EXCLUDED_FIELDS = {"genome", "gff_features"}

    def to_session_dict(self) -> dict:
        """
        Build a JSON-safe dict representing this session, excluding the
        genome sequence data itself and the rebuildable gff_features index.
        """
        data = {}
        for key, value in asdict_shallow(self).items():
            if key in self._SESSION_EXCLUDED_FIELDS:
                continue
            if key == "panel_design_result" and value is not None:
                data[key] = asdict(value)  # PanelResult is a real dataclass — recurse properly
            else:
                data[key] = value
        return data

    def save_session(self, path: str) -> None:
        """Write this session to a JSON file at `path`."""
        with open(path, "w", encoding="utf-8") as f:
            json.dump(self.to_session_dict(), f, indent=2, cls=_SessionJSONEncoder)

    @classmethod
    def load_session(cls, path: str) -> "AppState":
        """
        Load a session from a JSON file into a brand new AppState instance.
        Prefer load_session_into() when updating an AppState that other
        objects already hold a reference to (e.g. every panel in the GUI) —
        this classmethod is for cases where a fresh, standalone instance
        is genuinely what's wanted (e.g. the CLI, or tests).
        """
        state = cls()
        state.load_session_into(path)
        return state

    def load_session_into(self, path: str) -> None:
        """
        Load a session from a JSON file, updating THIS instance's fields
        in place rather than returning a new object. This is the version
        the GUI must use: main_window and every panel are constructed with
        a reference to one shared AppState — replacing that reference
        with a new object would leave every panel still pointing at the
        old, stale one.

        Re-loads the genome from its original file path (raises
        FileNotFoundError with a clear message if that file no longer
        exists or has moved — the caller should catch this and prompt the
        user to relocate it) and rebuilds gff_features from gff_path if
        one was set.
        """
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)

        for key, value in data.items():
            if key == "panel_design_result" and value is not None:
                value = _panel_result_from_dict(value)
            if hasattr(self, key):
                setattr(self, key, value)

        if self.genome_path:
            if not os.path.exists(self.genome_path):
                raise FileNotFoundError(
                    f"Genome file from this session is no longer at its original "
                    f"location: {self.genome_path}\n"
                    f"Please re-load the genome manually, then reload the rest of "
                    f"this session's results."
                )
            from core.fasta_loader import load_sequence_file
            self.genome = load_sequence_file(self.genome_path)

        if self.gff_path and os.path.exists(self.gff_path) and self.genome:
            from core.gff_parser import build_gff_index
            self.gff_features = build_gff_index(self.gff_path, genome=self.genome)


class _SessionJSONEncoder(json.JSONEncoder):
    """
    Defensive encoder for session saving. Every AppState field is SUPPOSED
    to hold plain dicts/lists/scalars by the time it gets here, but a
    pandas DataFrame (or a numpy scalar/array, which pandas often produces
    even inside otherwise-plain data) ending up in a field unexpectedly
    would otherwise crash the ENTIRE save with an unhelpful "Object of
    type X is not JSON serializable" and no indication of which field was
    responsible. This converts the common offenders correctly — regardless
    of how deeply nested they are in the overall structure, since
    json.dump() calls default() on whatever specific object it can't
    serialize natively — and for anything genuinely unrecognised, raises a
    clearer error naming the actual type involved instead of a bare crash.
    """
    def default(self, obj):
        try:
            import pandas as pd
            if isinstance(obj, pd.DataFrame):
                return obj.to_dict(orient="records")
            if isinstance(obj, pd.Series):
                return obj.tolist()
        except ImportError:
            pass
        try:
            import numpy as np
            if isinstance(obj, np.integer):
                return int(obj)
            if isinstance(obj, np.floating):
                return float(obj)
            if isinstance(obj, np.ndarray):
                return obj.tolist()
        except ImportError:
            pass
        raise TypeError(
            f"Session save failed: a field contains a {type(obj).__name__} "
            f"object, which isn't JSON-serializable and isn't one of the "
            f"types already handled here (DataFrame, Series, numpy scalar/"
            f"array). This usually means something stored a non-plain-data "
            f"object into AppState — please report this, including what "
            f"you were doing right before saving."
        )


def asdict_shallow(state: "AppState") -> dict:
    """
    Like dataclasses.asdict(), but does NOT recurse into nested
    dataclasses/objects for every field — only plain dicts/lists/scalars
    are needed for most AppState fields, and a full recursive asdict()
    would need special-casing for gff_features (a GFFIndex, not a
    dataclass, which asdict() cannot handle at all and would raise on).
    panel_design_result is handled separately by the caller, since it IS
    a real dataclass that does need proper recursive conversion.
    """
    from dataclasses import fields
    return {f.name: getattr(state, f.name) for f in fields(state)}


def _panel_result_from_dict(d: dict):
    """Reconstruct a PanelResult (with nested PanelMarker objects) from
    its plain-dict form, as produced by dataclasses.asdict()."""
    from core.panel_design import PanelResult, PanelMarker, PrimerInteraction
    markers = [PanelMarker(**m) for m in d.get("markers", [])]
    pools = [[PanelMarker(**m) for m in pool] for pool in d.get("pools", [])]
    return PanelResult(
        markers=markers,
        pools=pools,
        requested_count=d["requested_count"],
        achieved_count=d["achieved_count"],
        dyes_used=d["dyes_used"],
        per_dye_counts=d["per_dye_counts"],
        unplaced_reason=d.get("unplaced_reason"),
        interactions=[PrimerInteraction(**i) for i in d.get("interactions", [])],
        end_dg_threshold=d.get("end_dg_threshold"),
        dimer_dg_threshold=d.get("dimer_dg_threshold"),
    )