"""
panel_design_panel.py — Panel Design
Builds an optimised multiplex genotyping panel from PASS-specificity
primer pairs: dye assignment, Ta-based multiplex grouping, and same-dye
product-size separation to avoid fragment overlap.
"""
import time

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QGridLayout, QSpinBox, QDoubleSpinBox, QCheckBox,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView,
    QScrollArea, QFileDialog,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QFont

from ui.style import (
    ACCENT, SUCCESS, WARNING, TEXT_SECONDARY,
    FONT_UI, FONT_SIZE_LARGE, FONT_SIZE_SMALL, PANEL_PADDING,
    BG_MID, BORDER,
)
from ui.table_utils import NumericTableWidgetItem, secondary_label as _lbl
from core.panel_design import DEFAULT_END_DG_THRESHOLD, DEFAULT_DIMER_DG_THRESHOLD


class PanelWorker(QThread):
    # Primer-interaction checks make panel design take seconds rather than
    # milliseconds on larger candidate sets, so it runs off the UI thread.
    finished = pyqtSignal(object, float)
    error    = pyqtSignal(str)

    def __init__(self, specificity_results, params):
        super().__init__()
        self.specificity_results = specificity_results
        self.params = params

    def run(self):
        try:
            from core.panel_design import design_panel
            start = time.time()
            result = design_panel(self.specificity_results, **self.params)
            self.finished.emit(result, time.time() - start)
        except Exception as e:
            self.error.emit(str(e))


def _primer_label(ssr_id, pair_rank, primer):
    return f"SSR {ssr_id} (rank {pair_rank}) {primer}"


def _interaction_notes_by_marker(result) -> dict:
    """
    {(ssr_id, pair_rank): [note, ...]} for the CSV export. Near/flagged
    interactions inside a final pool are noted on both markers involved;
    a rejected placement is noted on the locus that was moved.
    """
    by_ssr = {m.ssr_id: (m.ssr_id, m.pair_rank) for m in result.markers}
    notes = {}
    for i in getattr(result, "interactions", None) or []:
        desc = f"{i.kind} {i.dg:.2f} kcal/mol (limit {i.threshold:.1f})"
        if i.status == "REJECTED":
            key = by_ssr.get(i.ssr_id_a)
            if key:
                notes.setdefault(key, []).append(
                    f"REJECTED from multiplex {i.multiplex_pool}: {i.primer_a} with "
                    f"{_primer_label(i.ssr_id_b, i.pair_rank_b, i.primer_b)}, {desc}"
                )
            continue
        notes.setdefault((i.ssr_id_a, i.pair_rank_a), []).append(
            f"{i.status}: {i.primer_a} with {_primer_label(i.ssr_id_b, i.pair_rank_b, i.primer_b)}, {desc}"
        )
        notes.setdefault((i.ssr_id_b, i.pair_rank_b), []).append(
            f"{i.status}: {i.primer_b} with {_primer_label(i.ssr_id_a, i.pair_rank_a, i.primer_a)}, {desc}"
        )
    return notes


class PanelDesignPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state = state
        self.mw    = main_window
        self._last_result = None
        self._worker = None
        self._build_ui()

    def _build_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        outer.addWidget(scroll)

        content = QWidget()
        scroll.setWidget(content)
        L = QVBoxLayout(content)
        L.setContentsMargins(PANEL_PADDING, PANEL_PADDING, PANEL_PADDING, PANEL_PADDING)
        L.setSpacing(16)

        title = QLabel("Panel Design")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setObjectName("role_title")
        sub = QLabel(
            "Build an optimised multiplex genotyping panel from your BLAST-validated primers — "
            "assigns dyes, groups compatible markers into multiplex PCRs, and keeps same-dye "
            "products apart in size to avoid fragment overlap."
        )
        sub.setWordWrap(True)
        sub.setObjectName("role_secondary")
        sub.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(title)
        L.addWidget(sub)

        # ---- Settings ----
        settings_group = QGroupBox("Panel Settings")
        sg = QGridLayout(settings_group)
        sg.setSpacing(10)

        sg.addWidget(_lbl(
            "Dyes to use",
            "Which fluorescent dye channels are available on your capillary "
            "electrophoresis setup. The panel spreads markers as evenly as "
            "possible across whichever dyes you tick."
        ), 0, 0)
        dye_row = QHBoxLayout()
        self.dye_checks = {}
        for dye in ["FAM", "HEX", "NED", "PET"]:
            cb = QCheckBox(dye)
            cb.setChecked(True)
            self.dye_checks[dye] = cb
            dye_row.addWidget(cb)
        dye_row.addStretch()
        dye_widget = QWidget()
        dye_widget.setLayout(dye_row)
        sg.addWidget(dye_widget, 0, 1, 1, 3)

        sg.addWidget(_lbl(
            "Target marker count",
            "Total number of SSRs you want in the final panel, across all "
            "multiplex PCRs combined. If this can't be fully achieved given "
            "your other settings, the best possible count is reported back "
            "with an explanation."
        ), 1, 0)
        self.target_count = QSpinBox()
        self.target_count.setRange(1, 500)
        self.target_count.setValue(20)
        sg.addWidget(self.target_count, 1, 1)

        sg.addWidget(_lbl(
            "Min size gap on same dye (bp)",
            "Minimum product-size difference required between two markers "
            "sharing the same dye channel, enforced across the WHOLE panel "
            "(not just within one multiplex) — different multiplex PCRs are "
            "commonly pooled together for one capillary run, so this needs "
            "to hold regardless of which PCR a marker came from."
        ), 1, 2)
        self.min_dye_gap = QSpinBox()
        self.min_dye_gap.setRange(1, 200)
        self.min_dye_gap.setValue(75)
        sg.addWidget(self.min_dye_gap, 1, 3)

        sg.addWidget(_lbl(
            "Max primers per multiplex",
            "Caps how many markers a single multiplex PCR reaction may contain. "
            "If a Ta-compatible group of loci is larger than this cap, it gets "
            "split across multiple multiplex pools rather than forced into one "
            "oversized reaction. Leave at 0 for no cap."
        ), 4, 0)
        self.max_per_multiplex = QSpinBox()
        self.max_per_multiplex.setRange(0, 200)
        self.max_per_multiplex.setValue(0)
        self.max_per_multiplex.setSpecialValueText("No cap")
        sg.addWidget(self.max_per_multiplex, 4, 1)

        sg.addWidget(_lbl(
            "Max 3\u2032-end interaction \u0394G (kcal/mol)",
            "A locus is kept out of a multiplex if either of its primers forms "
            "a 3\u2032-anchored interaction with any primer already in that "
            "multiplex more stable (more negative) than this. 3\u2032-anchored "
            "interactions can be extended by the polymerase into primer-dimer "
            "artefacts, so this is the stricter of the two limits. Calculated "
            "with Primer3's end-stability model at the same salt/dNTP "
            "concentrations used for primer design. The locus is tried in the "
            "next multiplex instead, or gets a new one."
        ), 4, 2)
        self.end_dg_threshold = QDoubleSpinBox()
        self.end_dg_threshold.setRange(-30.0, 0.0)
        self.end_dg_threshold.setSingleStep(0.5)
        self.end_dg_threshold.setDecimals(1)
        self.end_dg_threshold.setValue(DEFAULT_END_DG_THRESHOLD)
        sg.addWidget(self.end_dg_threshold, 4, 3)

        sg.addWidget(_lbl(
            "Max primer-dimer \u0394G (kcal/mol)",
            "A locus is kept out of a multiplex if either of its primers forms "
            "a dimer anywhere along its length with any primer already in that "
            "multiplex more stable (more negative) than this. Calculated with "
            "Primer3's heterodimer model at the same salt/dNTP concentrations "
            "used for primer design. Interactions within 1 kcal/mol of either "
            "limit are listed in the results so you can review them."
        ), 5, 0)
        self.dimer_dg_threshold = QDoubleSpinBox()
        self.dimer_dg_threshold.setRange(-30.0, 0.0)
        self.dimer_dg_threshold.setSingleStep(0.5)
        self.dimer_dg_threshold.setDecimals(1)
        self.dimer_dg_threshold.setValue(DEFAULT_DIMER_DG_THRESHOLD)
        sg.addWidget(self.dimer_dg_threshold, 5, 1)

        sg.addWidget(_lbl(
            "Max Ta range per multiplex (\u00b0C)",
            "Maximum annealing-temperature spread tolerated within a single "
            "multiplex PCR reaction. Loci are grouped into Ta-compatible "
            "pools based on this tolerance BEFORE markers are selected, "
            "then each pool is filled independently — a wider range means "
            "fewer, larger multiplexes; a tighter range means more, "
            "smaller ones. Your target marker count may end up split "
            "across several multiplexes if your PASS loci don't all share "
            "one narrow Ta band — e.g. 20 requested markers might come "
            "back as two pools of 10 and 10, or 9 and 11, depending on "
            "what's actually available."
        ), 2, 0)
        self.ta_tolerance = QDoubleSpinBox()
        self.ta_tolerance.setRange(0.5, 15.0)
        self.ta_tolerance.setSingleStep(0.5)
        self.ta_tolerance.setValue(1.0)
        sg.addWidget(self.ta_tolerance, 2, 1)

        self.contig_spread_cb = QCheckBox("Spread markers across different contigs")
        self.contig_spread_cb.setChecked(True)
        self.contig_spread_cb.setToolTip(
            "When enabled, GRACE tries to select at most one marker per contig — "
            "recommended for population genetics work where marker independence "
            "matters, since two SSRs on the same contig/chromosome aren't "
            "independently inherited. This is a preference, not a hard cutoff: "
            "if spreading across contigs can't fully reach your target marker "
            "count, GRACE falls back to reusing a contig rather than leaving "
            "the panel short — spreading is always tried first."
        )
        sg.addWidget(self.contig_spread_cb, 2, 2, 1, 2)

        sg.addWidget(_lbl(
            "Prioritise motif types",
            "Loci matching a ticked motif length are always considered first "
            "when filling the panel. This is a soft priority, not a hard "
            "filter — if there aren't enough of your preferred type(s) to "
            "reach the target marker count, other motif types fill the "
            "remainder rather than leaving the panel short. Leave all "
            "unticked (or tick none specifically) for no preference."
        ), 3, 0)
        motif_row = QHBoxLayout()
        self.motif_checks = {}
        for label, length in [("Mono", 1), ("Di", 2), ("Tri", 3), ("Tetra", 4), ("Penta", 5), ("Hexa", 6)]:
            cb = QCheckBox(label)
            cb.setChecked(False)
            self.motif_checks[length] = cb
            motif_row.addWidget(cb)
        motif_row.addStretch()
        motif_widget = QWidget()
        motif_widget.setLayout(motif_row)
        sg.addWidget(motif_widget, 3, 1, 1, 3)

        L.addWidget(settings_group)

        # ---- Run ----
        run_row = QHBoxLayout()
        self.run_btn = QPushButton("Generate Panel")
        self.run_btn.setObjectName("primary")
        self.run_btn.setFixedHeight(40)
        self.run_btn.clicked.connect(self._run)
        run_row.addWidget(self.run_btn)
        run_row.addStretch()
        L.addLayout(run_row)

        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        self.status_label.setObjectName("role_secondary")
        self.status_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(self.status_label)

        # ---- Results ----
        self.results_group = QGroupBox("Panel Results")
        self.results_group.setVisible(False)
        rl = QVBoxLayout(self.results_group)

        self.summary_label = QLabel("")
        self.summary_label.setWordWrap(True)
        self.summary_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        rl.addWidget(self.summary_label)

        self.table = QTableWidget()
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSortingEnabled(True)
        self.table.setMinimumHeight(400)
        rl.addWidget(self.table)

        self.interaction_label = QLabel("")
        self.interaction_label.setWordWrap(True)
        self.interaction_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        rl.addWidget(self.interaction_label)

        self.interaction_table = QTableWidget()
        self.interaction_table.setAlternatingRowColors(True)
        self.interaction_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.interaction_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.interaction_table.setSortingEnabled(True)
        self.interaction_table.setMinimumHeight(200)
        rl.addWidget(self.interaction_table)

        dl_row = QHBoxLayout()
        self.download_btn = QPushButton("Download panel (CSV)")
        self.download_btn.clicked.connect(self._download_csv)
        self.download_flanking_btn = QPushButton("Export Flanking Sequences (HTML)")
        self.download_flanking_btn.setToolTip(
            "For every marker in this panel, exports the primer + SSR region plus "
            "~10bp of flanking sequence on each side, as an HTML file with the "
            "primers and SSR region highlighted in colour for easy visual inspection."
        )
        self.download_flanking_btn.clicked.connect(self._download_flanking_html)
        self.download_flanking_fasta_btn = QPushButton("Export Flanking Sequences (FASTA)")
        self.download_flanking_fasta_btn.setToolTip(
            "Same region as the HTML export, as plain-text FASTA for pasting into "
            "other tools. Primer and SSR positions are given as numbers in the "
            "header line, since plain text can't show colour."
        )
        self.download_flanking_fasta_btn.clicked.connect(self._download_flanking_fasta)
        dl_row.addWidget(self.download_btn)
        dl_row.addWidget(self.download_flanking_btn)
        dl_row.addWidget(self.download_flanking_fasta_btn)
        dl_row.addStretch()
        rl.addLayout(dl_row)

        L.addWidget(self.results_group)
        L.addStretch()

    def _run(self):
        if self._worker is not None:
            return  # already running - ignore duplicate clicks
        if not self.state.has_specificity:
            self.status_label.setText("Run BLAST specificity checking first.")
            self.status_label.setObjectName("role_warning")
            self.status_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
            return

        dyes = [d for d, cb in self.dye_checks.items() if cb.isChecked()]
        if not dyes:
            self.status_label.setText("Select at least one dye.")
            self.status_label.setObjectName("role_warning")
            self.status_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
            return

        params = dict(
            dyes=dyes,
            target_count=self.target_count.value(),
            min_dye_gap=self.min_dye_gap.value(),
            ta_tolerance=self.ta_tolerance.value(),
            enforce_contig_spread=self.contig_spread_cb.isChecked(),
            preferred_motif_lengths=[length for length, cb in self.motif_checks.items() if cb.isChecked()] or None,
            max_per_multiplex=self.max_per_multiplex.value() or None,
            end_dg_threshold=self.end_dg_threshold.value(),
            dimer_dg_threshold=self.dimer_dg_threshold.value(),
        )
        self.run_btn.setEnabled(False)
        self.status_label.setText("Building panel and checking primer-primer interactions...")
        self.status_label.setObjectName("role_secondary")
        self.status_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        self.mw.set_status("Running panel design...")
        self._worker = PanelWorker(self.state.specificity_results, params)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _finish_worker(self):
        self.run_btn.setEnabled(True)
        if self._worker is not None:
            self._worker.wait()
            self._worker.deleteLater()
            self._worker = None

    def _on_error(self, msg):
        self._finish_worker()
        self.status_label.setText(f"Error: {msg}")
        self.status_label.setObjectName("role_error")
        self.status_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        self.mw.set_status("Panel design failed")

    def _on_done(self, result, elapsed):
        self._finish_worker()
        self._last_result = result
        self.state.panel_design_result = result

        if result.achieved_count < result.requested_count:
            self.status_label.setText(
                f"Achieved {result.achieved_count} of {result.requested_count} requested markers "
                f"in {elapsed:.1f}s."
            )
            self.status_label.setObjectName("role_warning")
            self.status_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        else:
            self.status_label.setText(
                f"Achieved all {result.achieved_count} requested markers in {elapsed:.1f}s."
            )
            self.status_label.setObjectName("role_success")
            self.status_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")

        self.mw.set_status(f"Panel design complete — {result.achieved_count} markers")
        self.mw.on_step_complete(6)
        self._populate_results(result)

    def _populate_results(self, result):
        self.results_group.setVisible(True)

        dye_summary = "  |  ".join(f"{d}: {c}" for d, c in result.per_dye_counts.items())
        summary = (
            f"<b>{result.achieved_count} of {result.requested_count} markers</b> across "
            f"<b>{len(result.pools)} multiplex pool(s)</b><br>"
            f"Per-dye counts — {dye_summary}"
        )
        if result.unplaced_reason:
            summary += f'<br><span style="color:{WARNING};">{result.unplaced_reason}</span>'
        self.summary_label.setText(summary)

        has_imported = any(m.record.get("source") == "imported" for m in result.markers)
        cols = ["multiplex_pool", "dye", "ssr_id"] + (["source"] if has_imported else []) + [
                "pair_rank", "contig", "motif", "product_size", "ta"]
        col_labels = {
            "multiplex_pool": "Multiplex",
            "dye":            "Dye",
            "ssr_id":         "SSR ID",
            "source":         "Source",
            "pair_rank":      "Pair rank",
            "contig":         "Contig",
            "motif":          "Motif",
            "product_size":   "Product size (bp)",
            "ta":             "Ta (\u00b0C)",
        }

        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(result.markers))
        self.table.setColumnCount(len(cols))
        self.table.setHorizontalHeaderLabels([col_labels[c] for c in cols])

        # Sort for display: by multiplex pool, then dye, then product size
        markers_sorted = sorted(result.markers, key=lambda m: (m.multiplex_pool, m.dye, m.product_size))

        for row_idx, m in enumerate(markers_sorted):
            values = {
                "multiplex_pool": m.multiplex_pool,
                "dye":            m.dye,
                "ssr_id":         m.ssr_id,
                "source":         m.record.get("source") or "designed",
                "pair_rank":      m.pair_rank,
                "contig":         m.contig,
                "motif":          m.motif,
                "product_size":   m.product_size,
                "ta":             round(m.ta, 2),
            }
            for col_idx, col in enumerate(cols):
                val = values[col]
                if isinstance(val, (int, float)):
                    item = NumericTableWidgetItem(str(val), numeric_value=val)
                else:
                    item = QTableWidgetItem(str(val))
                self.table.setItem(row_idx, col_idx, item)

        header = self.table.horizontalHeader()
        header.setStretchLastSection(True)
        for col_idx in range(len(cols)):
            header.setSectionResizeMode(col_idx, QHeaderView.ResizeMode.ResizeToContents)

        self.table.setSortingEnabled(True)
        self._populate_interactions(result)

    def _populate_interactions(self, result):
        interactions = getattr(result, "interactions", None) or []
        n_rejected = sum(1 for i in interactions if i.status == "REJECTED")
        n_flagged = sum(1 for i in interactions if i.status == "FLAGGED")
        n_near = sum(1 for i in interactions if i.status == "NEAR")

        if result.end_dg_threshold is None:
            # Panel from an older session, built before interaction checks existed.
            self.interaction_label.setText(
                "<b>Primer interactions</b> \u2014 not checked for this panel. "
                "Regenerate it to check primer-primer interactions."
            )
        elif not interactions:
            self.interaction_label.setText(
                f"<b>Primer interactions</b> \u2014 none at or near the limits "
                f"(3\u2032 end {result.end_dg_threshold:.1f}, dimer "
                f"{result.dimer_dg_threshold:.1f} kcal/mol)."
            )
        else:
            parts = []
            if n_rejected:
                parts.append(f"{n_rejected} placement(s) rejected and moved to another multiplex")
            if n_flagged:
                parts.append(f'<span style="color:{WARNING};">{n_flagged} over-threshold interaction(s) in final pools</span>')
            if n_near:
                parts.append(f"{n_near} near-threshold interaction(s) in final pools")
            self.interaction_label.setText(
                f"<b>Primer interactions</b> (3\u2032 end {result.end_dg_threshold:.1f}, dimer "
                f"{result.dimer_dg_threshold:.1f} kcal/mol) \u2014 " + "; ".join(parts) + "."
            )

        self.interaction_table.setVisible(bool(interactions))
        cols = ["Status", "Multiplex", "Type", "\u0394G (kcal/mol)", "Threshold",
                "Primer A", "Sequence A", "Primer B", "Sequence B"]
        self.interaction_table.setSortingEnabled(False)
        self.interaction_table.setRowCount(len(interactions))
        self.interaction_table.setColumnCount(len(cols))
        self.interaction_table.setHorizontalHeaderLabels(cols)
        for row_idx, i in enumerate(interactions):
            values = [
                i.status, i.multiplex_pool, i.kind, round(i.dg, 2), i.threshold,
                _primer_label(i.ssr_id_a, i.pair_rank_a, i.primer_a), i.seq_a,
                _primer_label(i.ssr_id_b, i.pair_rank_b, i.primer_b), i.seq_b,
            ]
            for col_idx, val in enumerate(values):
                if isinstance(val, (int, float)):
                    item = NumericTableWidgetItem(str(val), numeric_value=val)
                else:
                    item = QTableWidgetItem(str(val))
                self.interaction_table.setItem(row_idx, col_idx, item)
        header = self.interaction_table.horizontalHeader()
        header.setStretchLastSection(True)
        for col_idx in range(len(cols)):
            header.setSectionResizeMode(col_idx, QHeaderView.ResizeMode.ResizeToContents)
        self.interaction_table.setSortingEnabled(True)

    def _download_csv(self):
        if not self._last_result or not self._last_result.markers:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save panel", "genotyping_panel.csv", "CSV files (*.csv)")
        if not path:
            return
        import pandas as pd
        notes = _interaction_notes_by_marker(self._last_result)
        has_imported = any(m.record.get("source") == "imported" for m in self._last_result.markers)
        rows = []
        for m in sorted(self._last_result.markers, key=lambda m: (m.multiplex_pool, m.dye, m.product_size)):
            source = {"Source": m.record.get("source") or "designed"} if has_imported else {}
            rows.append({
                "Multiplex": m.multiplex_pool,
                "Dye": m.dye,
                "SSR ID": m.ssr_id,
                **source,
                "Pair rank": m.pair_rank,
                "Contig": m.contig,
                "Motif": m.motif,
                "Product size (bp)": m.product_size,
                "Ta (\u00b0C)": round(m.ta, 2),
                "Forward primer": m.record.get("left_primer", ""),
                "Reverse primer": m.record.get("right_primer", ""),
                "Primer interactions": "; ".join(notes.get((m.ssr_id, m.pair_rank), [])),
            })
        pd.DataFrame(rows).to_csv(path, index=False, encoding="utf-8-sig")
        self.mw.set_status(f"Saved to {path}")

    def _download_flanking_html(self):
        if not self._last_result or not self._last_result.markers:
            self.mw.set_status("No panel available to export")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save flanking sequences", "panel_flanking_sequences.html", "HTML files (*.html)"
        )
        if not path:
            return
        from core.primer_design import build_flanking_sequence_html
        records = [m.record for m in self._last_result.markers]
        html = build_flanking_sequence_html(self.state.genome, records, flank_bp=10)
        with open(path, "w", encoding="utf-8") as f:
            f.write(html)
        self.mw.set_status(f"Saved to {path}")

    def _download_flanking_fasta(self):
        if not self._last_result or not self._last_result.markers:
            self.mw.set_status("No panel available to export")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save flanking sequences", "panel_flanking_sequences.fasta", "FASTA files (*.fasta *.fa)"
        )
        if not path:
            return
        from core.primer_design import build_flanking_sequence_fasta
        records = [m.record for m in self._last_result.markers]
        fasta = build_flanking_sequence_fasta(self.state.genome, records, flank_bp=10)
        with open(path, "w", encoding="utf-8") as f:
            f.write(fasta)
        self.mw.set_status(f"Saved to {path}")

    def on_show(self):
        if self.state.panel_design_result and self._last_result is not self.state.panel_design_result:
            self._last_result = self.state.panel_design_result
            self._populate_results(self._last_result)