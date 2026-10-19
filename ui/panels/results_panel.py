"""
results_panel.py — BLAST Results
Pure display. Single scroll area. Reads state only.
"""
import os
import shutil
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QGridLayout, QTableWidget, QTableWidgetItem,
    QAbstractItemView, QScrollArea, QFileDialog,
    QCheckBox,
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QColor

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY,
    FONT_UI, FONT_MONO, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING,
)
from ui.table_utils import NumericTableWidgetItem, TABLE_DISPLAY_LIMIT, secondary_label as _lbl



def _clean_df(df):
    """Remove nested/messy columns before CSV export."""
    import pandas as pd
    exclude = {"amplicons", "left_hits", "right_hits", "offtarget_amplicons"}
    cols = [c for c in df.columns if c not in exclude]
    out = df[cols].copy()
    # Rename to human-readable
    rename = {
        "ssr_id": "SSR ID", "pair_rank": "Pair rank",
        "specificity_status": "Status", "contig": "Contig",
        "motif": "Motif", "canonical_motif": "Canonical motif",
        "repeat_count": "Repeat count", "start": "Start (bp)", "end": "End (bp)",
        "left_primer": "Forward primer", "right_primer": "Reverse primer",
        "product_size": "Product size (bp)",
        "left_tm": "Forward Tm (C)", "right_tm": "Reverse Tm (C)",
        "left_gc": "Forward GC (%)", "right_gc": "Reverse GC (%)",
        "left_3end_dg": "Forward 3' stability (kcal/mol)",
        "right_3end_dg": "Reverse 3' stability (kcal/mol)",
        "left_gc_clamp": "Forward GC clamp", "right_gc_clamp": "Reverse GC clamp",
    }
    return out.rename(columns={k: v for k, v in rename.items() if k in out.columns})

def _label_sources(df):
    """If any imported pairs are present, fill the source column so every
    row says whether it was designed or imported. Results with no imported
    pairs have no source column and are left exactly as before."""
    if "source" in df.columns:
        df = df.copy()
        df["source"] = df["source"].fillna("designed")
    return df


class ResultsPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state = state
        self.mw    = main_window
        self._last_rendered_version = -1
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

        # Title
        title = QLabel("BLAST Results")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setObjectName("role_title")
        sub = QLabel("Specificity check results — primers verified against the genome")
        sub.setObjectName("role_secondary")
        sub.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(title); L.addWidget(sub)

        # Summary
        self.summary_group = QGroupBox("Summary")
        sg = QGridLayout(self.summary_group)
        sg.setSpacing(12)
        self.metric_labels = {}
        metric_tooltips = {
            "Unique SSRs": (
                "Counts distinct SSR loci with at least one PASS primer pair, "
                "not just the total number of passing pairs — a locus with "
                "more than one ranked pair that each independently pass is "
                "still one unique SSR."
            ),
        }
        for col, key in enumerate(["PASS", "FAIL", "UNKNOWN", "Pass rate", "Unique SSRs"]):
            val_lbl = QLabel("—")
            val_lbl.setFont(QFont(FONT_MONO, 18, QFont.Weight.Bold))
            val_lbl.setObjectName("role_title")
            tip = metric_tooltips.get(key)
            key_lbl = _lbl(key, tip)
            key_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            val_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            sg.addWidget(val_lbl, 0, col)
            sg.addWidget(key_lbl, 1, col)
            self.metric_labels[key] = val_lbl

        # Pass/fail bar
        self.bar_label = QLabel("")
        self.bar_label.setFixedHeight(8)
        self.bar_label.setObjectName("role_muted_bg")
        sg.addWidget(self.bar_label, 2, 0, 1, 4)
        L.addWidget(self.summary_group)

        # Results table
        results_group = QGroupBox("Primer Results")
        rg = QVBoxLayout(results_group)

        # Filter checkbox
        filter_row = QHBoxLayout()
        self.pass_only_cb = QCheckBox("Show PASS primers only")
        self.pass_only_cb.setChecked(False)
        self.pass_only_cb.stateChanged.connect(self._populate_table)
        filter_row.addWidget(self.pass_only_cb)
        filter_row.addStretch()
        rg.addLayout(filter_row)

        self.table = QTableWidget()
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(False)
        self.table.setSortingEnabled(True)
        self.table.setFixedHeight(500)
        rg.addWidget(self.table)

        dl_row = QHBoxLayout()
        self.download_full_btn = QPushButton("Download full results (CSV)")
        self.download_full_btn.clicked.connect(self._download_full)
        self.download_pass_btn = QPushButton("Download PASS primers (CSV)")
        self.download_pass_btn.clicked.connect(self._download_pass)
        self.download_flanking_btn = QPushButton("Export Flanking Sequences (PASS, HTML)")
        self.download_flanking_btn.setToolTip(
            "For every PASS locus, exports the primer + SSR region plus ~10bp of "
            "flanking sequence on each side, as an HTML file with the primers and "
            "SSR region highlighted in colour for easy visual inspection in any browser."
        )
        self.download_flanking_btn.clicked.connect(self._download_flanking_html)
        self.download_flanking_fasta_btn = QPushButton("Export Flanking Sequences (PASS, FASTA)")
        self.download_flanking_fasta_btn.setToolTip(
            "Same region as the HTML export (primers + SSR + ~10bp flanking sequence), "
            "as plain-text FASTA for pasting into other tools. Primer and SSR positions "
            "are given as numbers in the header line, since plain text can't show colour."
        )
        self.download_flanking_fasta_btn.clicked.connect(self._download_flanking_fasta)
        dl_row.addWidget(self.download_full_btn)
        dl_row.addWidget(self.download_pass_btn)
        dl_row.addWidget(self.download_flanking_btn)
        dl_row.addWidget(self.download_flanking_fasta_btn)
        dl_row.addStretch()
        rg.addLayout(dl_row)
        L.addWidget(results_group)

        # Raw BLAST hits — written incrementally to disk during the run
        # itself (see core/primer_specificity_blast.py's hit_writer
        # mechanism), never held in memory or loaded into the app as a
        # table. For a large primer set against a repeat-dense genome,
        # every surviving hit's full alignment detail can be a very
        # large amount of data — this avoids ever building an in-memory
        # or in-app-table version of it at all, not just displaying it
        # efficiently once built.
        raw_group = QGroupBox("Raw BLAST Hits")
        raw_layout = QVBoxLayout(raw_group)
        raw_lbl = QLabel(
            "Every surviving BLAST hit (one row per primer-genome alignment) is written "
            "directly to a file during the run, not loaded into the app itself — for a "
            "large primer set this can be a substantial amount of data. Save a copy below."
        )
        raw_lbl.setObjectName("role_secondary")
        raw_lbl.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        raw_lbl.setWordWrap(True)
        raw_layout.addWidget(raw_lbl)

        self.raw_hits_status_lbl = QLabel("")
        self.raw_hits_status_lbl.setObjectName("role_secondary")
        self.raw_hits_status_lbl.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        raw_layout.addWidget(self.raw_hits_status_lbl)

        self.download_raw_btn = QPushButton("Save Raw BLAST Hits As\u2026")
        self.download_raw_btn.clicked.connect(self._download_raw)
        raw_layout.addWidget(self.download_raw_btn)
        L.addWidget(raw_group)

        L.addStretch()

        # No-results placeholder
        self.no_results_label = QLabel("No BLAST results yet — run Specificity Check first.")
        self.no_results_label.setObjectName("role_secondary")
        self.no_results_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        self.no_results_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        L.addWidget(self.no_results_label)

        self._refresh()

    # ---------------------------------------------------------
    # REFRESH
    # ---------------------------------------------------------
    def _refresh(self):
        has = self.state.has_specificity
        self.summary_group.setVisible(has)
        self.no_results_label.setVisible(not has)

        if not has:
            self._clear_tables()
            return

        # Skip re-render if data hasn't changed since last render
        if self.state.blast_version == self._last_rendered_version:
            return
        self._last_rendered_version = self.state.blast_version

        import pandas as pd
        spec_df   = pd.DataFrame(self.state.specificity_results)
        passed    = spec_df[spec_df["specificity_status"] == "PASS"]
        failed    = spec_df[spec_df["specificity_status"] == "FAIL"]
        unknown   = spec_df[spec_df["specificity_status"] == "UNKNOWN"]
        total     = max(len(spec_df), 1)
        pass_rate = 100 * len(passed) / total
        n_unique_pass_ssrs = passed["ssr_id"].nunique() if "ssr_id" in passed.columns else len(passed)

        self.metric_labels["PASS"].setText(str(len(passed)))
        self.metric_labels["PASS"].setObjectName("role_success")
        self.metric_labels["PASS"].setStyleSheet(f"font-family: {FONT_MONO}; font-size: 18pt; font-weight: bold;")
        self.metric_labels["FAIL"].setText(str(len(failed)))
        self.metric_labels["FAIL"].setObjectName("role_error")
        self.metric_labels["FAIL"].setStyleSheet(f"font-family: {FONT_MONO}; font-size: 18pt; font-weight: bold;")
        self.metric_labels["UNKNOWN"].setText(str(len(unknown)))
        self.metric_labels["UNKNOWN"].setObjectName("role_warning")
        self.metric_labels["UNKNOWN"].setStyleSheet(f"font-family: {FONT_MONO}; font-size: 18pt; font-weight: bold;")
        self.metric_labels["Pass rate"].setText(f"{pass_rate:.1f}%")
        self.metric_labels["Pass rate"].setObjectName("role_title")
        self.metric_labels["Pass rate"].setStyleSheet(f"font-family: {FONT_MONO}; font-size: 18pt; font-weight: bold;")
        self.metric_labels["Unique SSRs"].setText(str(n_unique_pass_ssrs))
        self.metric_labels["Unique SSRs"].setObjectName("role_title")
        self.metric_labels["Unique SSRs"].setStyleSheet(f"font-family: {FONT_MONO}; font-size: 18pt; font-weight: bold;")
        self.metric_labels["Unique SSRs"].setToolTip(
            f"{len(passed)} primer pairs PASS in total, across {n_unique_pass_ssrs} distinct SSR loci.\n"
            f"The PASS count can exceed this when a locus has more than one ranked primer pair "
            f"that each independently pass specificity — those are reported separately since panel "
            f"design can use whichever rank fits best, but they represent the same underlying locus."
        )

        # Colour bar
        bar_p = int(pass_rate); bar_f = int(100 * len(failed) / total)
        bar_u = max(0, 100 - bar_p - bar_f)
        self.bar_label.setStyleSheet(
            f"background: qlineargradient(x1:0, y1:0, x2:1, y2:0, "
            f"stop:0 {SUCCESS}, stop:{bar_p/100:.3f} {SUCCESS}, "
            f"stop:{bar_p/100:.3f} {ERROR}, stop:{(bar_p+bar_f)/100:.3f} {ERROR}, "
            f"stop:{(bar_p+bar_f)/100:.3f} {WARNING}, stop:1 {WARNING}); "
            f"border-radius: 3px;"
        )

        self._populate_table()
        self._update_raw_hits_status()

    def _populate_table(self):
        if not self.state.has_specificity:
            return
        import pandas as pd
        from PyQt6.QtGui import QColor
        # Strip nested list-of-dict fields (amplicons, hit records) before
        # building the DataFrame — pandas has real per-column overhead
        # materializing object-dtype columns holding arbitrary Python
        # structures, and the display table never uses these fields at
        # all (they're excluded a few lines below regardless). For a
        # large primer set this avoided cost is not trivial.
        _NESTED_FIELDS = {"amplicons", "left_hits", "right_hits", "offtarget_amplicons"}
        flat_records = [
            {k: v for k, v in rec.items() if k not in _NESTED_FIELDS}
            for rec in self.state.specificity_results
        ]
        # With imported pairs mixed in, many columns (repeat_count, motif
        # coordinates...) are missing for some rows; object dtype keeps each
        # cell's own type so those show blank rather than "nan", and designed
        # rows' integers aren't turned into floats.
        has_imported = any(r.get("source") == "imported" for r in flat_records)
        spec_df = _label_sources(pd.DataFrame(flat_records, dtype=object if has_imported else None))
        if self.pass_only_cb.isChecked():
            df = spec_df[spec_df["specificity_status"] == "PASS"].copy()
        else:
            df = spec_df.copy()

        COLS = {
            "ssr_id":             "SSR ID",
            "source":             "Source",
            "primer_name":        "Imported name",
            "pair_rank":          "Pair rank",
            "specificity_status": "Status",
            "contig":             "Contig",
            "motif":              "Motif",
            "repeat_count":       "Repeat count",
            "left_primer":        "Forward primer",
            "right_primer":       "Reverse primer",
            "product_size":       "Product size (bp)",
            "product_size_source": "Product size from",
            "expected_size":      "Expected size (bp)",
            "blast_amplicon_size": "BLAST amplicon size (bp)",
            "left_tm":            "Forward Tm (°C)",
            "right_tm":           "Reverse Tm (°C)",
            "left_gc":            "Forward GC (%)",
            "right_gc":           "Reverse GC (%)",
            "left_3end_dg":       "Forward 3' stability",
            "right_3end_dg":      "Reverse 3' stability",
        }
        exclude = {"amplicons", "left_hits", "right_hits", "offtarget_amplicons"}
        display_cols = [c for c in COLS if c in df.columns and c not in exclude]

        truncated = len(df) > TABLE_DISPLAY_LIMIT
        display_df = df.iloc[:TABLE_DISPLAY_LIMIT] if truncated else df

        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(display_df))
        self.table.setColumnCount(len(display_cols))
        self.table.setHorizontalHeaderLabels([COLS[c] for c in display_cols])

        for row_idx in range(len(display_df)):
            for col_idx, col in enumerate(display_cols):
                val = display_df.iat[row_idx, display_df.columns.get_loc(col)]
                
                # Use NumericTableWidgetItem for numeric columns
                if isinstance(val, float) and val != val:  # NaN — no value for this row
                    item = QTableWidgetItem("")
                elif isinstance(val, float):
                    text = f"{val:.2f}"
                    item = NumericTableWidgetItem(text, numeric_value=val)
                elif isinstance(val, int):
                    text = str(val)
                    item = NumericTableWidgetItem(text, numeric_value=val)
                elif val is None:
                    item = QTableWidgetItem("")
                else:
                    item = QTableWidgetItem(str(val))
                
                if col == "specificity_status":
                    if str(val) == "PASS":
                        item.setForeground(QColor(SUCCESS))
                    elif str(val) == "FAIL":
                        item.setForeground(QColor(ERROR))
                self.table.setItem(row_idx, col_idx, item)

        self.table.resizeColumnsToContents()
        self.table.setSortingEnabled(True)

        if truncated:
            # Update metrics label to note truncation
            current = self.metric_labels["Pass rate"].toolTip()
            self.table.setToolTip(
                f"Showing first {TABLE_DISPLAY_LIMIT:,} rows — download CSV for full results"
            )

    def _update_raw_hits_status(self):
        path = self.state.blast_raw_hits_path
        if not path or not os.path.exists(path):
            self.raw_hits_status_lbl.setText("No raw hits file available for this run.")
            self.download_raw_btn.setEnabled(False)
            return
        try:
            size_mb = os.path.getsize(path) / 1e6
            with open(path, encoding="utf-8") as f:
                n_rows = sum(1 for _ in f) - 1  # minus header
            self.raw_hits_status_lbl.setText(f"{n_rows:,} hit records \u2014 {size_mb:.1f} MB")
            self.download_raw_btn.setEnabled(True)
        except OSError:
            self.raw_hits_status_lbl.setText("Raw hits file could not be read.")
            self.download_raw_btn.setEnabled(False)

    def _clear_tables(self):
        self.table.setRowCount(0)

    # ---------------------------------------------------------
    # DOWNLOADS
    # ---------------------------------------------------------
    def _get_display_df(self):
        if not self.state.has_specificity:
            return None
        import pandas as pd
        records = self.state.specificity_results
        has_imported = any(r.get("source") == "imported" for r in records)
        # object dtype: see _populate_table — keeps designed rows' integers
        # as integers in the CSV when imported rows leave gaps.
        return _label_sources(pd.DataFrame(records, dtype=object if has_imported else None))

    def _clean_df(self, df):
        """Remove internal columns and rename for clean CSV output."""
        exclude = {"amplicons", "left_hits", "right_hits", "offtarget_amplicons",
                   "left_gc_clamp", "right_gc_clamp"}
        rename = {
            "ssr_id":             "SSR ID",
            "source":             "Source",
            "primer_name":        "Imported name",
            "expected_size":      "Expected size (bp)",
            "blast_amplicon_size": "BLAST amplicon size (bp)",
            "product_size_source": "Product size from",
            "pair_rank":          "Pair rank",
            "specificity_status": "Status",
            "contig":             "Contig",
            "start":              "Start (bp)",
            "end":                "End (bp)",
            "motif":              "Motif",
            "canonical_motif":    "Canonical motif",
            "repeat_count":       "Repeat count",
            "left_primer":        "Forward primer",
            "right_primer":       "Reverse primer",
            "product_size":       "Product size (bp)",
            "left_tm":            "Forward Tm (°C)",
            "right_tm":           "Reverse Tm (°C)",
            "left_gc":            "Forward GC (%)",
            "right_gc":           "Reverse GC (%)",
            "left_3end_dg":       "Forward 3' stability (kcal/mol)",
            "right_3end_dg":      "Reverse 3' stability (kcal/mol)",
        }
        cols = [col for col in df.columns if col not in exclude]
        return df[cols].rename(columns=rename)

    def _download_full(self):
        df = self._get_display_df()
        if df is None: return
        path, _ = QFileDialog.getSaveFileName(self, "Save full results", "specificity_results.csv", "CSV files (*.csv)")
        if path:
            self._clean_df(df).to_csv(path, index=False, encoding="utf-8-sig")
            self.mw.set_status(f"Saved to {path}")

    def _download_pass(self):
        df = self._get_display_df()
        if df is None: return
        path, _ = QFileDialog.getSaveFileName(self, "Save PASS primers", "pass_primers.csv", "CSV files (*.csv)")
        if path:
            self._clean_df(df[df["specificity_status"] == "PASS"]).to_csv(path, index=False, encoding="utf-8-sig")
            self.mw.set_status(f"Saved to {path}")

    def _download_flanking_html(self):
        if not self.state.has_specificity:
            self.mw.set_status("No BLAST specificity results available")
            return
        pass_records = [r for r in self.state.specificity_results if r.get("specificity_status") == "PASS"]
        if not pass_records:
            self.mw.set_status("No PASS primers available to export")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save flanking sequences", "flanking_sequences.html", "HTML files (*.html)"
        )
        if not path:
            return
        from core.primer_design import build_flanking_sequence_html
        html = build_flanking_sequence_html(self.state.genome, pass_records, flank_bp=10)
        with open(path, "w", encoding="utf-8") as f:
            f.write(html)
        self.mw.set_status(f"Saved to {path}")

    def _download_flanking_fasta(self):
        if not self.state.has_specificity:
            self.mw.set_status("No BLAST specificity results available")
            return
        pass_records = [r for r in self.state.specificity_results if r.get("specificity_status") == "PASS"]
        if not pass_records:
            self.mw.set_status("No PASS primers available to export")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save flanking sequences", "flanking_sequences.fasta", "FASTA files (*.fasta *.fa)"
        )
        if not path:
            return
        from core.primer_design import build_flanking_sequence_fasta
        fasta = build_flanking_sequence_fasta(self.state.genome, pass_records, flank_bp=10)
        with open(path, "w", encoding="utf-8") as f:
            f.write(fasta)
        self.mw.set_status(f"Saved to {path}")

    def _download_raw(self):
        src = self.state.blast_raw_hits_path
        if not src or not os.path.exists(src):
            self.mw.set_status("No raw hits file available for this run")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save raw BLAST hits", "blast_raw_hits.tsv", "Tab-separated files (*.tsv)"
        )
        if path:
            try:
                shutil.copy2(src, path)
                self.mw.set_status(f"Saved to {path}")
            except OSError as e:
                self.mw.set_status(f"Could not save raw hits file: {e}")

    def on_show(self):
        self._refresh()