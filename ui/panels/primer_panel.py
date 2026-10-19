"""
primer_panel.py — Primer Design
Single scroll area. Clean layout. Filters appear after design.
Supports Capillary Electrophoresis and Amplicon Sequencing modes.
"""
import time

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QGridLayout, QSpinBox, QDoubleSpinBox,
    QProgressBar, QTableWidget, QTableWidgetItem, QHeaderView,
    QAbstractItemView, QScrollArea, QFileDialog,
    QApplication,
)
from PyQt6.QtCore import QThread, pyqtSignal
from PyQt6.QtGui import QFont

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING, BG_MID, BG_LIGHT, BORDER,
)
from ui.table_utils import NumericTableWidgetItem, TABLE_DISPLAY_LIMIT, secondary_label as _lbl
from core.primer_design import PRESETS

PRIMER_DESIGN_WARN_THRESHOLD = 50_000


class PrimerWorker(QThread):
    progress = pyqtSignal(int, int)
    finished = pyqtSignal(dict, float)
    error    = pyqtSignal(str)

    def __init__(self, genome, ssr_list, params):
        super().__init__()
        self.genome   = genome
        self.ssr_list = ssr_list
        self.params   = params

    def run(self):
        try:
            from core.primer_design import design_primers_for_all_ssrs
            start = time.time()
            res = design_primers_for_all_ssrs(
                genome=self.genome,
                ssr_list=self.ssr_list,
                flank=self.params["flank"],
                product_size_range=(self.params["product_min"], self.params["product_max"]),
                preset=self.params["preset"],
                primer_opts=self.params["primer_opts"],
                num_pairs=self.params["num_pairs"],
                progress_callback=lambda d, t: self.progress.emit(d, t),
                amplicon_mode=self.params["amplicon_mode"],
            )
            self.finished.emit(res, time.time() - start)
        except Exception as e:
            self.error.emit(str(e))


class PrimerPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state = state
        self.mw    = main_window
        self._worker = None
        self._last_rendered_version = -1
        self._amplicon_mode = False
        self._filter_updating = False
        self._last_filter_key = None
        self._dg_empty_warned_key = None
        self._n_pass_dg = None  # pairs passing the 3' ΔG filter, before any fallback
        self._build_ui()

    # Maps every BASIC primer3 preset key onto its corresponding spinbox
    # attribute. Advanced settings (primer length, secondary structure,
    # thermodynamic parameters) are no longer exposed as widgets at all —
    # they're hardcoded to Primer3's own documented defaults directly in
    # core/primer_design.py's PRESETS, and applied via _run() reading the
    # active preset dict directly rather than any widget.
    _PRESET_KEY_TO_WIDGET = {
        "PRIMER_MIN_SIZE":           "primer_min_size",
        "PRIMER_OPT_SIZE":           "primer_opt_size",
        "PRIMER_MAX_SIZE":           "primer_max_size",
        "PRIMER_MIN_TM":             "tm_min",
        "PRIMER_OPT_TM":             "tm_opt",
        "PRIMER_MAX_TM":             "tm_max",
        "PRIMER_MIN_GC":             "gc_min",
        "PRIMER_MAX_GC":             "gc_max",
        "PRIMER_PAIR_MAX_DIFF_TM":   "max_tm_diff",
    }

    def _apply_preset_to_ui(self, preset_dict):
        """Push every BASIC value in a PRESETS entry onto its matching
        spinbox. Advanced values in the same preset dict are used
        directly from the dict itself at run time — see _run()."""
        for key, widget_name in self._PRESET_KEY_TO_WIDGET.items():
            if key in preset_dict:
                getattr(self, widget_name).setValue(preset_dict[key])

    def _on_mode_toggled(self, checked, mode):
        if not checked:
            return
        if mode == "capillary":
            self._btn_amplicon.setChecked(False)
            self._amplicon_mode = False
            self.state.workflow_mode = "capillary"
            self._mode_desc.setText(
                "Fragment analysis · Dye multiplexing\n"
                "Product sizes 100–500bp · Flank 50bp"
            )
            self.flank.setValue(50)
            self.product_min.setValue(100)
            self.product_max.setValue(500)
            self._apply_preset_to_ui(PRESETS["recommended"])
        else:  # amplicon
            self._btn_capillary.setChecked(False)
            self._amplicon_mode = True
            self.state.workflow_mode = "amplicon"
            self._mode_desc.setText(
                "Illumina amplicon sequencing\n"
                "Product sizes 80–200bp · Flank 50bp · Tighter Tm range"
            )
            self.flank.setValue(50)
            self.product_min.setValue(80)
            self.product_max.setValue(200)
            self._apply_preset_to_ui(PRESETS["amplicon"])
        QApplication.processEvents()

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

        title = QLabel("Primer Design")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setObjectName("role_title")
        sub = QLabel("Design PCR primers flanking each SSR using Primer3")
        sub.setObjectName("role_secondary")
        sub.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(title)
        L.addWidget(sub)

        # Mode selector
        mode_group = QGroupBox("Design Mode")
        mg = QVBoxLayout(mode_group)
        mg.setSpacing(8)

        mode_row = QHBoxLayout()
        mode_row.setSpacing(12)

        self._btn_capillary = QPushButton("Capillary Electrophoresis")
        self._btn_capillary.setCheckable(True)
        self._btn_capillary.setChecked(True)
        self._btn_capillary.setMinimumHeight(60)
        self._btn_capillary.setStyleSheet(f"""
            QPushButton {{
                background: {BG_MID};
                border: 1px solid {BORDER};
                border-radius: 8px;
                padding: 12px 20px;
                font-weight: bold;
                color: {TEXT_PRIMARY};
                text-align: left;
            }}
            QPushButton:checked {{
                background: qlineargradient(x1:0,y1:0,x2:1,y2:0,
                    stop:0 {ACCENT}33, stop:1 {ACCENT}11);
                border: 2px solid {ACCENT};
                color: {ACCENT};
            }}
            QPushButton:hover {{
                border: 1px solid {ACCENT}88;
                background: {BG_LIGHT};
            }}
        """)

        self._btn_amplicon = QPushButton("Amplicon Sequencing")
        self._btn_amplicon.setCheckable(True)
        self._btn_amplicon.setMinimumHeight(60)
        self._btn_amplicon.setStyleSheet(self._btn_capillary.styleSheet())

        self._mode_desc = QLabel(
            "Fragment analysis · M13 tails · Dye multiplexing\n"
            "Product sizes 100–500bp · Flank 50bp"
        )
        self._mode_desc.setWordWrap(True)
        self._mode_desc.setObjectName("role_secondary")
        self._mode_desc.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt; padding: 4px 0;")

        self._btn_capillary.toggled.connect(lambda checked: self._on_mode_toggled(checked, "capillary"))
        self._btn_amplicon.toggled.connect(lambda checked: self._on_mode_toggled(checked, "amplicon"))

        mode_row.addWidget(self._btn_capillary)
        mode_row.addWidget(self._btn_amplicon)
        mode_row.addStretch()
        mg.addLayout(mode_row)
        mg.addWidget(self._mode_desc)

        L.addWidget(mode_group)

        # Settings — no longer tabbed now that advanced settings are
        # hardcoded rather than user-configurable; a single-tab
        # QTabWidget would just be an unnecessary lone tab.
        basic_layout = QVBoxLayout()
        basic_layout.setContentsMargins(0, 0, 0, 0)

        basic_group = QGroupBox("PCR Parameters")
        bg = QGridLayout(basic_group); bg.setSpacing(10); bg.setContentsMargins(12, 12, 12, 12)
        bg.addWidget(_lbl(
            "Flanking region (bp)",
            "Minimum bases upstream/downstream of the SSR to search for primers.\n\n"
            "This is a floor, not a fixed value — if your Max product size "
            "requires more flanking sequence than this to be physically "
            "achievable, GRACE automatically searches a wider region. "
            "It will never search less than what you set here."
        ), 0, 0)
        self.flank = QSpinBox(); self.flank.setRange(50, 1000); self.flank.setValue(100); self.flank.setSingleStep(50)
        bg.addWidget(self.flank, 0, 1)
        bg.addWidget(_lbl(
            "Min product size (bp)",
            "Smallest acceptable amplicon size, including the SSR repeat "
            "itself plus both primer binding sites."
        ), 0, 2)
        self.product_min = QSpinBox(); self.product_min.setRange(50, 1000); self.product_min.setValue(100); self.product_min.setSingleStep(10)
        bg.addWidget(self.product_min, 0, 3)
        bg.addWidget(_lbl(
            "Max product size (bp)",
            "Largest acceptable amplicon size. GRACE automatically extends "
            "the flanking region searched per-locus if needed to make this "
            "size genuinely reachable, rather than silently capping it."
        ), 1, 0)
        self.product_max = QSpinBox(); self.product_max.setRange(50, 1000); self.product_max.setValue(250); self.product_max.setSingleStep(10)
        bg.addWidget(self.product_max, 1, 1)
        bg.addWidget(_lbl(
            "Primer pairs per SSR",
            "How many ranked candidate primer pairs to generate per locus, "
            "not just the single best-scoring one. Extra ranked pairs let "
            "specificity screening and panel design fall back to a "
            "different pair for the same locus if the top choice fails a "
            "later check or clashes with another marker."
        ), 1, 2)
        self.num_pairs = QSpinBox(); self.num_pairs.setRange(1, 5); self.num_pairs.setValue(2)
        bg.addWidget(self.num_pairs, 1, 3)
        bg.setColumnStretch(4, 1)
        basic_layout.addWidget(basic_group)

        primer_group = QGroupBox("Primer Length")
        pg = QGridLayout(primer_group); pg.setSpacing(10); pg.setContentsMargins(12, 12, 12, 12)
        pg.addWidget(_lbl(
            "Min (bp)",
            "Shortest primer Primer3 is allowed to design. Shorter primers "
            "bind less specifically."
        ), 0, 0)
        self.primer_min_size = QSpinBox(); self.primer_min_size.setRange(15, 35); self.primer_min_size.setValue(18)
        pg.addWidget(self.primer_min_size, 0, 1)
        pg.addWidget(_lbl(
            "Optimal (bp)",
            "Primer3's preferred length — it starts searching from here "
            "and only deviates toward Min/Max when needed to satisfy other "
            "constraints (Tm, GC content, secondary structure)."
        ), 0, 2)
        self.primer_opt_size = QSpinBox(); self.primer_opt_size.setRange(15, 35); self.primer_opt_size.setValue(20)
        pg.addWidget(self.primer_opt_size, 0, 3)
        pg.addWidget(_lbl(
            "Max (bp)",
            "Longest primer Primer3 is allowed to design. Longer primers "
            "cost more to synthesise and are more prone to secondary "
            "structure."
        ), 0, 4)
        self.primer_max_size = QSpinBox(); self.primer_max_size.setRange(15, 35); self.primer_max_size.setValue(23)
        pg.addWidget(self.primer_max_size, 0, 5)
        pg.setColumnStretch(6, 1)
        basic_layout.addWidget(primer_group)

        tm_group = QGroupBox("Melting Temperature (Tm)")
        tg = QGridLayout(tm_group); tg.setSpacing(10); tg.setContentsMargins(12, 12, 12, 12)
        tg.addWidget(_lbl(
            "Min (°C)",
            "Lowest acceptable primer melting temperature. Too low, and "
            "the primer may bind non-specifically at typical annealing "
            "temperatures."
        ), 0, 0)
        self.tm_min = QDoubleSpinBox(); self.tm_min.setRange(40, 75); self.tm_min.setValue(52); self.tm_min.setSingleStep(0.5)
        tg.addWidget(self.tm_min, 0, 1)
        tg.addWidget(_lbl(
            "Optimal (°C)",
            "Primer3's preferred melting temperature to design toward."
        ), 0, 2)
        self.tm_opt = QDoubleSpinBox(); self.tm_opt.setRange(40, 75); self.tm_opt.setValue(58); self.tm_opt.setSingleStep(0.5)
        tg.addWidget(self.tm_opt, 0, 3)
        tg.addWidget(_lbl(
            "Max (°C)",
            "Highest acceptable primer melting temperature."
        ), 0, 4)
        self.tm_max = QDoubleSpinBox(); self.tm_max.setRange(40, 75); self.tm_max.setValue(60); self.tm_max.setSingleStep(0.5)
        tg.addWidget(self.tm_max, 0, 5)
        tg.addWidget(_lbl("Max ΔTm (°C)", "Max difference between forward and reverse primer Tm.\nCapillary: 5.0°C, Amplicon: 1.0°C"), 1, 0)
        self.max_tm_diff = QDoubleSpinBox(); self.max_tm_diff.setRange(0, 100); self.max_tm_diff.setValue(5.0); self.max_tm_diff.setSingleStep(0.5)
        tg.addWidget(self.max_tm_diff, 1, 1)
        tg.setColumnStretch(6, 1)
        basic_layout.addWidget(tm_group)

        gc_group = QGroupBox("GC Content")
        gcg = QGridLayout(gc_group); gcg.setSpacing(10); gcg.setContentsMargins(12, 12, 12, 12)
        gcg.addWidget(_lbl(
            "Min GC %",
            "Lowest acceptable GC content. Too low reduces primer binding "
            "stability."
        ), 0, 0)
        self.gc_min = QDoubleSpinBox(); self.gc_min.setRange(0, 100); self.gc_min.setValue(40); self.gc_min.setSingleStep(5)
        gcg.addWidget(self.gc_min, 0, 1)
        gcg.addWidget(_lbl(
            "Max GC %",
            "Highest acceptable GC content. Too high increases risk of "
            "non-specific binding and secondary structure."
        ), 0, 2)
        self.gc_max = QDoubleSpinBox(); self.gc_max.setRange(0, 100); self.gc_max.setValue(60); self.gc_max.setSingleStep(5)
        gcg.addWidget(self.gc_max, 0, 3)
        gcg.setColumnStretch(4, 1)
        basic_layout.addWidget(gc_group)

        basic_layout.addStretch()
        L.addLayout(basic_layout)

        # Run controls
        run_group = QGroupBox("Run Primer Design")
        run_layout = QVBoxLayout(run_group)

        self.ssr_source_label = QLabel("No SSRs available")
        self.ssr_source_label.setObjectName("role_secondary")
        self.ssr_source_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        run_layout.addWidget(self.ssr_source_label)

        btn_row = QHBoxLayout()
        self.run_all_btn = QPushButton("Design for all SSRs")
        self.run_all_btn.setFixedHeight(40)
        self.run_all_btn.setEnabled(False)
        self.run_all_btn.clicked.connect(lambda: self._run("all"))
        self.run_sel_btn = QPushButton("Design for selected SSRs")
        self.run_sel_btn.setFixedHeight(40)
        self.run_sel_btn.setEnabled(False)
        self.run_sel_btn.clicked.connect(lambda: self._run("selected"))
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setFixedHeight(40)
        self.cancel_btn.setVisible(False)
        self.cancel_btn.clicked.connect(self._cancel)
        btn_row.addWidget(self.run_all_btn)
        btn_row.addWidget(self.run_sel_btn)
        btn_row.addWidget(self.cancel_btn)
        btn_row.addStretch()
        run_layout.addLayout(btn_row)

        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setFixedHeight(8)
        run_layout.addWidget(self.progress_bar)

        self.run_status = QLabel("")
        self.run_status.setObjectName("role_secondary")
        self.run_status.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        run_layout.addWidget(self.run_status)

        L.addWidget(run_group)

        # Design summary — kept separate from Results so it still shows
        # when a run designs zero pairs (Results is hidden in that case).
        self.summary_group = QGroupBox("Design Summary")
        self.summary_group.setVisible(False)
        sum_layout = QVBoxLayout(self.summary_group)

        cards_grid = QGridLayout(); cards_grid.setSpacing(8)
        self._cards = {}
        for col, (key, title) in enumerate([
            ("attempted", "SSRs attempted"),
            ("with",      "SSRs with primers"),
            ("pairs",     "Primer pairs"),
            ("without",   "SSRs without primers"),
            ("dg",        "Pass 3' ΔG filter"),
        ]):
            card = QWidget(); card.setObjectName("stat_card")
            cl = QVBoxLayout(card); cl.setSpacing(2)
            value = QLabel("—")
            value.setFont(QFont(FONT_MONO, FONT_SIZE_LARGE, QFont.Weight.Bold))
            value.setObjectName("role_title")
            caption = QLabel("")
            caption.setObjectName("role_secondary")
            caption.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
            cl.addWidget(value)
            cl.addWidget(QLabel(title))
            cl.addWidget(caption)
            cards_grid.addWidget(card, 0, col)
            self._cards[key] = (value, caption)
        sum_layout.addLayout(cards_grid)

        self.issues_label = QLabel("")
        self.issues_label.setWordWrap(True)
        self.issues_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        sum_layout.addWidget(self.issues_label)

        issues_btn_row = QHBoxLayout()
        self.issues_toggle_btn = QPushButton("Show affected loci")
        self.issues_toggle_btn.setCheckable(True)
        self.issues_toggle_btn.toggled.connect(self._on_issues_toggled)
        self.issues_export_btn = QPushButton("Export affected loci (CSV)")
        self.issues_export_btn.clicked.connect(self._download_issues_csv)
        issues_btn_row.addWidget(self.issues_toggle_btn)
        issues_btn_row.addWidget(self.issues_export_btn)
        issues_btn_row.addStretch()
        sum_layout.addLayout(issues_btn_row)

        self.issues_table = QTableWidget()
        self.issues_table.setAlternatingRowColors(True)
        self.issues_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.issues_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.issues_table.setSortingEnabled(True)
        self.issues_table.setFixedHeight(260)
        self.issues_table.setVisible(False)
        sum_layout.addWidget(self.issues_table)

        L.addWidget(self.summary_group)

        # Results
        self.results_group = QGroupBox("Results")
        self.results_group.setVisible(False)
        res_layout = QVBoxLayout(self.results_group)

        # Only used for the "showing first N rows" note — the headline
        # counts live in the Design Summary cards above.
        self.metrics_label = QLabel("")
        self.metrics_label.setObjectName("role_secondary")
        self.metrics_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        self.metrics_label.setVisible(False)
        res_layout.addWidget(self.metrics_label)

        filter_group = QGroupBox("Quality Filters")
        fg = QGridLayout(filter_group); fg.setSpacing(8); fg.setContentsMargins(10, 10, 10, 10)
        
        # Tooltip explaining 3' end stability
        dg_tooltip = (
            "3' end stability (ΔG in kcal/mol)\n\n"
            "More negative = stronger binding at primer's 3' end\n"
            "Less negative = weaker binding\n\n"
            "Recommended range: -9 to -5 kcal/mol\n"
            "• Too stable (< -9): Risk of non-specific priming\n"
            "• Too weak (> -5): Poor extension by polymerase\n\n"
            "The 3' end is critical because DNA polymerase extends from here.\n"
            "If binding is too weak, the primer dissociates before extension.\n"
            "If too strong, the primer may bind to partially matched sequences."
        )
        
        fg.addWidget(_lbl("Min 3' ΔG (kcal/mol)", dg_tooltip), 0, 0)
        self.dg_min = QDoubleSpinBox(); self.dg_min.setRange(-20, 0); self.dg_min.setValue(-9); self.dg_min.setSingleStep(0.5)
        self.dg_min.setToolTip(dg_tooltip)
        # Use editingFinished instead of valueChanged to prevent constant re-filtering
        self.dg_min.editingFinished.connect(self._on_filter_changed)
        fg.addWidget(self.dg_min, 0, 1)
        fg.addWidget(_lbl("Max 3' ΔG (kcal/mol)", dg_tooltip), 0, 2)
        self.dg_max = QDoubleSpinBox(); self.dg_max.setRange(-20, 0); self.dg_max.setValue(-5); self.dg_max.setSingleStep(0.5)
        self.dg_max.setToolTip(dg_tooltip)
        # Use editingFinished instead of valueChanged to prevent constant re-filtering
        self.dg_max.editingFinished.connect(self._on_filter_changed)
        fg.addWidget(self.dg_max, 0, 3)
        self.filter_status = QLabel("")
        self.filter_status.setObjectName("role_secondary")
        self.filter_status.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        fg.addWidget(self.filter_status, 1, 0, 1, 4)
        fg.setColumnStretch(4, 1)
        res_layout.addWidget(filter_group)

        self.table = QTableWidget()
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSortingEnabled(True)
        self.table.setFixedHeight(500)
        res_layout.addWidget(self.table)

        dl_row = QHBoxLayout()
        self.download_csv_btn = QPushButton("Download primer table (CSV)")
        self.download_csv_btn.clicked.connect(self._download_csv)
        self.download_fasta_btn = QPushButton("Download BLAST-ready FASTA")
        self.download_fasta_btn.clicked.connect(self._download_fasta)
        dl_row.addWidget(self.download_csv_btn)
        dl_row.addWidget(self.download_fasta_btn)
        dl_row.addStretch()
        res_layout.addLayout(dl_row)

        L.addWidget(self.results_group)
        L.addStretch()
        
        # Initialize Capillary mode defaults (button is already checked but callback wasn't triggered)
        self._on_mode_toggled(True, "capillary")
        
        self._refresh()

    def _run(self, mode):
        if not self.state.has_genome:
            self._set_status("Load a genome first", "role_warning"); return
        if not self.state.has_ssrs:
            self._set_status("Run SSR detection first", "role_warning"); return
        if mode == "selected":
            ssr_list = self.state.selected_ssrs
            if not ssr_list:
                self._set_status("No SSRs selected", "role_warning"); return
        else:
            ssr_list = self.state.ssrs

        if len(ssr_list) > PRIMER_DESIGN_WARN_THRESHOLD:
            from PyQt6.QtWidgets import QMessageBox
            msg = QMessageBox(self)
            msg.setWindowTitle("Large Dataset Warning")
            msg.setText(
                f"You are about to design primers for {len(ssr_list):,} SSRs.\n\n"
                f"This may take a long time.\n"
                f"Consider selecting a subset on the SSR Detection page first."
            )
            msg.addButton("Continue anyway", QMessageBox.ButtonRole.AcceptRole)
            cancel = msg.addButton("Cancel", QMessageBox.ButtonRole.RejectRole)
            msg.exec()
            if msg.clickedButton() == cancel:
                return

        # Start from the active preset — it already has every field,
        # including every now-hardcoded advanced setting — and override
        # only the basic fields that still have a widget on screen. This
        # keeps PRESETS in core/primer_design.py as the single source of
        # truth for advanced settings, exactly as before.
        active_preset = PRESETS["amplicon"] if self._amplicon_mode else PRESETS["recommended"]
        primer_opts = dict(active_preset)
        primer_opts.update({
            "PRIMER_MIN_SIZE":         self.primer_min_size.value(),
            "PRIMER_OPT_SIZE":         self.primer_opt_size.value(),
            "PRIMER_MAX_SIZE":         self.primer_max_size.value(),
            "PRIMER_MIN_TM":           self.tm_min.value(),
            "PRIMER_OPT_TM":           self.tm_opt.value(),
            "PRIMER_MAX_TM":           self.tm_max.value(),
            "PRIMER_MIN_GC":           self.gc_min.value(),
            "PRIMER_MAX_GC":           self.gc_max.value(),
            "PRIMER_PAIR_MAX_DIFF_TM": self.max_tm_diff.value(),
        })
        params = {
            "flank":        self.flank.value(),
            "product_min":  self.product_min.value(),
            "product_max":  self.product_max.value(),
            "num_pairs":    self.num_pairs.value(),
            "primer_opts":  primer_opts,
            "preset":       "amplicon" if self._amplicon_mode else "recommended",
            "amplicon_mode": self._amplicon_mode,
        }
        self.state.product_min = params["product_min"]
        self.state.product_max = params["product_max"]
        self.run_all_btn.setEnabled(False); self.run_sel_btn.setEnabled(False)
        self.cancel_btn.setVisible(True)
        self.progress_bar.setVisible(True); self.progress_bar.setValue(0)
        mode_tag = " [Amplicon]" if self._amplicon_mode else " [Capillary]"
        self._set_status(f"Running Primer3...{mode_tag}", "role_secondary")
        self.mw.set_status("Designing primers...")

        self._worker = PrimerWorker(self.state.genome, ssr_list, params)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_progress(self, done, total):
        self.progress_bar.setMaximum(total)
        self.progress_bar.setValue(done)
        if done % 100 == 0 or done == total:
            QApplication.processEvents()

    # fail_reason from core.primer_design → (Issue column text, fallback detail)
    _FAIL_REASONS = {
        "no_pairs":        ("No primer pair found", "Primer3 found no pair meeting the current settings"),
        "repetitive_only": ("Only repetitive primers",
                            "Primer3 returned pairs, but every one was rejected as a short-repeat primer"),
        "design_error":    ("Primer3 error", ""),
    }

    @staticmethod
    def _issue_row(ssr, severity, issue, detail):
        return {
            "ssr_id": ssr.get("ssr_id"), "contig": ssr.get("contig"),
            "start": ssr.get("start"), "end": ssr.get("end"),
            "motif": ssr.get("imperfect_structure") or ssr.get("motif"),
            "severity": severity, "issue": issue, "detail": detail,
        }

    def _on_done(self, result, elapsed):
        success      = result["success"]
        failed       = result.get("failed", [])
        size_limited = result.get("size_limited", [])
        params       = self._worker.params

        issues = []
        for ssr in failed:
            issue, fallback = self._FAIL_REASONS.get(
                ssr.get("fail_reason"), ("No primer pair found", ""))
            detail = ssr.get("fail_detail") or ssr.get("design_error") or fallback
            issues.append(self._issue_row(ssr, "failed", issue, detail))
        for ssr in size_limited:
            issues.append(self._issue_row(
                ssr, "note", "Near contig edge",
                "Too close to a contig edge to reach your max product size — "
                "primers were still designed within the available sequence"))

        self.state.primer_results = success
        self.state.filtered_primer_results = None
        self.state.primer_design_summary = {
            "n_attempted": len(self._worker.ssr_list),
            "num_pairs":   params["num_pairs"],
            "mode":        "Amplicon" if params["amplicon_mode"] else "Capillary",
            "elapsed":     elapsed,
            "issues":      issues,
        }
        self.state.primer_version += 1
        self.run_all_btn.setEnabled(True); self.run_sel_btn.setEnabled(True)
        self.cancel_btn.setVisible(False)
        self.progress_bar.setVisible(False)
        self.issues_toggle_btn.setChecked(False)
        n = len(success)
        mode = self.state.primer_design_summary["mode"]
        if n == 0:
            self._set_status(
                f"No primers could be designed ({elapsed:.1f}s, {mode} mode) — "
                f"see the Design Summary below for why", "role_error")
        else:
            self._set_status(f"Completed in {elapsed:.1f}s · {mode} mode", "role_success")
        self.mw.set_status(f"Primer design complete — {n:,} pairs designed")
        self.mw.on_step_complete(3)
        self._cleanup_worker()
        self._refresh()

    def _on_error(self, msg):
        self.run_all_btn.setEnabled(True); self.run_sel_btn.setEnabled(True)
        self.cancel_btn.setVisible(False)
        self.progress_bar.setVisible(False)
        self._set_status(f"Error: {msg}", "role_error")
        self.mw.set_status("Primer design failed")
        self._cleanup_worker()

    def _cancel(self):
        if self._worker:
            self._worker.terminate()
            self._worker.wait()
            self._cleanup_worker()
        self.run_all_btn.setEnabled(True); self.run_sel_btn.setEnabled(True)
        self.cancel_btn.setVisible(False)
        self.progress_bar.setVisible(False)
        self._set_status("Cancelled", "role_warning")
        self.mw.set_status("Primer design cancelled")

    def _cleanup_worker(self):
        if self._worker:
            try:
                self._worker.progress.disconnect()
                self._worker.finished.disconnect()
                self._worker.error.disconnect()
            except Exception:
                pass
            # See ssr_panel.py's _cleanup_worker() for why this matters
            # even though the thread has already finished by this point.
            self._worker.wait()
            self._worker = None

    def _set_status(self, msg, role):
        self.run_status.setText(msg)
        self.run_status.setObjectName(role)
        self.run_status.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")

    def _on_filter_changed(self):
        if self._filter_updating: return
        # editingFinished also fires on focus loss — skip if nothing changed.
        key = (self.state.primer_version, self.dg_min.value(), self.dg_max.value())
        if key == self._last_filter_key: return
        self._last_filter_key = key
        self._filter_updating = True
        try:
            self._apply_filters()
            self._populate_table()
        finally:
            self._filter_updating = False

    def _apply_filters(self):
        """Apply quality filters to primer results."""
        if not self.state.primer_results:
            return
        import pandas as pd

        try:
            dg_min = float(self.dg_min.value())
            dg_max = float(self.dg_max.value())
        except Exception:
            return

        if dg_min > dg_max:
            dg_min, dg_max = dg_max, dg_min

        df = pd.DataFrame(self.state.primer_results)
        total = len(df)

        # 3' stability filter - with error handling for missing or invalid values
        if "left_3end_dg" in df.columns and "right_3end_dg" in df.columns:
            try:
                # Convert to numeric, coercing errors to NaN
                df["left_3end_dg"] = pd.to_numeric(df["left_3end_dg"], errors='coerce')
                df["right_3end_dg"] = pd.to_numeric(df["right_3end_dg"], errors='coerce')
                
                # Only filter rows where both values are valid numbers
                mask = (
                    df["left_3end_dg"].notna() &
                    df["right_3end_dg"].notna() &
                    df["left_3end_dg"].between(dg_min, dg_max) &
                    df["right_3end_dg"].between(dg_min, dg_max)
                )
                df = df[mask].copy()
            except Exception as e:
                # If filtering fails, log the error and use all primers
                print(f"Warning: 3' stability filter failed: {e}")
                df = pd.DataFrame(self.state.primer_results)

        self._n_pass_dg = len(df)

        # Update filtered results
        if df.empty:
            self.state.filtered_primer_results = self.state.primer_results
            self.filter_status.setText("No primers pass current filters — using full set for BLAST")
            self.filter_status.setObjectName("role_warning")
            self.filter_status.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
            # Warn once per (primer run, ΔG range) — editingFinished also
            # fires on focus loss, and _refresh re-applies filters.
            warn_key = (self.state.primer_version, dg_min, dg_max)
            if warn_key != self._dg_empty_warned_key:
                self._dg_empty_warned_key = warn_key
                from PyQt6.QtWidgets import QMessageBox
                QMessageBox.warning(
                    self, "No Primers Pass 3' ΔG Filter",
                    f"None of the {total:,} primer pairs have both 3' ΔG values "
                    f"within {dg_min:g} to {dg_max:g} kcal/mol.\n\n"
                    f"The filter is being ignored: the full, unfiltered set of "
                    f"{total:,} pairs will be used for BLAST and downstream steps.\n\n"
                    f"Widen the 3' ΔG range if you want the filter applied."
                )
        else:
            self.state.filtered_primer_results = df.to_dict(orient="records")
            removed = total - len(df)
            if removed > 0:
                self.filter_status.setText(f"{removed:,} pairs removed — {len(df):,} remaining")
                self.filter_status.setObjectName("role_warning")
                self.filter_status.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
            else:
                self.filter_status.setText(f"All {len(df):,} pairs pass filters")
                self.filter_status.setObjectName("role_success")
                self.filter_status.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")

    def _refresh(self):
        if not self.state.has_genome or not self.state.has_ssrs:
            self.run_all_btn.setEnabled(False); self.run_sel_btn.setEnabled(False)
        else:
            self.run_all_btn.setEnabled(True)
            n_sel = len(self.state.selected_ssrs) if self.state.selected_ssrs else 0
            self.ssr_source_label.setText(f"{len(self.state.ssrs):,} SSRs available" + (f" — {n_sel:,} selected" if n_sel else ""))
        has_summary = self.state.primer_design_summary is not None
        if not self.state.has_primers and not has_summary:
            self.summary_group.setVisible(False)
            self.results_group.setVisible(False); return
        if self.state.primer_version == self._last_rendered_version: return
        self._last_rendered_version = self.state.primer_version
        self.summary_group.setVisible(True)
        self._populate_issues()
        if self.state.has_primers:
            self.results_group.setVisible(True)
            self._apply_filters()
            self._populate_table()   # also refreshes the summary cards
        else:
            self.results_group.setVisible(False)
            self._n_pass_dg = None
            self._update_summary_cards()

    def _populate_table(self):
        if not self.state.has_primers: return
        import pandas as pd
        primers = self.state.filtered_primer_results or self.state.primer_results
        df = pd.DataFrame(primers)
        truncated = len(df) > TABLE_DISPLAY_LIMIT
        display_df = df.iloc[:TABLE_DISPLAY_LIMIT] if truncated else df
        self._update_summary_cards()
        self.metrics_label.setVisible(truncated)
        if truncated:
            self.metrics_label.setText(
                f"Showing the first {TABLE_DISPLAY_LIMIT:,} of {len(df):,} pairs — "
                f"download the CSV for the full table")

        COLS = {
            "ssr_id": "SSR ID", "pair_rank": "Pair", "contig": "Contig",
            "motif": "Motif", "left_primer": "Forward", "right_primer": "Reverse",
            "product_size": "Size", "left_tm": "Fwd Tm", "right_tm": "Rev Tm",
            "left_gc": "Fwd GC", "right_gc": "Rev GC",
            "left_3end_dg": "Fwd 3' ΔG", "right_3end_dg": "Rev 3' ΔG",
        }

        # Add Chromosome column if mapping exists
        show_chromosome = hasattr(self.state, 'get_display_name') and self.state.chrom_names
        if show_chromosome:
            COLS["chromosome"] = "Chromosome"

        display_cols = [c for c in COLS if c in display_df.columns]

        # Insert Chromosome after Contig
        if show_chromosome and "chromosome" in display_cols and "contig" in display_cols:
            display_cols.remove("chromosome")
            contig_idx = display_cols.index("contig")
            display_cols.insert(contig_idx + 1, "chromosome")

        header = self.table.horizontalHeader()
        header.setStretchLastSection(True)

        # Sections left in ResizeToContents by a previous populate re-measure
        # every row on each setItem — O(n²), which froze the UI when the 3' ΔG
        # filter repopulated ~2,000 pairs. Drop to Interactive while filling
        # and size once at the end.
        header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        self.table.setSortingEnabled(False)
        self.table.setUpdatesEnabled(False)
        try:
            self.table.clearContents()
            self.table.setRowCount(len(display_df))
            self.table.setColumnCount(len(display_cols))
            self.table.setHorizontalHeaderLabels([COLS.get(c, c) for c in display_cols])

            columns = {c: display_df[c].tolist() for c in display_cols if c != "chromosome"}
            if "chromosome" in display_cols:
                columns["chromosome"] = [self.state.get_display_name(c)
                                         for c in display_df["contig"].tolist()]

            for col_idx, col in enumerate(display_cols):
                for row_idx, val in enumerate(columns[col]):
                    # Format text for display
                    if isinstance(val, float):
                        item = NumericTableWidgetItem(f"{val:.2f}", numeric_value=val)
                    elif isinstance(val, int):
                        item = NumericTableWidgetItem(str(val), numeric_value=val)
                    elif val is None:
                        item = QTableWidgetItem("")
                    else:
                        item = QTableWidgetItem(str(val))
                    self.table.setItem(row_idx, col_idx, item)

            # Column sizing runs AFTER population — sizing to contents needs
            # the actual cell text to measure against, otherwise longer values
            # get truncated. One-shot resize, not a ResizeToContents mode.
            for ci, col in enumerate(display_cols):
                if col == "contig":
                    self.table.setColumnWidth(ci, 140)
                elif col == "chromosome":
                    self.table.setColumnWidth(ci, 120)
                elif col in ("left_primer", "right_primer"):
                    self.table.setColumnWidth(ci, 200)
                else:
                    self.table.resizeColumnToContents(ci)
        finally:
            self.table.setUpdatesEnabled(True)
            self.table.setSortingEnabled(True)

    @staticmethod
    def _set_role(widget, role):
        widget.setObjectName(role)
        widget.style().unpolish(widget)
        widget.style().polish(widget)

    def _update_summary_cards(self):
        summary = self.state.primer_design_summary or {}
        primers = self.state.primer_results or []
        n_pairs = len(primers)
        n_with  = len({p.get("ssr_id") for p in primers})
        n_att   = summary.get("n_attempted")
        n_without = sum(1 for i in summary.get("issues", []) if i["severity"] == "failed")

        def card(key, value, caption="", role="role_title"):
            value_lbl, caption_lbl = self._cards[key]
            value_lbl.setText(value)
            caption_lbl.setText(caption)
            self._set_role(value_lbl, role)

        def pct(n):
            return f"{100 * n / n_att:.1f}% of attempted" if n_att else ""

        # Sessions saved before the summary existed only have the pairs.
        if n_att is None:
            card("attempted", "—", "not recorded for this session")
            card("without",   "—", "not recorded for this session")
            card("with",      f"{n_with:,}")
            card("pairs",     f"{n_pairs:,}")
        else:
            card("attempted", f"{n_att:,}", f"{summary['mode']} mode")
            card("without",   f"{n_without:,}", pct(n_without),
                 "role_warning" if n_without else "role_title")
            card("with",      f"{n_with:,}", pct(n_with),
                 "role_title" if n_with else "role_error")
            card("pairs",     f"{n_pairs:,}", f"up to {summary['num_pairs']} per SSR")

        if self._n_pass_dg is None or not n_pairs:
            card("dg", "—")
        else:
            card("dg", f"{self._n_pass_dg:,}", f"of {n_pairs:,} pairs",
                 "role_warning" if self._n_pass_dg == 0 else "role_title")

    def _populate_issues(self):
        summary = self.state.primer_design_summary
        issues = (summary or {}).get("issues", [])
        failed = [i for i in issues if i["severity"] == "failed"]
        notes  = [i for i in issues if i["severity"] == "note"]

        has_issues = bool(issues)
        self.issues_toggle_btn.setVisible(has_issues)
        self.issues_export_btn.setVisible(has_issues)
        if not has_issues:
            self.issues_table.setVisible(False)

        if summary is None:
            self.issues_label.setVisible(False)
            return
        self.issues_label.setVisible(True)

        if not has_issues:
            self.issues_label.setText("Every attempted SSR has at least one primer pair.")
            self._set_role(self.issues_label, "role_success")
            return

        # One line per category, e.g. "12 SSRs without primers: 10 no primer
        # pair found · 2 only repetitive primers"
        from collections import Counter
        def ssrs(n): return f"{n:,} SSR" + ("" if n == 1 else "s")
        lines = []
        if failed:
            by_issue = Counter(i["issue"] for i in failed)
            parts = " · ".join(f"{n:,} {issue.lower()}" for issue, n in by_issue.most_common())
            lines.append(f"<b>{ssrs(len(failed))} without primers:</b> {parts}")
        if notes:
            lines.append(
                f"<b>{ssrs(len(notes))} near a contig edge:</b> primers designed, "
                f"but your max product size wasn't reachable")
        self.issues_label.setText("<br>".join(lines))
        self._set_role(self.issues_label, "role_warning")

        self._on_issues_toggled(self.issues_toggle_btn.isChecked())

        cols = [("ssr_id", "SSR ID"), ("contig", "Contig"), ("start", "Start"),
                ("end", "End"), ("motif", "Motif"), ("issue", "Issue"), ("detail", "Details")]
        rows = issues[:TABLE_DISPLAY_LIMIT]
        t = self.issues_table
        t.setSortingEnabled(False)
        t.clearContents()
        t.setRowCount(len(rows))
        t.setColumnCount(len(cols))
        t.setHorizontalHeaderLabels([label for _, label in cols])
        for r, row in enumerate(rows):
            for c, (key, _) in enumerate(cols):
                val = row.get(key)
                if key == "contig" and val is not None:
                    val = self.state.get_display_name(val)
                if isinstance(val, int):
                    item = NumericTableWidgetItem(str(val), numeric_value=val)
                else:
                    item = QTableWidgetItem("" if val is None else str(val))
                    if key == "detail":
                        item.setToolTip(item.text())
                t.setItem(r, c, item)
        # Shrink to fit short lists; cap at 260px and scroll beyond that.
        row_h = t.verticalHeader().defaultSectionSize()
        t.setFixedHeight(min(260, t.horizontalHeader().sizeHint().height() + row_h * len(rows) + 4))
        header = t.horizontalHeader()
        for c in range(len(cols) - 1):
            header.setSectionResizeMode(c, QHeaderView.ResizeMode.ResizeToContents)
        header.setStretchLastSection(True)
        t.setSortingEnabled(True)

    def _on_issues_toggled(self, checked):
        self.issues_table.setVisible(checked)
        n = len((self.state.primer_design_summary or {}).get("issues", []))
        self.issues_toggle_btn.setText("Hide affected loci" if checked else f"Show affected loci ({n:,})")

    def _download_issues_csv(self):
        issues = (self.state.primer_design_summary or {}).get("issues", [])
        if not issues: return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save affected loci", "primer_design_issues.csv", "CSV files (*.csv)")
        if not path: return
        import pandas as pd
        df = pd.DataFrame(issues)
        df.insert(2, "Chromosome", df["contig"].apply(self.state.get_display_name))
        df["severity"] = df["severity"].map({"failed": "No primers", "note": "Primers designed"})
        df.rename(columns={
            "ssr_id": "SSR ID", "contig": "Contig", "start": "Start", "end": "End",
            "motif": "Motif", "severity": "Outcome", "issue": "Issue", "detail": "Details",
        }).to_csv(path, index=False, encoding="utf-8-sig")
        self.mw.set_status(f"Saved to {path}")

    def _download_csv(self):
        if not self.state.has_primers: return
        path, _ = QFileDialog.getSaveFileName(self, "Save primer table", "primer_table.csv", "CSV files (*.csv)")
        if not path: return
        import pandas as pd
        primers = self.state.filtered_primer_results or self.state.primer_results
        df = pd.DataFrame(primers)
        if hasattr(self.state, 'get_display_name'):
            df["Chromosome"] = df["contig"].apply(lambda c: self.state.get_display_name(c))
        df.rename(columns={
            "ssr_id": "SSR ID", "pair_rank": "Pair", "contig": "Contig",
            "motif": "Motif", "repeat_count": "Repeat count",
            "left_primer": "Forward", "right_primer": "Reverse",
            "product_size": "Size", "left_tm": "Fwd Tm", "right_tm": "Rev Tm",
            "left_gc": "Fwd GC", "right_gc": "Rev GC",
            "left_3end_dg": "Fwd 3' ΔG", "right_3end_dg": "Rev 3' ΔG",
        }).to_csv(path, index=False, encoding="utf-8-sig")
        self.mw.set_status(f"Saved to {path}")

    def _download_fasta(self):
        if not self.state.has_primers: return
        path, _ = QFileDialog.getSaveFileName(self, "Save BLAST FASTA", "primers_blast_ready.fasta", "FASTA files (*.fasta *.fa)")
        if not path: return
        from core.primer_design import primers_to_blast_fasta
        primers = self.state.filtered_primer_results or self.state.primer_results
        with open(path, "w") as f: f.write(primers_to_blast_fasta(primers))
        self.mw.set_status(f"Saved to {path}")

    def on_show(self):
        self._refresh()