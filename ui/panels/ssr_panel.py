"""
ssr_panel.py — SSR Detection
Single scroll area. Table has fixed height. No layout tricks.
Added adjacent same‑motif filter to reduce problematic loci.
"""
import time

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QGridLayout, QSpinBox, QCheckBox, QComboBox,
    QProgressBar, QTableWidget, QTableWidgetItem, QHeaderView,
    QAbstractItemView, QScrollArea, QFileDialog,
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QFont

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY,
    FONT_UI, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING, BG_MID, BORDER,
)
from ui.table_utils import NumericTableWidgetItem, TABLE_DISPLAY_LIMIT, secondary_label as _lbl


# ---------------------------------------------------------------------------
# SSR detection worker
# ---------------------------------------------------------------------------

class SSRWorker(QThread):
    progress = pyqtSignal(int, int)
    finished = pyqtSignal(list, float)
    error    = pyqtSignal(str)

    def __init__(self, genome, params):
        super().__init__()
        self.genome = genome
        self.params = params

    def run(self):
        try:
            from core.ssr_detection import find_ssrs
            start = time.time()
            ssrs = find_ssrs(
                self.genome,
                motif_lengths=self.params["motif_lengths"],
                min_repeats=self.params["min_repeats"],
                search_reverse=self.params["search_reverse"],
                motif_standardisation_level=self.params["std_level"],
                min_contig_len=self.params["min_contig_len"],
                n_workers=self.params.get("n_workers"),
                progress_callback=lambda d, t: self.progress.emit(d, t),
            )
            self.finished.emit(ssrs, time.time() - start)
        except Exception as e:
            self.error.emit(str(e))


# ---------------------------------------------------------------------------
# GFF annotation worker
# ---------------------------------------------------------------------------

class GFFAnnotationWorker(QThread):
    # annotated ssrs, gff_index, chrom_names dict
    finished = pyqtSignal(list, object, dict)
    error    = pyqtSignal(str)

    def __init__(self, ssrs, gff_path, gff_features):
        super().__init__()
        self.ssrs         = ssrs
        self.gff_path     = gff_path
        self.gff_features = gff_features

    def run(self):
        try:
            from core.gff_parser import build_gff_index, annotate_ssrs
            if self.gff_features is None:
                index = build_gff_index(self.gff_path)
            else:
                index = self.gff_features
            annotated = annotate_ssrs(self.ssrs, index)
            self.finished.emit(annotated, index, dict(index.chrom_names))
        except Exception as e:
            self.error.emit(str(e))


# ---------------------------------------------------------------------------
# Panel
# ---------------------------------------------------------------------------

class SSRPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state   = state
        self.mw      = main_window
        self._worker     = None
        self._gff_worker = None
        self._last_rendered_version = -1
        self._n_workers_auto_value = None   # last value WE auto-set, to detect manual overrides
        self._n_workers_genome_id  = None   # id() of the genome dict this recommendation was computed for
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
        layout = QVBoxLayout(content)
        layout.setContentsMargins(PANEL_PADDING, PANEL_PADDING, PANEL_PADDING, PANEL_PADDING)
        layout.setSpacing(16)

        # Title
        title = QLabel("SSR Detection")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setObjectName("role_title")
        sub = QLabel("Scan your genome for simple sequence repeats")
        sub.setObjectName("role_secondary")
        sub.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        layout.addWidget(title)
        layout.addWidget(sub)

        # ============================================================
        # GROUP 1 — What to search for
        # ============================================================
        search_group = QGroupBox("What to Search For")
        sg1 = QGridLayout(search_group)
        sg1.setSpacing(10)
        sg1.setColumnStretch(3, 1)

        sg1.addWidget(_lbl(
            "Motif length & minimum repeat count",
            "For each repeat unit size you want to search for, tick the box\n"
            "and set how many times it must repeat to count as an SSR.\n"
            "E.g. trinucleotide with min repeats 5 finds motifs of 3 bases\n"
            "repeated at least 5 times, like (CAG)5."
        ), 0, 0, 1, 4)

        MOTIF_NAMES = {
            1: "Mononucleotide", 2: "Dinucleotide", 3: "Trinucleotide",
            4: "Tetranucleotide", 5: "Pentanucleotide", 6: "Hexanucleotide",
        }
        self.motif_checks = {}
        self.min_repeats = {}
        self.min_repeat_labels = {}
        for row, k in enumerate([1, 2, 3, 4, 5, 6], start=1):
            cb = QCheckBox(MOTIF_NAMES[k])
            cb.setChecked(k in [3, 4])
            if k == 1:
                cb.setToolTip(
                    "Search for homopolymer runs (a single base repeated, e.g. "
                    "AAAAAA). Unticked by default. Note: a homopolymer run is "
                    "never counted as a fake repeat of a longer motif (e.g. "
                    "AAAA is never treated as dinucleotide \"AA\" repeated) "
                    "regardless of whether this option is ticked."
                )
            cb.stateChanged.connect(self._update_repeat_visibility)
            self.motif_checks[k] = cb
            sg1.addWidget(cb, row, 0)

            lbl = _lbl("Min repeats:")
            sp = QSpinBox(); sp.setRange(2, 20); sp.setValue(6 if k in (1, 2) else 5)
            sp.setMinimumWidth(70)
            self.min_repeats[k] = sp
            self.min_repeat_labels[k] = lbl
            sg1.addWidget(lbl, row, 1)
            sg1.addWidget(sp, row, 2)

        gap_row = 7
        sg1.addWidget(_lbl(
            "SSR search mode",
            "Perfect repeats only: standard single-motif SSRs (e.g. (CA)8).\n\n"
            "Compound loci only: merges adjacent perfect repeats of DIFFERENT\n"
            "motifs into one locus, e.g. (CA)8 immediately followed by (GA)6.\n"
            "Standalone perfect SSRs that don't have a compound partner are\n"
            "excluded from the results in this mode.\n\n"
            "Both: shows compound loci AND keeps every individual perfect\n"
            "repeat (including the ones that make up a compound locus) as\n"
            "its own row too, so you can compare designing primers against\n"
            "the compound locus vs. its individual components."
        ), gap_row, 0, 1, 2)
        self.ssr_mode_combo = QComboBox()
        self.ssr_mode_combo.addItem("Perfect repeats only", "perfect")
        self.ssr_mode_combo.addItem("Compound loci only", "compound")
        self.ssr_mode_combo.addItem("Both (perfect + compound)", "both")
        self.ssr_mode_combo.currentIndexChanged.connect(self._on_ssr_mode_changed)
        sg1.addWidget(self.ssr_mode_combo, gap_row, 2, 1, 2)

        self.compound_max_gap_lbl = _lbl(
            "Max gap between motifs (bp)",
            "Maximum bp between one repeat's end and the next repeat's start\n"
            "for them to be merged into one compound locus."
        )
        sg1.addWidget(self.compound_max_gap_lbl, gap_row + 1, 0, 1, 2)
        self.compound_max_gap = QSpinBox()
        self.compound_max_gap.setRange(0, 50)
        self.compound_max_gap.setValue(10)
        sg1.addWidget(self.compound_max_gap, gap_row + 1, 2)
        self._on_ssr_mode_changed()  # set initial enabled state

        # ----- Imperfect / interrupted repeat detection -----
        imperfect_row = gap_row + 2
        self.detect_imperfect_cb = QCheckBox("Detect imperfect / interrupted repeats")
        self.detect_imperfect_cb.setChecked(False)
        self.detect_imperfect_cb.setToolTip(
            "Extend perfect repeats through a short interruption, e.g.\n"
            "(CA)8-TT-(CA)5 becomes one imperfect locus instead of being\n"
            "missed or reported as two separate short repeats.\n\n"
            "This is a heuristic extension, not an alignment — in rare cases\n"
            "a different, unrelated repeat sitting immediately next to one\n"
            "being extended could be absorbed as if it were an interruption."
        )
        self.detect_imperfect_cb.stateChanged.connect(self._on_imperfect_toggled)
        sg1.addWidget(self.detect_imperfect_cb, imperfect_row, 0, 1, 4)

        self.imperfect_gap_lbl = _lbl(
            "Max interruption (bp)",
            "Maximum length of the non-matching gap tolerated inside an\n"
            "otherwise perfect repeat array."
        )
        sg1.addWidget(self.imperfect_gap_lbl, imperfect_row + 1, 0, 1, 2)
        self.imperfect_max_gap = QSpinBox()
        self.imperfect_max_gap.setRange(1, 20)
        self.imperfect_max_gap.setValue(3)
        sg1.addWidget(self.imperfect_max_gap, imperfect_row + 1, 2)

        self.imperfect_min_ext_lbl = _lbl(
            "Min repeats after interruption",
            "How many more full repeat units must immediately follow the\n"
            "interruption for it to be absorbed into the locus, rather than\n"
            "treated as the end of the repeat."
        )
        sg1.addWidget(self.imperfect_min_ext_lbl, imperfect_row + 2, 0, 1, 2)
        self.imperfect_min_ext = QSpinBox()
        self.imperfect_min_ext.setRange(1, 10)
        self.imperfect_min_ext.setValue(2)
        sg1.addWidget(self.imperfect_min_ext, imperfect_row + 2, 2)
        self._on_imperfect_toggled()  # set initial enabled state

        layout.addWidget(search_group)

        # ============================================================
        # GROUP 2 — Filtering & quality
        # ============================================================
        filter_group = QGroupBox("Filtering & Quality")
        sg2 = QGridLayout(filter_group)
        sg2.setSpacing(10)
        sg2.setColumnStretch(4, 1)

        sg2.addWidget(_lbl(
            "Min contig length (bp)", "Skip contigs shorter than this length."
        ), 0, 0)
        self.min_contig = QSpinBox()
        self.min_contig.setRange(0, 100000)
        self.min_contig.setValue(50)
        self.min_contig.setSingleStep(50)
        sg2.addWidget(self.min_contig, 0, 1)

        sg2.addWidget(_lbl(
            "Motif standardisation",
            "Controls how equivalent motifs are grouped.\n"
            "Level 4 groups all rotations and strand orientations."
        ), 0, 2)
        self.std_level = QComboBox()
        for lvl, desc in [
            (0, "0 — None"),
            (1, "1 — Smallest rotation"),
            (2, "2 — + reverse complement"),
            (3, "3 — Extended"),
            (4, "4 — Full canonical (recommended)"),
        ]:
            self.std_level.addItem(desc, lvl)
        self.std_level.setCurrentIndex(4)
        sg2.addWidget(self.std_level, 0, 3)

        self.search_rev = QCheckBox("Search reverse complement strand")
        self.search_rev.setChecked(False)
        sg2.addWidget(self.search_rev, 1, 0, 1, 2)

        import os as _os
        _cpu_count = _os.cpu_count() or 1
        sg2.addWidget(_lbl(
            "CPU cores to use",
            "Number of parallel worker processes for SSR scanning.\n\n"
            "GRACE automatically recommends a safe value based on your "
            "available memory and this genome's size once a genome is "
            "loaded — lower this manually if you still hit memory issues "
            "on very large genomes: each worker process holds its own "
            "copy of whatever contig it's currently scanning, so fewer "
            "workers means lower peak memory use at the cost of speed."
        ), 2, 0)
        self.n_workers = QSpinBox()
        self.n_workers.setRange(1, _cpu_count)
        self.n_workers.setValue(_cpu_count)
        sg2.addWidget(self.n_workers, 2, 1)

        self.n_workers_info = QLabel("")
        self.n_workers_info.setWordWrap(True)
        self.n_workers_info.setObjectName("role_secondary")
        self.n_workers_info.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        sg2.addWidget(self.n_workers_info, 3, 0, 1, 4)

        layout.addWidget(filter_group)

        # Run
        run_row = QHBoxLayout()
        self.run_btn = QPushButton("Run SSR Detection")
        self.run_btn.setObjectName("primary")
        self.run_btn.clicked.connect(self._run)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setVisible(False)
        self.cancel_btn.clicked.connect(self._cancel)
        self.status_label = QLabel("")
        self.status_label.setObjectName("role_secondary")
        self.status_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        run_row.addWidget(self.run_btn)
        run_row.addWidget(self.cancel_btn)
        run_row.addWidget(self.status_label)
        run_row.addStretch()
        layout.addLayout(run_row)

        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False); self.progress_bar.setFixedHeight(6)
        layout.addWidget(self.progress_bar)

        # Results
        self.results_group = QGroupBox("Results")
        self.results_group.setVisible(False)
        res_layout = QVBoxLayout(self.results_group)

        self.metrics_label = QLabel("")
        self.metrics_label.setObjectName("role_secondary")
        self.metrics_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        res_layout.addWidget(self.metrics_label)

        self.table = QTableWidget()
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSortingEnabled(True)
        self.table.setFixedHeight(500)
        res_layout.addWidget(self.table)

        sel_row = QHBoxLayout()
        self.sel_all_btn   = QPushButton("Select all")
        self.desel_all_btn = QPushButton("Deselect all")
        self.sel_all_btn.clicked.connect(lambda: self.table.selectAll())
        self.desel_all_btn.clicked.connect(lambda: self.table.clearSelection())
        self.sel_count_label = QLabel("")
        self.sel_count_label.setObjectName("role_secondary")
        self.sel_count_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        self.table.itemSelectionChanged.connect(self._on_selection_changed)
        self.download_btn = QPushButton("Download SSR results (CSV)")
        self.download_btn.clicked.connect(self._download_csv)
        sel_row.addWidget(self.sel_all_btn)
        sel_row.addWidget(self.desel_all_btn)
        sel_row.addWidget(self.sel_count_label)
        sel_row.addStretch()
        sel_row.addWidget(self.download_btn)
        res_layout.addLayout(sel_row)

        layout.addWidget(self.results_group)
        layout.addStretch()

        self._update_repeat_visibility()
        self._refresh()

    def _update_repeat_visibility(self):
        for k in [1, 2, 3, 4, 5, 6]:
            v = self.motif_checks[k].isChecked()
            self.min_repeat_labels[k].setVisible(v)
            self.min_repeats[k].setVisible(v)

    def _on_ssr_mode_changed(self):
        compound_relevant = self.ssr_mode_combo.currentData() != "perfect"
        self.compound_max_gap.setEnabled(compound_relevant)
        self.compound_max_gap_lbl.setEnabled(compound_relevant)

    def _on_imperfect_toggled(self):
        enabled = self.detect_imperfect_cb.isChecked()
        self.imperfect_max_gap.setEnabled(enabled)
        self.imperfect_gap_lbl.setEnabled(enabled)
        self.imperfect_min_ext.setEnabled(enabled)
        self.imperfect_min_ext_lbl.setEnabled(enabled)

    # ---------------------------------------------------------
    # RUN
    # ---------------------------------------------------------
    def _run(self):
        if self._worker is not None:
            # A scan is already in progress — ignore duplicate clicks rather
            # than reassigning self._worker and orphaning the still-running
            # QThread (which crashes with "QThread: Destroyed while thread
            # is still running").
            return
        if not self.state.has_genome:
            self._set_status("Load a genome first", "role_warning"); return
        motif_lengths = tuple(k for k, cb in self.motif_checks.items() if cb.isChecked())
        if not motif_lengths:
            self._set_status("Select at least one motif length", "role_warning"); return

        total_genome_size = sum(len(s) for s in self.state.genome.values())
        LARGE_GENOME_THRESHOLD = 1_000_000_000  # 1Gb
        if total_genome_size > LARGE_GENOME_THRESHOLD:
            from PyQt6.QtWidgets import QMessageBox
            msg = QMessageBox(self)
            msg.setWindowTitle("Large Genome Warning")
            msg.setText(
                f"This genome is ~{total_genome_size / 1e9:.1f}Gb.\n\n"
                f"SSR detection on genomes this size can take a while and use "
                f"significant memory. If it runs out of memory partway through, "
                f"try lowering 'CPU cores to use' below — fewer workers means "
                f"less memory used at once, at the cost of some speed."
            )
            msg.addButton("Continue anyway", QMessageBox.ButtonRole.AcceptRole)
            cancel = msg.addButton("Cancel", QMessageBox.ButtonRole.RejectRole)
            msg.exec()
            if msg.clickedButton() == cancel:
                return

        params = {
            "motif_lengths":        motif_lengths,
            "min_repeats":          {k: self.min_repeats[k].value() for k in motif_lengths},
            "search_reverse":       self.search_rev.isChecked(),
            "std_level":            self.std_level.currentData(),
            "min_contig_len":       self.min_contig.value(),
            "n_workers":            self.n_workers.value(),
        }
        self.run_btn.setEnabled(False)
        self.cancel_btn.setVisible(True)
        self.progress_bar.setVisible(True); self.progress_bar.setValue(0)
        self._set_status("Scanning...", "role_secondary")
        self.mw.set_status("Running SSR detection...")
        self._worker = SSRWorker(self.state.genome, params)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _cancel(self):
        if self._worker is not None:
            self._worker.terminate()
            self._worker.wait()
            self._cleanup_worker()
        if self._gff_worker is not None:
            self._gff_worker.terminate()
            self._gff_worker.wait()
            self._gff_worker = None
        self.run_btn.setEnabled(True)
        self.cancel_btn.setVisible(False)
        self.progress_bar.setVisible(False)
        self._set_status("Cancelled", "role_warning")
        self.mw.set_status("SSR detection cancelled")

    def _on_progress(self, done, total):
        self.progress_bar.setMaximum(total); self.progress_bar.setValue(done)
        self._set_status(f"Scanning contig {done:,} of {total:,}...", "role_secondary")

    def _on_done(self, ssrs, elapsed):
        try:
            mode = self.ssr_mode_combo.currentData()  # "perfect" | "compound" | "both"

            # Imperfect/interrupted extension runs FIRST — it only changes a
            # single locus's own boundaries, and compound merging's gap
            # logic (below) should see the final, possibly-extended
            # boundaries rather than the original tight perfect-seed ones.
            if self.detect_imperfect_cb.isChecked():
                from core.ssr_detection import extend_imperfect_ssrs
                ssrs = extend_imperfect_ssrs(
                    ssrs,
                    self.state.genome,
                    max_interruption=self.imperfect_max_gap.value(),
                    min_extension_repeats=self.imperfect_min_ext.value(),
                )

            if mode != "perfect":
                from core.ssr_detection import merge_compound_ssrs
                max_gap = self.compound_max_gap.value()
                ssrs = merge_compound_ssrs(ssrs, self.state.genome, max_gap=max_gap)

            if mode == "compound":
                # Only compound loci — standalone perfect SSRs that never found
                # a compound partner are excluded from the results entirely.
                ssrs = [s for s in ssrs if s.get("is_compound")]
            elif mode == "both":
                # Keep every compound locus AND re-expand its individual perfect
                # components back into their own rows too, so both views are
                # available side by side. Both the compound row and its
                # components share a "compound_group" number so they can be
                # traced back to each other in the table/CSV.
                expanded = []
                group_id = 0
                for s in ssrs:
                    if s.get("is_compound"):
                        group_id += 1
                        s = dict(s)
                        s["compound_group"] = group_id
                        expanded.append(s)
                        for comp in s["components"]:
                            comp = dict(comp)
                            comp["compound_group"] = group_id
                            expanded.append(comp)
                    else:
                        expanded.append(s)
                ssrs = expanded

            # Reassign ssr_id sequentially now that all processing is done
            ssrs.sort(key=lambda x: (x["contig"], x["start"]))
            for i, s in enumerate(ssrs, start=1):
                s["ssr_id"] = i
        except Exception as e:
            import traceback
            traceback.print_exc()
            self._cleanup_worker()
            self.run_btn.setEnabled(True)
            self.cancel_btn.setVisible(False)
            self.progress_bar.setVisible(False)
            self._set_status(f"Error processing SSR results: {e}", "role_error")
            self.mw.set_status("SSR detection failed during post-processing")
            return

        self.state.clear_downstream_of_ssrs()
        self.state.ssrs             = ssrs
        self.state.ssr_detection_time = elapsed
        self.state.selected_ssrs    = None
        self.state.ssr_version     += 1
        self.run_btn.setEnabled(True)
        self.cancel_btn.setVisible(False)
        self.progress_bar.setVisible(False)
        n_actual_compound = sum(1 for s in ssrs if s.get("is_compound"))
        n_actual_imperfect = sum(1 for s in ssrs if s.get("is_imperfect"))
        status = f"Found {len(ssrs):,} SSRs in {elapsed:.1f}s"
        if n_actual_compound:
            status += f" — {n_actual_compound:,} compound loci"
        if n_actual_imperfect:
            status += f" — {n_actual_imperfect:,} imperfect loci"
        self._set_status(status, "role_success")
        self.mw.set_status(f"SSR detection complete — {len(ssrs):,} SSRs found")
        self.mw.on_step_complete(1)
        self._cleanup_worker()
        self._refresh()
        if self.state.has_gff and ssrs:
            self._run_gff_annotation(ssrs)

    def _run_gff_annotation(self, ssrs):
        self._set_status(
            f"Found {len(ssrs):,} SSRs — annotating with GFF features...",
            "role_secondary"
        )
        self.mw.set_status("Annotating SSRs with genomic features...")
        self._gff_worker = GFFAnnotationWorker(
            ssrs, self.state.gff_path, self.state.gff_features
        )
        self._gff_worker.finished.connect(self._on_gff_done)
        self._gff_worker.error.connect(self._on_gff_error)
        self._gff_worker.start()

    def _on_gff_done(self, annotated_ssrs, gff_index, chrom_names):
        self.state.ssrs         = annotated_ssrs
        self.state.gff_features = gff_index
        self.state.chrom_names  = chrom_names if chrom_names else {}
        self.state.ssr_version += 1
        n = len(annotated_ssrs)
        chrom_note = f" — {len(chrom_names):,} chromosome names mapped" if chrom_names else ""
        self._set_status(f"Found {n:,} SSRs — genomic features annotated ✓{chrom_note}", "role_success")
        self.mw.set_status(f"SSR annotation complete — {n:,} SSRs annotated")
        if self._gff_worker:
            try:
                self._gff_worker.finished.disconnect()
                self._gff_worker.error.disconnect()
            except Exception:
                pass
            self._gff_worker = None
        self._refresh()

    def _on_gff_error(self, msg):
        self._set_status(f"SSRs found but GFF annotation failed: {msg}", "role_warning")
        self.mw.set_status("GFF annotation failed — SSRs available without feature labels")
        if self._gff_worker:
            try:
                self._gff_worker.finished.disconnect()
                self._gff_worker.error.disconnect()
            except Exception:
                pass
            self._gff_worker = None

    def _on_error(self, msg):
        self.run_btn.setEnabled(True)
        self.cancel_btn.setVisible(False)
        self.progress_bar.setVisible(False)
        self._set_status(f"Error: {msg}", "role_error")
        self.mw.set_status("SSR detection failed")
        self._cleanup_worker()

    def _cleanup_worker(self):
        if self._worker is not None:
            try:
                self._worker.progress.disconnect()
                self._worker.finished.disconnect()
                self._worker.error.disconnect()
            except Exception:
                pass
            # The thread has already finished by the time we're called
            # (from a finished/error signal handler) — this returns near-
            # instantly, but without it there's a real race: Qt's signal
            # firing and the underlying OS thread being FULLY marked as
            # stopped aren't the same instant, and releasing the last
            # Python reference (below) right in that window can trigger
            # "QThread: Destroyed while thread is still running".
            self._worker.wait()
            self._worker.genome = None
            self._worker.params = None
            self._worker = None

    def _set_status(self, msg, role):
        self.status_label.setText(msg)
        self.status_label.setObjectName(role)
        self.status_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")

    # ---------------------------------------------------------
    # REFRESH / TABLE
    # ---------------------------------------------------------
    def _update_worker_recommendation(self):
        if not self.state.has_genome:
            self.n_workers_info.setText("")
            return

        genome_id = id(self.state.genome)
        if genome_id == self._n_workers_genome_id:
            return  # already computed for this exact genome object
        self._n_workers_genome_id = genome_id

        from core.ssr_detection import recommend_worker_count, get_available_memory_bytes
        genome_size = sum(len(s) for s in self.state.genome.values())
        available = get_available_memory_bytes()
        recommended = recommend_worker_count(genome_size, available_memory_bytes=available)

        # Only auto-apply if the user hasn't manually changed the value
        # away from whatever we last auto-set (don't clobber their choice)
        if self._n_workers_auto_value is None or self.n_workers.value() == self._n_workers_auto_value:
            self.n_workers.setValue(recommended)
            self._n_workers_auto_value = recommended

        genome_gb = genome_size / 1e9
        if available is not None:
            available_gb = available / 1e9
            self.n_workers_info.setText(
                f"Auto-recommended: {recommended} core(s) — {genome_gb:.1f}Gb genome, "
                f"{available_gb:.1f}Gb memory available. Adjust manually if needed."
            )
        else:
            self.n_workers_info.setText(
                f"Recommended: {recommended} core(s) — {genome_gb:.1f}Gb genome "
                f"(install psutil for memory-aware recommendations)."
            )

    def _refresh(self):
        self._update_worker_recommendation()
        if not self.state.has_ssrs:
            self.results_group.setVisible(False); return
        if self.state.ssr_version == self._last_rendered_version:
            return
        self._last_rendered_version = self.state.ssr_version
        self.results_group.setVisible(True)

        import pandas as pd
        ssrs = self.state.ssrs
        df   = pd.DataFrame(ssrs)
        has_compound = "is_compound" in df.columns and df["is_compound"].any()
        has_imperfect = "imperfect_structure" in df.columns and df["imperfect_structure"].notna().any()
        if has_compound:
            df["is_compound"] = df["is_compound"].apply(lambda v: "Yes" if v is True else "")
        if has_imperfect:
            df["imperfect_structure"] = df["imperfect_structure"].fillna("")
        t    = self.state.ssr_detection_time
        truncated  = len(df) > TABLE_DISPLAY_LIMIT
        display_df = df.iloc[:TABLE_DISPLAY_LIMIT] if truncated else df

        metrics = (
            f"{len(ssrs):,} SSRs   |   {df['motif'].nunique():,} unique motifs   |   "
            f"{df['contig'].nunique():,} contigs" + (f"   |   {t:.1f}s" if t else "")
        )
        if truncated:
            metrics += f"   |   showing first {TABLE_DISPLAY_LIMIT:,} rows — download CSV for full results"
        self.metrics_label.setText(metrics)

        has_annotation = "genomic_feature" in df.columns
        has_group = "compound_group" in df.columns
        cols = ["ssr_id", "contig", "start", "end", "motif", "canonical_motif", "repeat_count"]
        if has_compound:
            cols.append("is_compound")
        if has_imperfect:
            cols.append("imperfect_structure")
        if has_group:
            cols.append("compound_group")
        if has_annotation:
            cols.append("genomic_feature")

        col_labels = {
            "ssr_id":          "SSR ID",
            "contig":          "Contig",
            "start":           "Start (bp)",
            "end":             "End (bp)",
            "motif":           "Motif",
            "canonical_motif": "Canonical motif",
            "repeat_count":    "Repeat count",
            "genomic_feature": "Genomic feature",
            "is_compound":     "Compound locus",
            "imperfect_structure": "Repeat structure (if imperfect)",
            "compound_group":  "Compound group",
        }
        display_cols = [c for c in cols if c in df.columns]

        # Add Chromosome column if we have mapping
        show_chromosome = hasattr(self.state, 'get_display_name') and self.state.chrom_names
        if show_chromosome:
            display_cols.append("chromosome")  # will be added virtually
            col_labels["chromosome"] = "Chromosome"

        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(display_df))
        self.table.setColumnCount(len(display_cols))
        self.table.setHorizontalHeaderLabels([col_labels.get(c, c) for c in display_cols])

        header = self.table.horizontalHeader()
        header.setStretchLastSection(True)

        for row_idx in range(len(display_df)):
            for col_idx, col in enumerate(display_cols):
                if col == "chromosome":
                    contig = display_df.iat[row_idx, display_df.columns.get_loc("contig")]
                    val = self.state.get_display_name(contig)
                else:
                    val = display_df.iat[row_idx, display_df.columns.get_loc(col)]
                
                # Use NumericTableWidgetItem for numeric columns
                if isinstance(val, (int, float)):
                    text = str(val)
                    item = NumericTableWidgetItem(text, numeric_value=val)
                else:
                    item = QTableWidgetItem(str(val))
                
                if col_idx == 0:
                    item.setData(
                        Qt.ItemDataRole.UserRole,
                        int(display_df.iat[row_idx, display_df.columns.get_loc("ssr_id")])
                    )
                self.table.setItem(row_idx, col_idx, item)

        # Column sizing runs AFTER population: ResizeToContents needs the
        # actual cell text to measure against, not the empty table — doing
        # this before population sizes columns against nothing, which is
        # exactly what was truncating longer values like the imperfect
        # repeat structure string (e.g. "(CA)4-t-(CA)8-tt-(CA)3").
        for col_idx, col in enumerate(display_cols):
            if col == "contig":
                header.setSectionResizeMode(col_idx, QHeaderView.ResizeMode.Interactive)
                self.table.setColumnWidth(col_idx, 140)
            elif col == "chromosome":
                header.setSectionResizeMode(col_idx, QHeaderView.ResizeMode.Interactive)
                self.table.setColumnWidth(col_idx, 120)
            elif col in ("left_primer", "right_primer"):
                header.setSectionResizeMode(col_idx, QHeaderView.ResizeMode.Interactive)
                self.table.setColumnWidth(col_idx, 200)
            else:
                header.setSectionResizeMode(col_idx, QHeaderView.ResizeMode.ResizeToContents)

        self.table.setSortingEnabled(True)

    def _on_selection_changed(self):
        rows = set(i.row() for i in self.table.selectedItems())
        if not rows or not self.state.has_ssrs:
            self.sel_count_label.setText("")
            self.state.selected_ssrs = None
            return
        ssr_id_set = set()
        for row in rows:
            item = self.table.item(row, 0)
            if item is not None:
                uid = item.data(Qt.ItemDataRole.UserRole)
                if uid is not None:
                    ssr_id_set.add(uid)
        if ssr_id_set:
            ssr_by_id = {s["ssr_id"]: s for s in self.state.ssrs}
            self.state.selected_ssrs = [
                ssr_by_id[i] for i in sorted(ssr_id_set) if i in ssr_by_id
            ]
        else:
            self.state.selected_ssrs = None
        count = len(self.state.selected_ssrs) if self.state.selected_ssrs else 0
        self.sel_count_label.setText(f"{count:,} SSRs selected" if count else "")

    def _download_csv(self):
        if not self.state.has_ssrs: return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save SSR results", "ssr_detection_results.csv", "CSV files (*.csv)"
        )
        if path:
            import pandas as pd
            df = pd.DataFrame(self.state.ssrs)
            if "components" in df.columns:
                df = df.drop(columns=["components"])
            if "interruptions" in df.columns:
                df = df.drop(columns=["interruptions"])
            if "is_compound" in df.columns:
                df["is_compound"] = df["is_compound"].apply(lambda v: "Yes" if v is True else "")
            if "is_imperfect" in df.columns:
                df = df.drop(columns=["is_imperfect"])
            if "imperfect_structure" in df.columns:
                df["imperfect_structure"] = df["imperfect_structure"].fillna("")
            if hasattr(self.state, 'get_display_name'):
                df["Chromosome"] = df["contig"].apply(lambda c: self.state.get_display_name(c))
            df.rename(columns={
                "ssr_id":          "SSR ID",
                "contig":          "Contig",
                "start":           "Start (bp)",
                "end":             "End (bp)",
                "motif":           "Motif",
                "canonical_motif": "Canonical motif",
                "repeat_count":    "Repeat count",
                "genomic_feature": "Genomic feature",
                "is_compound":     "Compound locus",
                "imperfect_structure": "Repeat structure (if imperfect)",
                "compound_group":  "Compound group",
            }).to_csv(path, index=False, encoding="utf-8-sig")
            self.mw.set_status(f"Saved to {path}")

    def on_show(self):
        self._refresh()