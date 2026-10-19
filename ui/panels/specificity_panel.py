"""
specificity_panel.py — Specificity Check
Single scroll area. Runs BLAST. Stores results. Shows pass count only.
BLAST binaries are now bundled with the app — no path input needed.
"""
import os, sys, shutil, tempfile

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QGridLayout, QSpinBox, QDoubleSpinBox, QComboBox,
    QProgressBar, QScrollArea, QCheckBox, QFileDialog, QMessageBox,
)
from PyQt6.QtCore import QThread, pyqtSignal, QTimer
from PyQt6.QtGui import QFont

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_SIZE_LARGE, FONT_SIZE_SMALL, PANEL_PADDING,
)


from ui.table_utils import secondary_label as _lbl
from core.primer_import import MIN_PRIMER_LEN, MAX_PRIMER_LEN


class BLASTWorker(QThread):
    finished = pyqtSignal(list, str)  # (results, path to the raw-hits tab-separated file)
    error    = pyqtSignal(str)
    progress = pyqtSignal(int, int)  # (primers_done, primers_total)

    def __init__(self, genome, primers, blast_params, specificity_params):
        super().__init__()
        self.genome             = genome
        self.primers            = primers
        self.blast_params       = blast_params
        self.specificity_params = specificity_params

    def run(self):
        try:
            from core.primer_specificity_blast import check_specificity_blast, _get_or_write_genome_fasta

            # Always reconstruct a genuine FASTA file from the already-
            # parsed genome dict, rather than assuming the ORIGINAL file
            # the user loaded is FASTA — it might be GBFF or FASTQ, both
            # supported load formats that look nothing like FASTA and
            # would make BLAST's makeblastdb fail with "does NOT appear
            # to be FASTA" if passed straight through. Written to a
            # STABLE, content-keyed path rather than a fresh temp
            # directory each run, so repeat runs on the same genome can
            # actually reuse the cached BLAST database built from it.
            genome_fasta_path = _get_or_write_genome_fasta(self.genome)

            # A stable, predictable location for this run's raw hit
            # detail, generated here (rather than left to
            # check_specificity_blast()'s own auto-generated fallback)
            # specifically so this worker KNOWS where the file ended up
            # afterward without needing it handed back across the thread
            # boundary at all — only the path string travels with the
            # finished signal, not any of the actual hit data.
            raw_hits_path = os.path.join(
                tempfile.gettempdir(), f"grace_raw_blast_hits_{os.getpid()}_{id(self)}.tsv"
            )

            results = check_specificity_blast(
                genome_fasta=genome_fasta_path,
                primer_results=self.primers,
                blast_params=self.blast_params,
                specificity_params=self.specificity_params,
                progress_callback=lambda done, total: self.progress.emit(done, total),
                raw_hits_path=raw_hits_path,
                # blast_bin_dir omitted — resolved automatically from bundled blast/
            )
            from core.primer_import import apply_blast_amplicon_sizes
            apply_blast_amplicon_sizes(results)
            # results never carries raw hit records or full amplicon
            # detail at all now — every surviving hit was written
            # straight to raw_hits_path incrementally, batch by batch, as
            # the run progressed, rather than accumulated in memory for
            # the whole run only to be processed all at once at the end.
            self.finished.emit(results, raw_hits_path)
        except Exception as e:
            self.error.emit(str(e))



class SpecificityPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state   = state
        self.mw      = main_window
        self._worker = None
        self._run_start_time = None
        self._primers_done = 0
        self._primers_total = 0
        self._elapsed_timer = QTimer(self)
        self._elapsed_timer.setInterval(1000)  # tick every second
        self._elapsed_timer.timeout.connect(self._update_elapsed_time)
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
        title = QLabel("Specificity Check")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setObjectName("role_title")
        sub = QLabel("BLAST primers against the genome to verify unique binding")
        sub.setObjectName("role_secondary")
        sub.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(title)
        L.addWidget(sub)

        # BLAST status (bundled — read only, no path input)
        blast_group = QGroupBox("BLAST+ Status")
        bg = QVBoxLayout(blast_group)
        self.blast_detected_label = QLabel("")
        self.blast_detected_label.setObjectName("role_secondary")
        self.blast_detected_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        bg.addWidget(self.blast_detected_label)
        L.addWidget(blast_group)

        # Primer queue
        queue_group = QGroupBox("Primer Queue")
        qg = QVBoxLayout(queue_group)
        self.queue_label = QLabel("No primers available")
        self.queue_label.setObjectName("role_secondary")
        self.queue_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        qg.addWidget(self.queue_label)

        import_row = QHBoxLayout()
        self.import_btn = QPushButton("Import primers\u2026")
        self.import_btn.setToolTip(
            "Screen existing primer pairs without running SSR detection or "
            "primer design. Accepts:\n"
            "  \u2022 CSV or TSV with a header row: name, forward, reverse, and "
            "optionally expected_size (bp)\n"
            "  \u2022 FASTA with each pair as two entries named <name>_F and <name>_R\n\n"
            f"Sequences must be A/C/G/T only and {MIN_PRIMER_LEN}-{MAX_PRIMER_LEN} bp. "
            "Imported pairs are BLASTed with the same settings as designed ones, "
            "show up on the BLAST Results page with IDs starting 'IMP-', and can "
            "be used in panel design. Where no expected_size is given, the size "
            "of the single on-target product BLAST predicts is used instead.\n\n"
            "Importing replaces any previously imported pairs and clears "
            "existing BLAST results, since the screening queue has changed."
        )
        self.import_btn.clicked.connect(self._import_primers)
        self.clear_import_btn = QPushButton("Remove imported")
        self.clear_import_btn.setToolTip("Remove all imported primer pairs from the queue.")
        self.clear_import_btn.clicked.connect(self._clear_imported)
        import_row.addWidget(self.import_btn)
        import_row.addWidget(self.clear_import_btn)
        import_row.addStretch()
        qg.addLayout(import_row)

        self.import_info = QLabel("")
        self.import_info.setWordWrap(True)
        self.import_info.setObjectName("role_secondary")
        self.import_info.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        qg.addWidget(self.import_info)
        L.addWidget(queue_group)

        # Specificity criteria — one unified set of controls instead of a
        # binary Amplicons/Hits preset choice. The base check (exactly one
        # correctly-sized PCR product) always applies; the two checkboxes
        # below are independent, optional strictness layers on top of it.
        mode_group = QGroupBox("Specificity Criteria")
        mg = QVBoxLayout(mode_group)
        mg.setSpacing(10)

        base_desc = QLabel(
            "Base check (always applied): exactly one PCR product must form "
            "within your product size range."
        )
        base_desc.setWordWrap(True)
        base_desc.setObjectName("role_secondary")
        base_desc.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        mg.addWidget(base_desc)

        self.require_unique_cb = QCheckBox("Require unique primer binding (strict)")
        self.require_unique_cb.setChecked(True)
        self.require_unique_cb.setToolTip(
            "Adds an extra requirement on top of the base check: each primer "
            "must map to at most 'Max locations per primer' genomic sites, "
            "checked at the same mismatch/identity thresholds below.\n\n"
            "This is the stricter of the two toggles here — since forming any "
            "off-target amplicon inherently means a primer has a second binding "
            "site, turning this ON means the 'tolerate off-target' option below "
            "rarely gets a chance to matter: this check will usually reject a "
            "secondary site before size is even considered."
        )
        self.require_unique_cb.stateChanged.connect(self._on_unique_toggled)
        mg.addWidget(self.require_unique_cb)

        unique_row = QHBoxLayout()
        self.max_hits_label = _lbl("Max locations per primer",
            "Maximum genomic locations each primer may map to and still pass.")
        unique_row.addWidget(self.max_hits_label)
        self.max_hits = QSpinBox(); self.max_hits.setRange(1, 10); self.max_hits.setValue(1)
        unique_row.addWidget(self.max_hits)
        unique_row.addStretch()
        mg.addLayout(unique_row)

        self.allow_offtarget_cb = QCheckBox("Tolerate large off-target products")
        self.allow_offtarget_cb.setChecked(True)
        self.allow_offtarget_cb.setToolTip(
            "If a secondary product could theoretically form elsewhere in the "
            "genome, allow it as long as it's larger than the size below — a "
            "product that much bigger essentially never co-amplifies alongside "
            "your intended short product under standard PCR conditions.\n\n"
            "When 'Require unique primer binding' above is also ON, this mostly "
            "only matters for a primer's extra binding site that doesn't pair "
            "into any amplicon at all — most genuine off-target products get "
            "caught by the unique-binding check first, regardless of size."
        )
        self.allow_offtarget_cb.stateChanged.connect(self._on_offtarget_toggled)
        mg.addWidget(self.allow_offtarget_cb)

        offtarget_row = QHBoxLayout()
        offtarget_row.addWidget(_lbl(
            "Max tolerated off-target size (bp)",
            "An off-target product larger than this is assumed unlikely to "
            "compete with the much shorter intended amplicon under standard "
            "PCR conditions, and won't cause the locus to fail specificity "
            "on its own."
        ))
        self.max_offtarget = QSpinBox(); self.max_offtarget.setRange(100, 20000)
        self.max_offtarget.setSingleStep(100); self.max_offtarget.setValue(1500)
        offtarget_row.addWidget(self.max_offtarget)
        offtarget_row.addStretch()
        mg.addLayout(offtarget_row)

        self.require_3prime_cb = QCheckBox("Require mismatches at primer's 3' end to ignore a hit")
        self.require_3prime_cb.setChecked(True)
        self.require_3prime_cb.setToolTip(
            "Whole-primer identity/mismatch counts can't tell you WHERE a "
            "mismatch sits — and that's what actually determines whether a "
            "hit is a real PCR competitor, since a mismatch at or near the 3' "
            "end is what stops polymerase extension, while a mismatch further "
            "away often doesn't matter at all. By default, a hit is only considered safe to "
            "ignore once it has AT LEAST the number of mismatches set below "
            "within the terminal window — a single 3' mismatch alone is NOT "
            "enough to clear a hit."
        )
        self.require_3prime_cb.stateChanged.connect(self._on_3prime_toggled)
        mg.addWidget(self.require_3prime_cb)

        threeprime_row = QHBoxLayout()
        self.three_prime_label = _lbl("3' end window (bp)")
        threeprime_row.addWidget(self.three_prime_label)
        self.three_prime_window = QSpinBox(); self.three_prime_window.setRange(1, 15)
        self.three_prime_window.setValue(5)
        threeprime_row.addWidget(self.three_prime_window)

        self.three_prime_min_mismatches_label = _lbl(
            "Min mismatches in window to ignore a hit",
            "Default: a hit needs at least 2 mismatches within "
            "the 3' end window before it's considered safe to ignore. A hit "
            "with only 0 or 1 mismatch there still counts as a real competitor."
        )
        threeprime_row.addWidget(self.three_prime_min_mismatches_label)
        self.three_prime_min_mismatches = QSpinBox(); self.three_prime_min_mismatches.setRange(1, 5)
        self.three_prime_min_mismatches.setValue(2)
        threeprime_row.addWidget(self.three_prime_min_mismatches)
        threeprime_row.addStretch()
        mg.addLayout(threeprime_row)

        L.addWidget(mode_group)

        # Settings
        settings_group = QGroupBox("BLAST Settings")
        sg = QGridLayout(settings_group)
        sg.setSpacing(10); sg.setContentsMargins(12, 12, 12, 12)

        sg.addWidget(_lbl("Max mismatches",
            "Maximum mismatches allowed in a BLAST alignment to count as a hit."), 0, 0)
        self.max_mismatches = QSpinBox(); self.max_mismatches.setRange(0, 10); self.max_mismatches.setValue(5)
        sg.addWidget(self.max_mismatches, 0, 1)

        sg.addWidget(_lbl("Min alignment length (bp)",
            "Minimum alignment length to count as a significant hit."), 0, 2)
        self.min_align_len = QSpinBox(); self.min_align_len.setRange(5, 30); self.min_align_len.setValue(18)
        sg.addWidget(self.min_align_len, 0, 3)

        sg.addWidget(_lbl("Min % identity",
            "Minimum percentage identity for a hit to be counted."), 1, 0)
        self.min_identity = QDoubleSpinBox(); self.min_identity.setRange(50, 100); self.min_identity.setValue(70.0); self.min_identity.setSingleStep(1)
        sg.addWidget(self.min_identity, 1, 1)

        sg.addWidget(_lbl("BLAST task",
            "blastn-short: optimised for short sequences like primers (recommended).\nblastn: standard nucleotide BLAST."), 1, 2)
        self.blast_task = QComboBox()
        self.blast_task.addItem("blastn-short (recommended for primers)", "blastn-short")
        self.blast_task.addItem("blastn", "blastn")
        sg.addWidget(self.blast_task, 1, 3)

        sg.addWidget(_lbl("Word size",
            "BLAST word size. 7 is recommended for short primer sequences — "
            "raising it reduces sensitivity to real off-target sites rather "
            "than you deciding to tolerate them."), 2, 0)
        self.word_size = QSpinBox(); self.word_size.setRange(4, 28); self.word_size.setValue(7)
        sg.addWidget(self.word_size, 2, 1)

        sg.addWidget(_lbl("E-value",
            "Maximum E-value for a BLAST hit to be reported. The looser default (0.01) "
            "lets BLAST report candidate hits broadly; precision is enforced by the "
            "mismatch/identity thresholds above, not by starving BLAST's own search."), 2, 2)
        self.evalue = QDoubleSpinBox(); self.evalue.setDecimals(4); self.evalue.setRange(0.0001, 30000.0); self.evalue.setSingleStep(1.0); self.evalue.setValue(10.0)
        sg.addWidget(self.evalue, 2, 3)

        sg.addWidget(_lbl("Dust filter",
            "DUST masks low-complexity regions before BLAST alignment.\n"
            "Recommended: OFF — primers are short (18-25bp); masking even a few "
            "bases of one is a large fraction of it, risking primers being wrongly "
            "flagged as unable to find their own true binding site."), 3, 0)
        self.dust_filter = QComboBox()
        self.dust_filter.addItem("OFF (recommended)", "no")
        self.dust_filter.addItem("ON", "yes")
        sg.addWidget(self.dust_filter, 3, 1, 1, 3)

        L.addWidget(settings_group)

        self.threads_info = QLabel("")
        self.threads_info.setWordWrap(True)
        self.threads_info.setObjectName("role_secondary")
        self.threads_info.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(self.threads_info)

        # Run button
        run_row = QHBoxLayout()
        self.run_btn = QPushButton("Run BLAST")
        self.run_btn.setFixedHeight(40)
        self.run_btn.clicked.connect(self._run)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setFixedHeight(40)
        self.cancel_btn.setVisible(False)
        self.cancel_btn.clicked.connect(self._cancel)
        run_row.addWidget(self.run_btn)
        run_row.addWidget(self.cancel_btn)
        run_row.addStretch()
        L.addLayout(run_row)

        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setFixedHeight(6)
        self.progress_bar.setRange(0, 0)  # indeterminate
        L.addWidget(self.progress_bar)

        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)  # import/size-mismatch notes can run long
        self.status_label.setObjectName("role_secondary")
        self.status_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        L.addWidget(self.status_label)

        L.addStretch()

        self._on_unique_toggled()      # set initial enabled state
        self._on_offtarget_toggled()   # set initial enabled state
        self._on_3prime_toggled()      # set initial enabled state
        self._refresh()

    # ---------------------------------------------------------
    # CRITERIA TOGGLE HANDLERS
    # ---------------------------------------------------------
    def _on_unique_toggled(self):
        enabled = self.require_unique_cb.isChecked()
        self.max_hits.setEnabled(enabled)
        self.max_hits_label.setEnabled(enabled)

    def _on_offtarget_toggled(self):
        enabled = self.allow_offtarget_cb.isChecked()
        self.max_offtarget.setEnabled(enabled)

    def _on_3prime_toggled(self):
        enabled = self.require_3prime_cb.isChecked()
        self.three_prime_window.setEnabled(enabled)
        self.three_prime_label.setEnabled(enabled)
        self.three_prime_min_mismatches.setEnabled(enabled)
        self.three_prime_min_mismatches_label.setEnabled(enabled)

    # ---------------------------------------------------------
    # RUN
    # ---------------------------------------------------------
    def _on_progress(self, done, total):
        self._primers_done = done
        self._primers_total = total

    def _update_elapsed_time(self):
        if self._run_start_time is None:
            return
        import time
        elapsed = int(time.time() - self._run_start_time)
        mins, secs = divmod(elapsed, 60)
        self._set_status(
            f"Running BLAST\u2026 {mins}m {secs}s elapsed. "
            f"Primers complete: {self._primers_done:,} / {self._primers_total:,}",
            "role_secondary"
        )

    def _run(self):
        if not self.state.has_genome:
            self._set_status("Load a genome first", "role_warning"); return
        if not self.state.has_primers and not self.state.has_imported_primers:
            self._set_status("Run primer design or import primers first", "role_warning"); return

        # Snapshot primers — local variable, immune to any state changes during run
        _primers = list(self.state.blast_primers)
        if not _primers:
            self._set_status("No primers to BLAST", "role_warning"); return

        genome_size = sum(len(s) for s in self.state.genome.values())
        LARGE_GENOME_THRESHOLD = 1_000_000_000  # 1Gb
        MANY_PRIMERS_THRESHOLD = 200            # primer pairs
        if genome_size > LARGE_GENOME_THRESHOLD or len(_primers) > MANY_PRIMERS_THRESHOLD:
            from PyQt6.QtWidgets import QMessageBox
            msg = QMessageBox(self)
            msg.setWindowTitle("BLAST May Take a While")
            msg.setText(
                f"This run has {len(_primers):,} primer pairs against a "
                f"~{genome_size / 1e9:.1f}Gb genome.\n\n"
                f"For reference, a similar-scale run (~530 primer pairs) took "
                f"roughly 28 minutes via an external BLAST service. Expect a "
                f"comparable or longer runtime here, depending on your hardware.\n\n"
                f"Elapsed time will be shown throughout the run."
            )
            msg.addButton("Continue anyway", QMessageBox.ButtonRole.AcceptRole)
            cancel = msg.addButton("Cancel", QMessageBox.ButtonRole.RejectRole)
            msg.exec()
            if msg.clickedButton() == cancel:
                return

        from core.primer_specificity_blast import BlastParams, SpecificityParams

        blast_params = BlastParams(
            task=self.blast_task.currentData(),
            word_size=self.word_size.value(),
            evalue=self.evalue.value(),
            dust=self.dust_filter.currentData(),
        )

        sp = SpecificityParams(
            min_product=self.state.product_min,
            max_product=self.state.product_max,
            max_mismatches=self.max_mismatches.value(),
            min_align_len=self.min_align_len.value(),
            min_identity=self.min_identity.value(),
            allow_offtarget_amplicons=self.allow_offtarget_cb.isChecked(),
            max_offtarget_amplicon_size=self.max_offtarget.value(),
            require_unique_hits=self.require_unique_cb.isChecked(),
            max_total_hits=self.max_hits.value(),
            require_3prime_match=self.require_3prime_cb.isChecked(),
            three_prime_window=self.three_prime_window.value(),
            three_prime_min_mismatches_to_ignore=self.three_prime_min_mismatches.value(),
        )

        self.run_btn.setEnabled(False)
        self.cancel_btn.setVisible(True)
        self.progress_bar.setVisible(True)
        self._set_status(f"Running BLAST on {len(_primers):,} primers...", "role_secondary")
        self.mw.set_status("Running BLAST...")

        import time
        self._run_start_time = time.time()
        self._primers_done = 0
        self._primers_total = len(_primers)
        self._elapsed_timer.start()

        self._worker = BLASTWorker(
            self.state.genome, _primers, blast_params, sp
        )
        self._worker.finished.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.progress.connect(self._on_progress)
        self._worker.start()

    def _on_done(self, results, raw_hits_path):
        self.state.specificity_results = results
        self.state.blast_raw_hits_path = raw_hits_path
        self.state.blast_version += 1
        self.mw.refresh_sidebar()

        n_pass = sum(1 for r in results if r.get("specificity_status") == "PASS")
        n_fail = len(results) - n_pass
        n_unique_pass = len({r.get("ssr_id") for r in results if r.get("specificity_status") == "PASS"})

        self.run_btn.setEnabled(True)
        self.cancel_btn.setVisible(False)
        self.progress_bar.setVisible(False)
        # Imported pairs whose stated expected_size doesn't match the product
        # BLAST predicts — worth a look before using them in a panel.
        size_mismatch = [
            r.get("primer_name") for r in results
            if r.get("source") == "imported" and r.get("expected_size")
            and r.get("blast_amplicon_size") and r["expected_size"] != r["blast_amplicon_size"]
        ]
        mismatch_note = (
            f" {len(size_mismatch):,} imported pair(s) have an expected_size that differs from "
            f"the BLAST-predicted product size (see the BLAST Results page)."
            if size_mismatch else ""
        )
        self._set_status(
            f"BLAST complete — {n_pass:,} PASS ({n_unique_pass:,} unique SSRs), {n_fail:,} FAIL. "
            f"Navigate to BLAST Results to view.{mismatch_note}",
            "role_warning" if size_mismatch else "role_success"
        )
        self.mw.set_status(f"Specificity check complete — {n_unique_pass:,} unique SSRs passed")
        self.mw.on_step_complete(4)
        self._cleanup_worker()

    def _on_error(self, msg):
        self.run_btn.setEnabled(True)
        self.cancel_btn.setVisible(False)
        self.progress_bar.setVisible(False)
        if "not found" in msg.lower() or "no such file" in msg.lower():
            self._set_status("BLAST+ executables not found — please report this issue.", "role_error")
        else:
            self._set_status(f"BLAST failed: {msg}", "role_error")
        self.mw.set_status("Specificity check failed")
        self._cleanup_worker()

    def _cancel(self):
        if self._worker is not None:
            self._worker.terminate()
            self._worker.wait()
            self._cleanup_worker()
        self.run_btn.setEnabled(True)
        self.cancel_btn.setVisible(False)
        self.progress_bar.setVisible(False)
        self._set_status("Cancelled", "role_warning")
        self.mw.set_status("Specificity check cancelled")

    # ---------------------------------------------------------
    # IMPORT
    # ---------------------------------------------------------
    def _import_primers(self):
        if self._worker is not None:
            return  # don't change the queue mid-run
        path, _ = QFileDialog.getOpenFileName(
            self, "Import primer pairs", "",
            "Primer files (*.csv *.tsv *.txt *.fa *.fasta *.fas *.fna);;All files (*)"
        )
        if not path:
            return
        from core.primer_import import import_primer_file
        result = import_primer_file(path)
        name = os.path.basename(path)

        if not result.ok:
            shown = result.errors[:15]
            more = len(result.errors) - len(shown)
            msg = QMessageBox(self)
            msg.setIcon(QMessageBox.Icon.Warning)
            msg.setWindowTitle("Import Failed")
            msg.setText(
                f"{name} wasn't imported \u2014 {len(result.errors)} problem(s) found. "
                f"Fix them and import the file again; nothing was changed."
            )
            msg.setInformativeText("\n".join(shown) + (f"\n\u2026and {more} more (see details)." if more else ""))
            msg.setDetailedText("\n".join(result.errors + result.warnings))
            msg.exec()
            self.mw.set_status(f"Primer import failed \u2014 {len(result.errors)} problem(s) in {name}")
            return

        self.state.imported_primers = result.pairs
        self.state.clear_downstream_of_primers()
        self.mw.refresh_sidebar()

        warnings = list(result.warnings)
        lo, hi = self.state.product_min, self.state.product_max
        outside = [p["primer_name"] for p in result.pairs
                   if p["expected_size"] and not lo <= p["expected_size"] <= hi]
        if outside:
            listed = ", ".join(outside[:10]) + ("\u2026" if len(outside) > 10 else "")
            warnings.append(
                f"{len(outside)} pair(s) have an expected_size outside the {lo}-{hi} bp "
                f"product size range (set on the Primer Design page). Screening only "
                f"recognises an on-target product inside that range, so these pairs "
                f"will fail unless BLAST finds some other product inside it: {listed}"
            )
        if warnings:
            QMessageBox.information(
                self, "Imported With Warnings",
                f"Imported {len(result.pairs):,} primer pair(s) from {name}.\n\n" + "\n".join(warnings)
            )
        self._import_source = name
        self._set_status("", "role_secondary")
        self.mw.set_status(f"Imported {len(result.pairs):,} primer pair(s) from {name}")
        self._refresh()

    def _clear_imported(self):
        if self._worker is not None or not self.state.has_imported_primers:
            return
        self.state.imported_primers = None
        self.state.clear_downstream_of_primers()
        self.mw.refresh_sidebar()
        self._set_status("", "role_secondary")
        self.mw.set_status("Imported primer pairs removed")
        self._refresh()

    def _cleanup_worker(self):
        """Disconnect signals and null large refs to prevent idle crashes."""
        self._elapsed_timer.stop()
        self._run_start_time = None
        if self._worker is not None:
            try:
                self._worker.finished.disconnect()
                self._worker.error.disconnect()
                self._worker.progress.disconnect()
            except Exception:
                pass
            # See ssr_panel.py's _cleanup_worker() for why this matters
            # even though the thread has already finished by this point.
            self._worker.wait()
            self._worker.genome             = None
            self._worker.primers            = None
            self._worker.blast_params       = None
            self._worker.specificity_params = None
            self._worker = None

    def _set_status(self, msg, role):
        self.status_label.setText(msg)
        self.status_label.setObjectName(role)
        self.status_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")

    # ---------------------------------------------------------
    # REFRESH
    # ---------------------------------------------------------
    def _update_threads_recommendation(self):
        if not self.state.has_genome:
            self.threads_info.setText("")
            return

        from core.primer_specificity_blast import recommend_blast_threads
        from core.ssr_detection import get_available_memory_bytes
        genome_size = sum(len(s) for s in self.state.genome.values())
        available = get_available_memory_bytes()
        recommended = recommend_blast_threads(genome_size, available_memory_bytes=available)

        genome_gb = genome_size / 1e9
        if available is not None:
            available_gb = available / 1e9
            self.threads_info.setText(
                f"BLAST will automatically use {recommended} thread(s) — {genome_gb:.1f}Gb genome, "
                f"{available_gb:.1f}Gb memory available. Scales down automatically on tighter memory."
            )
        else:
            self.threads_info.setText(
                f"BLAST will automatically use {recommended} thread(s) for this {genome_gb:.1f}Gb genome "
                f"(install psutil for memory-aware recommendations)."
            )

    def _refresh(self):
        self._update_threads_recommendation()

        # Show bundled BLAST status
        from core.primer_specificity_blast import get_blast_bin_dir
        blast_dir = get_blast_bin_dir()
        exe_suffix = ".exe" if sys.platform == "win32" else ""
        blastn_path = os.path.join(blast_dir, f"blastn{exe_suffix}") if blast_dir else ""

        if blast_dir and os.path.isfile(blastn_path):
            self.blast_detected_label.setText(f"✓ Bundled BLAST+ ready — {blastn_path}")
            self.blast_detected_label.setObjectName("role_success")
            self.blast_detected_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        elif shutil.which("blastn"):
            found = shutil.which("blastn")
            self.blast_detected_label.setText(f"✓ BLAST+ detected on system PATH — {found}")
            self.blast_detected_label.setObjectName("role_success")
            self.blast_detected_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        else:
            self.blast_detected_label.setText("✗ BLAST+ not found — please report this issue")
            self.blast_detected_label.setObjectName("role_error")
            self.blast_detected_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")

        # Queue info
        n_imported = len(self.state.imported_primers or [])
        imported_note = f" + {n_imported:,} imported" if n_imported else ""
        if self.state.has_primers:
            _filtered = self.state.filtered_primer_results
            _all      = self.state.primer_results or []
            n_blast   = len(_filtered) if _filtered else len(_all)
            n_total   = len(_all)
            if n_blast == n_total:
                self.queue_label.setText(f"{n_blast:,} primer pairs queued for BLAST{imported_note} — no quality filters applied")
            else:
                self.queue_label.setText(
                    f"{n_blast:,} primer pairs queued{imported_note} "
                    f"({n_total - n_blast:,} excluded by quality filters on Primer Design page)"
                )
            self.queue_label.setObjectName("role_primary_text")
            self.queue_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        elif n_imported:
            self.queue_label.setText(f"{n_imported:,} imported primer pairs queued for BLAST")
            self.queue_label.setObjectName("role_primary_text")
            self.queue_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        else:
            self.queue_label.setText("No primers available — run Primer Design or import primers")
            self.queue_label.setObjectName("role_secondary")
            self.queue_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")

        self.clear_import_btn.setVisible(bool(n_imported))
        if n_imported:
            n_sized = sum(1 for p in self.state.imported_primers if p.get("expected_size"))
            source = getattr(self, "_import_source", None)
            from_note = f" from {source}" if source else ""
            self.import_info.setText(
                f"{n_imported:,} imported pair(s){from_note} ({n_sized:,} with expected_size). "
                f"BLAST counts a product as on-target only if it is "
                f"{self.state.product_min}-{self.state.product_max} bp "
                f"(the product size range on the Primer Design page)."
            )
        else:
            self.import_info.setText("")

        # Show last run result if available
        if self.state.has_specificity and not self.status_label.text():
            n_pass = sum(1 for r in self.state.specificity_results if r.get("specificity_status") == "PASS")
            n_unique_pass = len({r.get("ssr_id") for r in self.state.specificity_results if r.get("specificity_status") == "PASS"})
            self._set_status(
                f"Last run: {n_pass:,} primers passed ({n_unique_pass:,} unique SSRs) — "
                f"navigate to BLAST Results to view, or run again.",
                "role_secondary"
            )

    def on_show(self):
        self._refresh()