"""
home_panel.py — Load Genome + Optional GFF Annotation + Session Management

Supports FASTA, FASTQ, and GBFF (GenBank flat file) formats.
GBFF files contain both sequence and annotation — loading one automatically
populates both the genome and the annotation index.
"""
import os, json

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QGridLayout, QFileDialog, QMessageBox, QScrollArea,
)
from PyQt6.QtCore import QThread, pyqtSignal
from PyQt6.QtGui import QFont

from ui.style import (
    ACCENT, SUCCESS, ERROR, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING,
)

SAVEABLE_KEYS = [
    "ssrs", "primer_results", "specificity_results",
    "genome_filename", "product_min", "product_max",
    "filtered_primer_results",
    "gff_filename", "gff_path",
]


def _state_to_dict(state):
    payload = {}
    for k in SAVEABLE_KEYS:
        v = getattr(state, k, None)
        if v is not None:
            payload[k] = v.to_dict(orient="records") if hasattr(v, "to_dict") else v
    return payload


def _dict_to_state(state, payload):
    for k in SAVEABLE_KEYS:
        if k in payload:
            setattr(state, k, payload[k])
    if state.ssrs:
        state.ssr_version = 1
    if state.primer_results:
        state.filtered_primer_results = state.primer_results
        state.primer_version = 1
    if state.specificity_results:
        state.blast_version = 1


# ---------------------------------------------------------------------------
# Workers
# ---------------------------------------------------------------------------

class GenomeLoadWorker(QThread):
    """Loads FASTA / FASTQ files — sequence only."""
    finished = pyqtSignal(dict, str)
    error    = pyqtSignal(str)

    def __init__(self, path):
        super().__init__()
        self.path = path

    def run(self):
        try:
            from core.fasta_loader import load_sequence_file
            genome = load_sequence_file(self.path)
            self.finished.emit(genome, self.path)
        except Exception as e:
            self.error.emit(str(e))


class GBFFLoadWorker(QThread):
    """
    Loads a GBFF file — returns both genome dict AND GFFIndex.
    This is the key difference from GenomeLoadWorker: one file gives
    both sequence and annotation.
    """
    finished = pyqtSignal(dict, object, str)   # genome, gff_index, path
    error    = pyqtSignal(str)

    def __init__(self, path):
        super().__init__()
        self.path = path

    def run(self):
        try:
            from core.fasta_loader import load_gbff
            genome, gff_index = load_gbff(self.path, build_annotation=True)
            self.finished.emit(genome, gff_index, self.path)
        except Exception as e:
            self.error.emit(str(e))


class GFFLoadWorker(QThread):
    """Loads and parses a GFF3/GTF file into a GFFIndex."""
    finished = pyqtSignal(object, str)   # gff_index, path
    error    = pyqtSignal(str)

    def __init__(self, path, genome):
        super().__init__()
        self.path = path
        self.genome = genome

    def run(self):
        try:
            from core.gff_parser import build_gff_index
            gff_index = build_gff_index(self.path, genome=self.genome)
            self.finished.emit(gff_index, self.path)
        except Exception as e:
            self.error.emit(str(e))


# ---------------------------------------------------------------------------
# Panel
# ---------------------------------------------------------------------------

class HomePanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state   = state
        self.mw      = main_window
        self._worker            = None
        self._gff_worker        = None
        self._selected_path     = None
        self._selected_gff_path = None
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
        title = QLabel("Load Genome")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setObjectName("role_title")
        sub = QLabel(
            "Select a genome file — FASTA (.fa / .fasta / .fna), "
            "FASTQ (.fastq / .fq), or GenBank (.gb / .gbff / .gbk) supported"
        )
        sub.setObjectName("role_secondary")
        sub.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        sub.setWordWrap(True)
        L.addWidget(title)
        L.addWidget(sub)

        # ── Load genome ───────────────────────────────────
        load_group = QGroupBox("Genome File")
        lg = QVBoxLayout(load_group)
        lg.setSpacing(10)

        path_row = QHBoxLayout()
        self.path_label = QLabel("No file selected")
        self.path_label.setStyleSheet(
            f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;"
        )
        self.path_label.setWordWrap(True)
        self.browse_btn = QPushButton("Browse...")
        self.browse_btn.setFixedWidth(100)
        self.browse_btn.clicked.connect(self._browse)
        path_row.addWidget(self.path_label, 1)
        path_row.addWidget(self.browse_btn)
        lg.addLayout(path_row)

        # GBFF info note — shown when a GBFF file is selected
        self._gbff_note = QLabel(
            "✦ GenBank file detected — sequence and annotation will be loaded together. "
            "No separate GFF file needed."
        )
        self._gbff_note.setStyleSheet(
            f"color: {ACCENT}; font-size: {FONT_SIZE_SMALL}pt;"
        )
        self._gbff_note.setWordWrap(True)
        self._gbff_note.setVisible(False)
        lg.addWidget(self._gbff_note)

        self.load_btn = QPushButton("Load Genome")
        self.load_btn.setObjectName("primary")
        self.load_btn.setEnabled(False)
        self.load_btn.clicked.connect(self._load)
        lg.addWidget(self.load_btn)

        self.load_status = QLabel("")
        self.load_status.setStyleSheet(
            f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;"
        )
        lg.addWidget(self.load_status)
        L.addWidget(load_group)

        # ── Load GFF annotation (optional) ───────────────
        self._gff_group = QGroupBox("Genome Annotation (optional — not needed for GBFF)")
        gg = QVBoxLayout(self._gff_group)
        gg.setSpacing(10)

        gff_desc = QLabel(
            "Load a GFF3 or GTF annotation file to enable genomic context analysis "
            "(exon / intron / intergenic classification). "
            "Not required if you loaded a GBFF file — annotation is extracted automatically."
        )
        gff_desc.setObjectName("role_secondary")
        gff_desc.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        gff_desc.setWordWrap(True)
        gg.addWidget(gff_desc)

        gff_path_row = QHBoxLayout()
        self.gff_path_label = QLabel("No annotation file loaded")
        self.gff_path_label.setStyleSheet(
            f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;"
        )
        self.gff_path_label.setWordWrap(True)
        self.gff_browse_btn = QPushButton("Browse...")
        self.gff_browse_btn.setFixedWidth(100)
        self.gff_browse_btn.clicked.connect(self._browse_gff)
        gff_path_row.addWidget(self.gff_path_label, 1)
        gff_path_row.addWidget(self.gff_browse_btn)
        gg.addLayout(gff_path_row)

        # GFF file mismatch warning
        self._gff_mismatch_warning = QLabel(
            "⚠️ <b>Important: Matching Sequence and Annotation Files</b><br>"
            "The contig identifiers in your sequence file must match those in your "
            "annotation file for chromosome names and genomic features to be displayed "
            "correctly. Using mismatched files is a common source of confusion.<br><br>"
            "For best results, always download the sequence and annotation files together "
            "from the same assembly version. If you are using NCBI data, we recommend "
            "using the <b>RefSeq assembly</b> (accession starting with 'GCF_') for complete "
            "compatibility."
        )
        self._gff_mismatch_warning.setStyleSheet(
            f"background-color: {WARNING}22; color: {WARNING}; padding: 12px; "
            "border-radius: 4px; margin-top: 8px;"
        )
        self._gff_mismatch_warning.setWordWrap(True)
        self._gff_mismatch_warning.setVisible(True)  # always show as helpful reminder
        gg.addWidget(self._gff_mismatch_warning)

        gff_btn_row = QHBoxLayout()
        self.gff_load_btn = QPushButton("Load Annotation")
        self.gff_load_btn.setEnabled(False)
        self.gff_load_btn.clicked.connect(self._load_gff)
        self.gff_clear_btn = QPushButton("Clear")
        self.gff_clear_btn.setVisible(False)
        self.gff_clear_btn.clicked.connect(self._clear_gff)
        gff_btn_row.addWidget(self.gff_load_btn)
        gff_btn_row.addWidget(self.gff_clear_btn)
        gff_btn_row.addStretch()
        gg.addLayout(gff_btn_row)

        self.gff_status = QLabel("")
        self.gff_status.setStyleSheet(
            f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;"
        )
        gg.addWidget(self.gff_status)
        L.addWidget(self._gff_group)

        # ── Genome summary ────────────────────────────────
        self.summary_group = QGroupBox("Genome Summary")
        self.summary_group.setVisible(False)
        self.summary_layout = QGridLayout(self.summary_group)
        self.summary_layout.setSpacing(12)
        L.addWidget(self.summary_group)

        # ── Session management ────────────────────────────
        session_group = QGroupBox("Session Management")
        sg = QVBoxLayout(session_group)

        save_restore_row = QHBoxLayout()

        save_col = QVBoxLayout()
        save_lbl = QLabel("Save session")
        save_lbl.setFont(QFont(FONT_UI, FONT_SIZE_SMALL, QFont.Weight.Bold))
        save_lbl.setObjectName("role_primary_text")
        save_desc = QLabel("Write all results to a JSON file")
        save_desc.setObjectName("role_secondary")
        save_desc.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        self.save_btn = QPushButton("Save session...")
        self.save_btn.clicked.connect(self._save_session)
        self.save_btn.setToolTip(
            "Saves SSR results, primer results, and BLAST results to a JSON file.\n"
            "The genome file itself is not saved — you will need to reload it."
        )
        save_col.addWidget(save_lbl)
        save_col.addWidget(save_desc)
        save_col.addWidget(self.save_btn)
        save_restore_row.addLayout(save_col)
        save_restore_row.addSpacing(24)

        restore_col = QVBoxLayout()
        restore_lbl = QLabel("Restore session")
        restore_lbl.setFont(QFont(FONT_UI, FONT_SIZE_SMALL, QFont.Weight.Bold))
        restore_lbl.setObjectName("role_primary_text")
        restore_desc = QLabel("Load results from a previously saved JSON file")
        restore_desc.setObjectName("role_secondary")
        restore_desc.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        self.restore_btn = QPushButton("Restore session...")
        self.restore_btn.clicked.connect(self._restore_session)
        self.restore_btn.setToolTip(
            "Load a previously saved session.\n"
            "Load your genome first, then restore to continue where you left off."
        )
        restore_col.addWidget(restore_lbl)
        restore_col.addWidget(restore_desc)
        restore_col.addWidget(self.restore_btn)
        save_restore_row.addLayout(restore_col)
        save_restore_row.addStretch()
        sg.addLayout(save_restore_row)

        self.session_status = QLabel("")
        self.session_status.setStyleSheet(
            f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;"
        )
        self.session_status.setWordWrap(True)
        sg.addWidget(self.session_status)
        L.addWidget(session_group)

        # ── Reset ─────────────────────────────────────────
        self.reset_btn = QPushButton("Reset — clear all data and start fresh")
        self.reset_btn.setObjectName("role_error")
        self.reset_btn.setStyleSheet(f"border-color: {ERROR};")
        self.reset_btn.clicked.connect(self._reset)
        self.reset_btn.setVisible(False)
        L.addWidget(self.reset_btn)

        L.addStretch()
        self._refresh()

    # ---------------------------------------------------------
    # GENOME LOAD
    # ---------------------------------------------------------

    def _is_gbff(self, path: str) -> bool:
        return path.lower().endswith((".gb", ".gbff", ".gbk"))

    def _browse(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select genome file", "",
            "All supported (*.fa *.fasta *.fna *.fastq *.fq *.gb *.gbff *.gbk);;"
            "FASTA files (*.fa *.fasta *.fna);;"
            "FASTQ files (*.fastq *.fq);;"
            "GenBank files (*.gb *.gbff *.gbk);;"
            "All files (*)"
        )
        if path:
            self._selected_path = path
            self.path_label.setText(path)
            self.path_label.setStyleSheet(
                f"color: {TEXT_PRIMARY}; font-size: {FONT_SIZE_SMALL}pt;"
            )
            self.load_btn.setEnabled(True)
            # Show GBFF note when applicable
            is_gbff = self._is_gbff(path)
            self._gbff_note.setVisible(is_gbff)
            # Dim GFF section when GBFF selected (not needed)
            self._gff_group.setTitle(
                "Genome Annotation (not needed — annotation extracted from GBFF)"
                if is_gbff else
                "Genome Annotation (optional — not needed for GBFF)"
            )

    def _load(self):
        if not self._selected_path:
            return
        self.load_btn.setEnabled(False)
        self.browse_btn.setEnabled(False)
        self._set_load_status("Loading...", "role_secondary")
        self.mw.set_status("Loading genome...")

        if self._is_gbff(self._selected_path):
            self._worker = GBFFLoadWorker(self._selected_path)
            self._worker.finished.connect(self._on_gbff_load_done)
            self._worker.error.connect(self._on_load_error)
        else:
            self._worker = GenomeLoadWorker(self._selected_path)
            self._worker.finished.connect(self._on_load_done)
            self._worker.error.connect(self._on_load_error)

        self._worker.start()

    def _on_load_done(self, genome, path):
        """Callback for FASTA/FASTQ loads."""
        self.state.clear_downstream_of_genome()
        self.state.genome          = genome
        self.state.genome_path     = path
        self.state.genome_filename = os.path.basename(path)
        self._set_load_status(f"Loaded — {len(genome):,} sequences", "role_success")
        self.browse_btn.setEnabled(True)
        self.load_btn.setEnabled(True)
        self.mw.set_status(f"Genome loaded — {len(genome):,} sequences")
        self.mw.on_step_complete(0)
        self._refresh()

    def _on_gbff_load_done(self, genome, gff_index, path):
        """
        Callback for GBFF loads.
        Populates both genome AND annotation from the single file.
        """
        self.state.clear_downstream_of_genome()
        self.state.clear_gff()

        self.state.genome          = genome
        self.state.genome_path     = path
        self.state.genome_filename = os.path.basename(path)

        n_chrom = 0
        if gff_index and gff_index.n_features > 0:
            self.state.gff_path     = path   # same file
            self.state.gff_filename = os.path.basename(path)
            self.state.gff_features = gff_index
            self.state.chrom_names  = dict(gff_index.chrom_names)
            n_chrom = len(gff_index.chrom_names)

            self.gff_status.setText(
                f"✓ Annotation extracted from GBFF — "
                f"{gff_index.n_features:,} features, "
                f"{n_chrom:,} chromosome name{'s' if n_chrom != 1 else ''} mapped"
            )
            self.gff_status.setStyleSheet(
                f"color: {SUCCESS}; font-size: {FONT_SIZE_SMALL}pt;"
            )
            self.gff_clear_btn.setVisible(True)
        else:
            self.gff_status.setText(
                "GBFF loaded but no annotation features found in file."
            )
            self.gff_status.setStyleSheet(
                f"color: {WARNING}; font-size: {FONT_SIZE_SMALL}pt;"
            )

        n_seq = len(genome)
        msg   = f"GBFF loaded — {n_seq:,} sequences"
        if n_chrom:
            msg += f", {n_chrom:,} chromosomes named"
        self._set_load_status(msg, "role_success")
        self.browse_btn.setEnabled(True)
        self.load_btn.setEnabled(True)
        self.mw.set_status(msg)
        self.mw.on_step_complete(0)
        self._refresh()

    def _on_load_error(self, msg):
        self._set_load_status(f"Error: {msg}", "role_error")
        self.load_btn.setEnabled(True)
        self.browse_btn.setEnabled(True)
        self.mw.set_status("Failed to load genome")

    def _set_load_status(self, msg, role):
        self.load_status.setText(msg)
        self.load_status.setObjectName(role)
        self.load_status.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")

    # ---------------------------------------------------------
    # GFF LOAD
    # ---------------------------------------------------------
    def _browse_gff(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select annotation file", "",
            "Annotation files (*.gff *.gff3 *.gtf);;All files (*)"
        )
        if path:
            self._selected_gff_path = path
            self.gff_path_label.setText(path)
            self.gff_path_label.setStyleSheet(
                f"color: {TEXT_PRIMARY}; font-size: {FONT_SIZE_SMALL}pt;"
            )
            self.gff_load_btn.setEnabled(True)

    def _load_gff(self):
        if not self._selected_gff_path:
            return
        if not self.state.has_genome:
            self.gff_status.setText("Load a genome first before loading annotation")
            self.gff_status.setObjectName("role_warning")
            self.gff_status.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
            return

        self.gff_load_btn.setEnabled(False)
        self.gff_status.setText("Parsing annotation file...")
        self.gff_status.setObjectName("role_secondary")
        self.gff_status.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        self.mw.set_status("Loading annotation...")

        self._gff_worker = GFFLoadWorker(self._selected_gff_path, self.state.genome)
        self._gff_worker.finished.connect(self._on_gff_load_done)
        self._gff_worker.error.connect(self._on_gff_load_error)
        self._gff_worker.start()

    def _on_gff_load_done(self, gff_index, path):
        self.state.clear_gff()
        self.state.gff_path     = path
        self.state.gff_filename = os.path.basename(path)
        self.state.gff_features = gff_index
        self.state.chrom_names  = dict(gff_index.chrom_names) if gff_index else {}

        n_chrom = len(self.state.chrom_names)
        self.gff_status.setText(
            f"✓ Annotation loaded — {gff_index.n_features:,} features, "
            f"{n_chrom:,} chromosome name{'s' if n_chrom != 1 else ''} mapped"
        )
        self.gff_status.setObjectName("role_success")
        self.gff_status.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        self.gff_clear_btn.setVisible(True)
        self.gff_load_btn.setEnabled(True)
        self.mw.set_status(f"Annotation loaded — {self.state.gff_filename}")
        self._refresh()

    def _on_gff_load_error(self, msg):
        self.gff_status.setText(f"Error: {msg}")
        self.gff_status.setObjectName("role_error")
        self.gff_status.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
        self.gff_load_btn.setEnabled(True)
        self.mw.set_status("Failed to load annotation")

    def _clear_gff(self):
        self.state.clear_gff()
        self._selected_gff_path = None
        self.gff_path_label.setText("No annotation file loaded")
        self.gff_path_label.setStyleSheet(
            f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;"
        )
        self.gff_load_btn.setEnabled(False)
        self.gff_clear_btn.setVisible(False)
        self.gff_status.setText("")
        self.mw.set_status("Annotation cleared")
        self._refresh()

    # ---------------------------------------------------------
    # SESSION SAVE / RESTORE
    # ---------------------------------------------------------
    def _save_session(self):
        default = os.path.join(
            os.path.dirname(self.state.genome_path)
            if self.state.genome_path else os.path.expanduser("~"),
            "grace_session.json"
        )
        path, _ = QFileDialog.getSaveFileName(
            self, "Save session", default, "JSON files (*.json)"
        )
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as f:
                json.dump(_state_to_dict(self.state), f, default=str, indent=2)
            self._set_session_status(f"Session saved to {path}", "role_success")
            self.mw.set_status("Session saved")
        except Exception as e:
            self._set_session_status(f"Failed to save: {e}", "role_error")

    def _restore_session(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Restore session", os.path.expanduser("~"), "JSON files (*.json)"
        )
        if not path:
            return
        self._load_session_file(path)

    def _load_session_file(self, path):
        try:
            with open(path, "r", encoding="utf-8") as f:
                payload = json.load(f)
            _dict_to_state(self.state, payload)
            parts = []
            if self.state.ssrs:                parts.append(f"{len(self.state.ssrs):,} SSRs")
            if self.state.primer_results:      parts.append(f"{len(self.state.primer_results):,} primer pairs")
            if self.state.specificity_results: parts.append(f"{len(self.state.specificity_results):,} specificity results")
            summary = ", ".join(parts) if parts else "no result data"
            self._set_session_status(f"Session restored — {summary}", "role_success")
            if self.state.has_genome:      self.mw.on_step_complete(0)
            if self.state.has_ssrs:        self.mw.on_step_complete(1)
            if self.state.has_primers:     self.mw.on_step_complete(2)
            if self.state.has_specificity: self.mw.on_step_complete(3)
            self.mw.refresh_sidebar()
            self.mw.set_status("Session restored")
            self._refresh()
        except Exception as e:
            self._set_session_status(f"Failed to restore: {e}", "role_error")

    def _set_session_status(self, msg, role):
        self.session_status.setText(msg)
        self.session_status.setObjectName(role)
        self.session_status.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")

    # ---------------------------------------------------------
    # RESET
    # ---------------------------------------------------------
    def _reset(self):
        reply = QMessageBox.question(
            self, "Reset", "Clear all data and start fresh?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if reply == QMessageBox.StandardButton.Yes:
            self.state.clear_downstream_of_genome()
            self.state.clear_gff()
            self.state.genome = self.state.genome_path = self.state.genome_filename = None
            self._selected_path     = None
            self._selected_gff_path = None
            self.path_label.setText("No file selected")
            self.path_label.setStyleSheet(
                f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;"
            )
            self.gff_path_label.setText("No annotation file loaded")
            self.gff_path_label.setStyleSheet(
                f"color: {TEXT_SECONDARY}; font-size: {FONT_SIZE_SMALL}pt;"
            )
            self._gbff_note.setVisible(False)
            self.load_btn.setEnabled(False)
            self.gff_load_btn.setEnabled(False)
            self.gff_clear_btn.setVisible(False)
            self.load_status.setText("")
            self.gff_status.setText("")
            self.mw.set_status("Reset")
            self.mw.refresh_sidebar()
            self._refresh()

    # ---------------------------------------------------------
    # REFRESH
    # ---------------------------------------------------------
    def _refresh(self):
        has = self.state.has_genome
        self.summary_group.setVisible(has)
        self.reset_btn.setVisible(has)

        if self.state.has_gff and not self.gff_status.text():
            n_chrom = len(self.state.chrom_names or {})
            self.gff_status.setText(
                f"✓ Annotation loaded — {self.state.gff_filename} "
                f"({n_chrom:,} chromosome names mapped)"
            )
            self.gff_status.setObjectName("role_success")
            self.gff_status.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
            self.gff_clear_btn.setVisible(True)

        if has:
            genome   = self.state.genome
            total_bp = sum(len(s) for s in genome.values())
            while self.summary_layout.count():
                item = self.summary_layout.takeAt(0)
                if item.widget(): item.widget().deleteLater()

            is_gbff = self.state.genome_filename and \
                      self.state.genome_filename.lower().endswith((".gb", ".gbff", ".gbk"))

            metrics = [
                ("File",            self.state.genome_filename or "—"),
                ("Format",          "GenBank (GBFF)" if is_gbff else "FASTA/FASTQ"),
                ("Sequences",       f"{len(genome):,}"),
                ("Total length",    f"{total_bp:,} bp"),
                ("Shortest",        f"{min(len(s) for s in genome.values()):,} bp"),
                ("Longest",         f"{max(len(s) for s in genome.values()):,} bp"),
                ("Annotation",      self.state.gff_filename or "None loaded"),
                ("Chr names",       f"{len(self.state.chrom_names or {}):,} mapped"
                                    if self.state.chrom_names else "None"),
            ]
            for row, (label, value) in enumerate(metrics):
                lbl = QLabel(label)
                lbl.setObjectName("role_secondary")
                lbl.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
                val = QLabel(value)
                val.setFont(QFont(FONT_MONO, FONT_SIZE_NORMAL))
                val.setObjectName("role_primary_text")
                self.summary_layout.addWidget(lbl, row // 2, (row % 2) * 2)
                self.summary_layout.addWidget(val, row // 2, (row % 2) * 2 + 1)

    def on_show(self):
        self._refresh()