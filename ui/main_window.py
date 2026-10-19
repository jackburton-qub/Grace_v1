"""
main_window.py
--------------
QMainWindow — sidebar navigation + stacked panel area.
"""

import os
from datetime import datetime

from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QHBoxLayout, QVBoxLayout,
    QListWidget, QListWidgetItem, QStackedWidget,
    QLabel, QStatusBar,
    QButtonGroup, QRadioButton,
    QPushButton, QFileDialog, QMessageBox,
)
from PyQt6.QtCore import Qt, QSize
from PyQt6.QtGui import QFont, QPixmap, QShortcut, QKeySequence

from ui.style import (
    STYLESHEET, build_stylesheet,
    SIDEBAR_WIDTH, WINDOW_MIN_W, WINDOW_MIN_H,
    ACCENT, SUCCESS, MUTED, TEXT_SECONDARY,
    BG_MID, BORDER,
    FONT_UI, FONT_MONO, FONT_SIZE_NORMAL,
    FONT_SIZE_PRESET_DEFAULT, DEFAULT_THEME,
)
from ui.panels.home_panel         import HomePanel
from ui.panels.ssr_panel          import SSRPanel
from ui.panels.ssr_summary_panel  import SSRSummaryPanel
from ui.panels.primer_panel       import PrimerPanel
from ui.panels.specificity_panel  import SpecificityPanel
from ui.panels.results_panel      import ResultsPanel
from ui.panels.panel_design_panel import PanelDesignPanel

STEPS = [
    (0,  "01  Load Genome"),
    (1,  "02  Detect SSRs"),
    (2,  "03  SSR Summary"),
    (3,  "04  Design Primers"),
    (4,  "05  Check Specificity"),
    (5,  "06  BLAST Results"),
    (6,  "07  Panel Design"),
]


class MainWindow(QMainWindow):
    def __init__(self, state):
        super().__init__()
        self.state = state
        self.setWindowTitle("GRACE — Genomic Repeat Analysis and Characterisation Engine")
        # Increased minimum height to accommodate 10 sidebar items
        self.setMinimumSize(WINDOW_MIN_W, max(WINDOW_MIN_H, 800))
        self.setStyleSheet(STYLESHEET)
        self._build_ui()
        self._connect_signals()
        self.refresh_sidebar()

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # ── Sidebar ──────────────────────────────────────
        sidebar_container = QWidget()
        sidebar_container.setFixedWidth(SIDEBAR_WIDTH)
        sidebar_container.setObjectName("sidebar_container")
        sidebar_layout = QVBoxLayout(sidebar_container)
        sidebar_layout.setContentsMargins(0, 0, 0, 0)
        sidebar_layout.setSpacing(0)

        header = QWidget()
        header.setFixedHeight(64)
        header.setObjectName("sidebar_header")
        header_layout = QVBoxLayout(header)
        header_layout.setContentsMargins(16, 12, 16, 12)
        header_layout.setSpacing(2)

        title_label = QLabel("GRACE")
        title_label.setFont(QFont(FONT_UI, 13, QFont.Weight.Bold))
        title_label.setObjectName("role_title")

        sub_label = QLabel("Genomic Repeat Analysis")
        sub_label.setFont(QFont(FONT_MONO, 7))
        sub_label.setObjectName("role_secondary")

        header_layout.addWidget(title_label)
        header_layout.addWidget(sub_label)

        self.sidebar = QListWidget()
        self.sidebar.setObjectName("sidebar")
        self.sidebar.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.sidebar.setSpacing(2)

        for _, label in STEPS:
            item = QListWidgetItem(label)
            item.setSizeHint(QSize(SIDEBAR_WIDTH, 48))
            item.setFont(QFont(FONT_UI, FONT_SIZE_NORMAL))
            self.sidebar.addItem(item)

        # Font size selector
        font_size_widget = QWidget()
        font_size_widget.setObjectName("sidebar_footer_section")
        font_size_layout = QVBoxLayout(font_size_widget)
        font_size_layout.setContentsMargins(16, 8, 16, 4)
        font_size_layout.setSpacing(6)

        font_size_label = QLabel("TEXT SIZE")
        font_size_label.setObjectName("role_section_label")
        font_size_layout.addWidget(font_size_label)

        btn_row = QWidget()
        btn_row.setStyleSheet("background: transparent; border: none;")
        btn_row_layout = QHBoxLayout(btn_row)
        btn_row_layout.setContentsMargins(0, 0, 0, 0)
        btn_row_layout.setSpacing(4)

        self._font_btn_group = QButtonGroup(self)
        self._font_preset = FONT_SIZE_PRESET_DEFAULT

        for preset in ("Medium", "Large"):
            btn = QRadioButton(preset)
            btn.setObjectName("sidebar_radio")
            if preset == FONT_SIZE_PRESET_DEFAULT:
                btn.setChecked(True)
            self._font_btn_group.addButton(btn)
            btn_row_layout.addWidget(btn)

        font_size_layout.addWidget(btn_row)

        theme_label = QLabel("THEME")
        theme_label.setObjectName("role_section_label")
        font_size_layout.addWidget(theme_label)

        theme_btn_row = QWidget()
        theme_btn_row.setStyleSheet("background: transparent; border: none;")
        theme_btn_row_layout = QHBoxLayout(theme_btn_row)
        theme_btn_row_layout.setContentsMargins(0, 0, 0, 0)
        theme_btn_row_layout.setSpacing(4)

        self._theme_btn_group = QButtonGroup(self)
        self._theme = DEFAULT_THEME

        for theme_name, label in (("dark", "Dark"), ("light", "Light")):
            btn = QRadioButton(label)
            btn.setObjectName("sidebar_radio")
            btn.setProperty("theme_name", theme_name)
            if theme_name == DEFAULT_THEME:
                btn.setChecked(True)
            self._theme_btn_group.addButton(btn)
            theme_btn_row_layout.addWidget(btn)

        self._theme_btn_group.buttonClicked.connect(self._on_theme_changed)
        font_size_layout.addWidget(theme_btn_row)

        session_widget = QWidget()
        session_widget.setObjectName("sidebar_footer_section")
        session_layout = QVBoxLayout(session_widget)
        session_layout.setContentsMargins(16, 8, 16, 4)
        session_layout.setSpacing(6)

        session_label = QLabel("SESSION")
        session_label.setObjectName("role_section_label")
        session_layout.addWidget(session_label)

        self.save_session_btn = QPushButton("Save Session")
        self.save_session_btn.clicked.connect(self._on_save_session)
        self.load_session_btn = QPushButton("Load Session")
        self.load_session_btn.clicked.connect(self._on_load_session)
        session_layout.addWidget(self.save_session_btn)
        session_layout.addWidget(self.load_session_btn)

        report_widget = QWidget()
        report_widget.setObjectName("sidebar_footer_section")
        report_layout = QVBoxLayout(report_widget)
        report_layout.setContentsMargins(16, 8, 16, 4)
        report_layout.setSpacing(6)

        report_label = QLabel("REPORT")
        report_label.setObjectName("role_section_label")
        report_layout.addWidget(report_label)

        self.export_report_btn = QPushButton("Export Report (PDF)")
        self.export_report_btn.setToolTip(
            "Generates a PDF summary of everything completed so far — genome "
            "stats, SSR charts, primer design, BLAST specificity, and panel "
            "design, whichever of these you've actually run. Available at "
            "any stage, not just after a full pipeline run."
        )
        self.export_report_btn.clicked.connect(self._on_export_report)
        report_layout.addWidget(self.export_report_btn)

        version_label = QLabel("v1.0.0")
        version_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        version_label.setObjectName("role_muted")
        version_label.setStyleSheet("font-size: 7pt; padding: 8px;")

        sidebar_layout.addWidget(header)
        sidebar_layout.addWidget(self.sidebar)
        sidebar_layout.addStretch()
        sidebar_layout.addWidget(font_size_widget)
        sidebar_layout.addWidget(session_widget)
        sidebar_layout.addWidget(report_widget)
        sidebar_layout.addWidget(version_label)

        # ── Main panel stack ─────────────────────────────
        self.stack = QStackedWidget()

        self.home_panel        = HomePanel(self.state, self)
        self.ssr_panel         = SSRPanel(self.state, self)
        self.ssr_summary_panel = SSRSummaryPanel(self.state, self)
        self.primer_panel      = PrimerPanel(self.state, self)
        self.specificity_panel = SpecificityPanel(self.state, self)
        self.results_panel     = ResultsPanel(self.state, self)
        self.panel_design_panel = PanelDesignPanel(self.state, self)

        for panel in [
            self.home_panel, self.ssr_panel, self.ssr_summary_panel,
            self.primer_panel, self.specificity_panel, self.results_panel,
            self.panel_design_panel,
        ]:
            self.stack.addWidget(panel)

        root.addWidget(sidebar_container)
        root.addWidget(self.stack, 1)

        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready")

    def _connect_signals(self):
        self.sidebar.currentRowChanged.connect(self._on_tab_changed)
        self.sidebar.setCurrentRow(0)
        self._font_btn_group.buttonClicked.connect(self._on_font_size_changed)
        # Hidden shortcut: high-resolution screenshot of the whole window.
        QShortcut(QKeySequence("Ctrl+Shift+P"), self, activated=self._save_screenshot)

    def _on_tab_changed(self, index):
        self.stack.setCurrentIndex(index)
        panel = self.stack.currentWidget()
        if hasattr(panel, "on_show"):
            panel.on_show()

    def _on_font_size_changed(self, btn):
        preset = btn.text()
        if preset == self._font_preset:
            return
        self._font_preset = preset
        self.setStyleSheet(build_stylesheet(preset, self._theme))

    def _on_theme_changed(self, btn):
        theme_name = btn.property("theme_name")
        if theme_name == self._theme:
            return
        self._theme = theme_name
        self.setStyleSheet(build_stylesheet(self._font_preset, theme_name))
        # QPainter-based widgets (the SSR Summary charts) don't respond to
        # the QSS stylesheet at all — they need to be told directly.
        if hasattr(self, "ssr_summary_panel"):
            self.ssr_summary_panel.set_theme(theme_name)

    def _save_screenshot(self):
        scale = 3
        pixmap = QPixmap(self.size() * scale)
        pixmap.setDevicePixelRatio(scale)
        pixmap.fill(Qt.GlobalColor.transparent)
        self.render(pixmap)

        out_dir = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "screenshots"
        )
        os.makedirs(out_dir, exist_ok=True)
        path = os.path.join(
            out_dir, f"grace_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
        )
        if pixmap.save(path, "PNG"):
            self.set_status(f"Screenshot saved to {path}")
        else:
            self.set_status(f"Could not save screenshot to {path}")

    def set_status(self, message: str):
        self.status_bar.showMessage(message)

    def _on_save_session(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Session", "grace_session.json", "GRACE Session (*.json)"
        )
        if not path:
            return
        try:
            self.state.save_session(path)
            self.set_status(f"Session saved to {path}")
        except Exception as e:
            QMessageBox.critical(self, "Save Failed", f"Could not save session:\n\n{e}")

    def _on_load_session(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Session", "", "GRACE Session (*.json)"
        )
        if not path:
            return
        try:
            # Mutates self.state in place — every panel already holds a
            # reference to this exact object, so we must not replace it
            # with a new instance.
            self.state.load_session_into(path)
        except FileNotFoundError as e:
            QMessageBox.warning(self, "Genome File Missing", str(e))
            return
        except Exception as e:
            QMessageBox.critical(self, "Load Failed", f"Could not load session:\n\n{e}")
            return

        # Force every panel to re-render against the newly loaded state.
        # Several panels only re-render when their own "last rendered
        # version" stamp differs from state's current version — since a
        # loaded session's version numbers could coincidentally match
        # whatever a panel already has cached from earlier in this same
        # app run, that check alone can't be trusted here, so refresh
        # everything unconditionally instead.
        for panel in [
            self.home_panel, self.ssr_panel, self.ssr_summary_panel,
            self.primer_panel, self.specificity_panel, self.results_panel,
            self.panel_design_panel,
        ]:
            if hasattr(panel, "_refresh"):
                panel._refresh()
            elif hasattr(panel, "on_show"):
                panel.on_show()

        self.refresh_sidebar()
        self.set_status(f"Session loaded from {path}")

    def _on_export_report(self):
        from ui.report_generator import export_report_dialog
        export_report_dialog(self.state, self)

    def refresh_sidebar(self):
        from PyQt6.QtGui import QColor
        s = self.state

        statuses = [
            s.has_genome,       # 0 Load Genome
            s.has_ssrs,         # 1 Detect SSRs
            s.has_ssrs,         # 2 SSR Summary
            s.has_primers,      # 3 Design Primers
            s.has_specificity,  # 4 Check Specificity
            s.has_specificity,  # 5 BLAST Results
            s.has_panel_design, # 6 Panel Design
        ]

        for i, (_, label) in enumerate(STEPS):
            item = self.sidebar.item(i)
            if item is None:
                continue
            is_current = self.sidebar.currentRow() == i
            done = statuses[i]

            if done:
                item.setForeground(QColor(SUCCESS))
            elif is_current:
                item.setForeground(QColor(ACCENT))
            else:
                item.setForeground(QColor(TEXT_SECONDARY))

    def navigate_to(self, step: int):
        self.sidebar.setCurrentRow(step)

    def on_step_complete(self, step: int):
        self.refresh_sidebar()
        self.set_status(f"Step {step + 1} complete")