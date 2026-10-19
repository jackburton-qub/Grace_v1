"""
style.py
--------
GRACE colour palette and Qt stylesheet.
IBM Plex Sans + IBM Plex Mono. Two themes: "dark" (GitHub dark palette
with Irish green accent, the original/default) and "light" (equivalent
GitHub-light-style palette, same green accent identity).

Theme switching mechanism:
    Every panel imports colour constants by name, e.g.
    `from ui.style import ACCENT, BG_MID`. Because Python copies the
    value at import time, reassigning ui.style.ACCENT later does NOT
    update files that already did `from ui.style import ACCENT` — so a
    runtime toggle can't work by mutating these constants after the fact.

    Instead, build_stylesheet(preset, theme) takes the theme explicitly
    and builds the ENTIRE Qt stylesheet string using that theme's colours
    baked directly into the returned string. Re-applying that string via
    self.setStyleSheet(...) on the main window re-themes every widget
    styled through the global QSS (backgrounds, buttons, inputs, tables,
    borders, tooltips, etc.) immediately, without needing any file to
    re-import anything.

    Known limitation: a handful of individual widgets set their OWN
    per-widget stylesheet directly at construction time using a literal
    f-string with one of these colour constants baked in (e.g.
    `label.setStyleSheet(f"color: {ACCENT};")` for a title or status
    message). Those widgets' colours were fixed at construction time and
    will not retroactively update when the theme toggles — this affects
    a small number of accent-coloured labels, not the overall look of the
    app, and is resolved by restarting GRACE after changing theme. Full
    correction would require converting every such call site to
    object-name-based QSS selectors instead of inline f-strings, which is
    a larger, separate cleanup.
"""

# ── Theme palettes ────────────────────────────────────────
THEMES = {
    "dark": {
        "BG_DARK":    "#0d1117",   # main window background
        "BG_MID":     "#161b22",   # panel / sidebar background
        "BG_LIGHT":   "#21262d",   # card / input / button background
        "BG_HOVER":   "#30363d",   # hover state

        "ACCENT":        "#2ea043",   # primary green accent
        "ACCENT_DIM":    "#238636",   # darker green for borders / pressed
        "ACCENT_BRIGHT": "#3fb950",   # bright green for highlights

        "SUCCESS": "#3fb950",   # PASS
        "ERROR":   "#f85149",   # FAIL
        "WARNING": "#d29922",   # amber
        "MUTED":   "#484f58",   # very muted

        "TEXT_PRIMARY":   "#e6edf3",   # main text
        "TEXT_SECONDARY": "#8b949e",   # labels / help text
        "TEXT_MONO":      "#79c0ff",   # monospace data values

        "BORDER":        "#21262d",    # standard border
        "BORDER_ACTIVE":  "#2ea043",   # focused / active border

        # Derived shades used directly in the stylesheet (table selection
        # tint, success-button background/hover, primary-button pressed)
        "SELECTION_BG":     "#1f3d2a",
        "SUCCESS_BTN_BG":   "#0f2d16",
        "SUCCESS_BTN_HOVER": "#1a4d26",
        "PRIMARY_PRESSED":  "#196127",
    },
    "light": {
        "BG_DARK":  "#fafbfc",   # main window background (soft off-white, not harsh pure white)
        "BG_MID":   "#f0f2f5",   # panel / sidebar background
        "BG_LIGHT": "#ffffff",   # card / input / button background (white surfaces read fine against the off-white canvas)
        "BG_HOVER": "#e5e8eb",   # hover state

        "ACCENT":        "#1a7f37",   # primary green accent (darker for contrast on light bg)
        "ACCENT_DIM":    "#116329",   # darker green for borders / pressed
        "ACCENT_BRIGHT": "#1a7f37",   # bright green for highlights

        "SUCCESS": "#1a7f37",   # PASS
        "ERROR":   "#cf222e",   # FAIL
        "WARNING": "#9a6700",   # amber (darkened for contrast on light bg)
        "MUTED":   "#8c959f",   # very muted

        "TEXT_PRIMARY":   "#1f2328",   # main text (near-black)
        "TEXT_SECONDARY": "#59636e",   # labels / help text
        "TEXT_MONO":      "#0969da",   # monospace data values

        "BORDER":        "#d0d7de",    # standard border
        "BORDER_ACTIVE":  "#1a7f37",   # focused / active border

        "SELECTION_BG":     "#d7f0dd",
        "SUCCESS_BTN_BG":   "#dafbe1",
        "SUCCESS_BTN_HOVER": "#c3f0cb",
        "PRIMARY_PRESSED":  "#0e4c1e",
    },
}

DEFAULT_THEME = "dark"

# ── Module-level constants ───────────────────────────────
# Reflect DEFAULT_THEME so existing `from ui.style import ACCENT`-style
# imports elsewhere in the app keep working unchanged. These do NOT
# update on theme switch (see module docstring) — build_stylesheet()
# below is the actual mechanism a runtime toggle must use.
_active = THEMES[DEFAULT_THEME]
BG_DARK, BG_MID, BG_LIGHT, BG_HOVER = _active["BG_DARK"], _active["BG_MID"], _active["BG_LIGHT"], _active["BG_HOVER"]
ACCENT, ACCENT_DIM, ACCENT_BRIGHT = _active["ACCENT"], _active["ACCENT_DIM"], _active["ACCENT_BRIGHT"]
SUCCESS, ERROR, WARNING, MUTED = _active["SUCCESS"], _active["ERROR"], _active["WARNING"], _active["MUTED"]
TEXT_PRIMARY, TEXT_SECONDARY, TEXT_MONO = _active["TEXT_PRIMARY"], _active["TEXT_SECONDARY"], _active["TEXT_MONO"]
BORDER, BORDER_ACTIVE = _active["BORDER"], _active["BORDER_ACTIVE"]

# ── Fonts ─────────────────────────────────────────────────
FONT_UI   = "IBM Plex Sans"
FONT_MONO = "IBM Plex Mono"

# Font size presets — base size drives all other sizes proportionally
FONT_PRESETS = {
    "Medium": {"normal": 11,  "small": 10,  "large": 13, "title": 16},
    "Large":  {"normal": 13, "small": 12, "large": 15, "title": 18},
}
FONT_SIZE_PRESET_DEFAULT = "Medium"

# Resolved at import time using the default preset — used by widgets that
# import these constants directly (sidebar font, etc.)
FONT_SIZE_NORMAL = FONT_PRESETS[FONT_SIZE_PRESET_DEFAULT]["normal"]
FONT_SIZE_SMALL  = FONT_PRESETS[FONT_SIZE_PRESET_DEFAULT]["small"]
FONT_SIZE_LARGE  = FONT_PRESETS[FONT_SIZE_PRESET_DEFAULT]["large"]
FONT_SIZE_TITLE  = FONT_PRESETS[FONT_SIZE_PRESET_DEFAULT]["title"]

# ── Dimensions ───────────────────────────────────────────
SIDEBAR_WIDTH  = 210
WINDOW_MIN_W   = 1200
WINDOW_MIN_H   = 800
PANEL_PADDING  = 24
CARD_RADIUS    = 6
BUTTON_RADIUS  = 6
BUTTON_HEIGHT  = 32


# ── Stylesheet builder ────────────────────────────────────
def build_stylesheet(preset: str = FONT_SIZE_PRESET_DEFAULT, theme: str = DEFAULT_THEME) -> str:
    """
    Build and return the full Qt stylesheet for the given font size preset
    and theme.
    preset must be one of: "Medium", "Large"
    theme must be one of: "dark", "light"
    """
    sizes = FONT_PRESETS.get(preset, FONT_PRESETS[FONT_SIZE_PRESET_DEFAULT])
    _n = sizes["normal"]   # general UI text
    _s = sizes["small"]    # tables, headers, tabs, status bar
    _l = sizes["large"]    # section headings
    _t = sizes["title"]    # not used in stylesheet directly but available

    c = THEMES.get(theme, THEMES[DEFAULT_THEME])
    BG_DARK, BG_MID, BG_LIGHT, BG_HOVER = c["BG_DARK"], c["BG_MID"], c["BG_LIGHT"], c["BG_HOVER"]
    ACCENT, ACCENT_DIM, ACCENT_BRIGHT = c["ACCENT"], c["ACCENT_DIM"], c["ACCENT_BRIGHT"]
    SUCCESS, ERROR, WARNING, MUTED = c["SUCCESS"], c["ERROR"], c["WARNING"], c["MUTED"]
    TEXT_PRIMARY, TEXT_SECONDARY = c["TEXT_PRIMARY"], c["TEXT_SECONDARY"]
    BORDER = c["BORDER"]
    SELECTION_BG = c["SELECTION_BG"]
    SUCCESS_BTN_BG, SUCCESS_BTN_HOVER = c["SUCCESS_BTN_BG"], c["SUCCESS_BTN_HOVER"]
    PRIMARY_PRESSED = c["PRIMARY_PRESSED"]

    return f"""
/* ── Global ── */
QMainWindow, QWidget {{
    background-color: {BG_DARK};
    color: {TEXT_PRIMARY};
    font-family: "{FONT_UI}", "Segoe UI", system-ui, sans-serif;
    font-size: {_n}pt;
}}

/* ── Sidebar ── */
QListWidget#sidebar {{
    background-color: {BG_MID};
    border: none;
    border-right: 1px solid {BORDER};
    outline: none;
    padding: 4px 0;
}}
QListWidget#sidebar::item {{
    padding: 10px 16px 10px 20px;
    border-left: 3px solid transparent;
    color: {TEXT_SECONDARY};
    font-family: "{FONT_UI}", sans-serif;
    font-size: {_n}pt;
}}
QListWidget#sidebar::item:selected {{
    background-color: {BG_LIGHT};
    border-left: 3px solid {ACCENT};
    color: {TEXT_PRIMARY};
}}
QListWidget#sidebar::item:hover:!selected {{
    background-color: {BG_HOVER};
    color: {TEXT_PRIMARY};
    border-left: 3px solid {MUTED};
}}

/* ── Buttons ── */
QPushButton {{
    background-color: {BG_LIGHT};
    color: {TEXT_PRIMARY};
    border: 1px solid {BORDER};
    border-radius: {BUTTON_RADIUS}px;
    padding: 5px 16px;
    min-height: {BUTTON_HEIGHT}px;
    font-family: "{FONT_UI}", sans-serif;
    font-size: {_n}pt;
    font-weight: 500;
}}
QPushButton:hover {{
    background-color: {BG_HOVER};
    border-color: {ACCENT};
    color: {TEXT_PRIMARY};
}}
QPushButton:pressed {{
    background-color: {ACCENT_DIM};
    border-color: {ACCENT_DIM};
}}
QPushButton:disabled {{
    color: {MUTED};
    border-color: {BORDER};
    background-color: {BG_MID};
}}
QPushButton#primary {{
    background-color: {ACCENT_DIM};
    border-color: {ACCENT_DIM};
    color: #ffffff;
    font-weight: 600;
}}
QPushButton#primary:hover {{
    background-color: {ACCENT};
    border-color: {ACCENT};
}}
QPushButton#primary:pressed {{
    background-color: {PRIMARY_PRESSED};
    border-color: {PRIMARY_PRESSED};
}}
QPushButton#success {{
    background-color: {SUCCESS_BTN_BG};
    border-color: {SUCCESS};
    color: {SUCCESS};
    font-weight: 500;
}}
QPushButton#success:hover {{
    background-color: {SUCCESS_BTN_HOVER};
}}

/* ── Inputs ── */
QLineEdit, QTextEdit {{
    background-color: {BG_LIGHT};
    border: 1px solid {BORDER};
    border-radius: {BUTTON_RADIUS}px;
    padding: 5px 10px;
    color: {TEXT_PRIMARY};
    font-family: "{FONT_UI}", sans-serif;
    selection-background-color: {ACCENT_DIM};
    min-height: 28px;
}}
QLineEdit:focus, QTextEdit:focus {{
    border-color: {ACCENT};
    outline: none;
}}
QSpinBox, QDoubleSpinBox {{
    background-color: {BG_LIGHT};
    border: 1px solid {BORDER};
    border-radius: {BUTTON_RADIUS}px;
    padding: 4px 8px;
    padding-right: 20px;
    color: {TEXT_PRIMARY};
    font-family: "{FONT_UI}", sans-serif;
    min-height: 28px;
    min-width: 80px;
}}
QSpinBox:focus, QDoubleSpinBox:focus {{
    border-color: {ACCENT};
}}
QSpinBox::up-button, QDoubleSpinBox::up-button {{
    subcontrol-origin: border;
    subcontrol-position: top right;
    width: 18px;
    border-left: 1px solid {BORDER};
    background: {BG_HOVER};
    border-top-right-radius: {BUTTON_RADIUS}px;
}}
QSpinBox::down-button, QDoubleSpinBox::down-button {{
    subcontrol-origin: border;
    subcontrol-position: bottom right;
    width: 18px;
    border-left: 1px solid {BORDER};
    background: {BG_HOVER};
    border-bottom-right-radius: {BUTTON_RADIUS}px;
}}
QSpinBox::up-arrow, QDoubleSpinBox::up-arrow {{
    width: 0; height: 0;
    border-left: 4px solid transparent;
    border-right: 4px solid transparent;
    border-bottom: 5px solid {TEXT_SECONDARY};
}}
QSpinBox::down-arrow, QDoubleSpinBox::down-arrow {{
    width: 0; height: 0;
    border-left: 4px solid transparent;
    border-right: 4px solid transparent;
    border-top: 5px solid {TEXT_SECONDARY};
}}
QSpinBox::up-button:hover, QDoubleSpinBox::up-button:hover,
QSpinBox::down-button:hover, QDoubleSpinBox::down-button:hover {{
    background: {ACCENT_DIM};
}}
QSpinBox:disabled, QDoubleSpinBox:disabled {{
    color: {MUTED};
    border-color: {BORDER};
    background-color: {BG_MID};
}}
QComboBox:disabled {{
    color: {MUTED};
    border-color: {BORDER};
    background-color: {BG_MID};
}}
QLabel:disabled {{
    color: {MUTED};
}}
QComboBox {{
    background-color: {BG_LIGHT};
    border: 1px solid {BORDER};
    border-radius: {BUTTON_RADIUS}px;
    padding: 5px 10px;
    color: {TEXT_PRIMARY};
    font-family: "{FONT_UI}", sans-serif;
    min-height: 28px;
    selection-background-color: {ACCENT_DIM};
}}
QComboBox:focus {{
    border-color: {ACCENT};
}}
QComboBox::drop-down {{
    border: none;
    width: 24px;
}}
QComboBox QAbstractItemView {{
    background-color: {BG_LIGHT};
    border: 1px solid {BORDER};
    border-radius: 0px;
    selection-background-color: {ACCENT_DIM};
    color: {TEXT_PRIMARY};
    padding: 2px;
}}

/* ── Labels ── */
QLabel {{
    font-family: "{FONT_UI}", sans-serif;
    color: {TEXT_PRIMARY};
    background: transparent;
}}

/* ── Group box ── */
QGroupBox {{
    border: 1px solid {BORDER};
    border-radius: {CARD_RADIUS}px;
    margin-top: 14px;
    padding: 12px 16px 16px 16px;
    font-family: "{FONT_UI}", sans-serif;
    font-size: {_s}pt;
    font-weight: 600;
    color: {TEXT_SECONDARY};
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 12px;
    padding: 0 6px;
    color: {TEXT_SECONDARY};
    background-color: {BG_DARK};
    font-size: {_s}pt;
    font-weight: 600;
    letter-spacing: 0.5px;
    text-transform: uppercase;
}}

/* ── Tables ── */
QTableWidget, QTableView {{
    background-color: {BG_MID};
    alternate-background-color: {BG_DARK};
    border: 1px solid {BORDER};
    border-radius: {CARD_RADIUS}px;
    gridline-color: {BORDER};
    color: {TEXT_PRIMARY};
    font-family: "{FONT_UI}", sans-serif;
    font-size: {_s}pt;
    selection-background-color: {SELECTION_BG};
    selection-color: {TEXT_PRIMARY};
}}
QTableWidget::item, QTableView::item {{
    padding: 5px 10px;
    border: none;
}}
QTableWidget::item:selected, QTableView::item:selected {{
    background-color: {SELECTION_BG};
    color: {TEXT_PRIMARY};
}}
QHeaderView {{
    background-color: {BG_MID};
}}
QHeaderView::section {{
    background-color: {BG_LIGHT};
    color: {TEXT_SECONDARY};
    border: none;
    border-right: 1px solid {BORDER};
    border-bottom: 1px solid {BORDER};
    padding: 6px 10px;
    font-family: "{FONT_UI}", sans-serif;
    font-weight: 600;
    font-size: {_s}pt;
    letter-spacing: 0.3px;
}}
QHeaderView::section:first {{
    border-top-left-radius: {CARD_RADIUS}px;
}}

/* ── Progress bar ── */
QProgressBar {{
    background-color: {BG_LIGHT};
    border: none;
    border-radius: 3px;
    height: 6px;
    text-align: center;
    color: transparent;
}}
QProgressBar::chunk {{
    background-color: {ACCENT};
    border-radius: 3px;
}}

/* ── Scrollbars ── */
QScrollBar:vertical {{
    background: {BG_MID};
    width: 8px;
    border: none;
    margin: 0;
}}
QScrollBar::handle:vertical {{
    background: {BG_HOVER};
    border-radius: 4px;
    min-height: 24px;
    margin: 2px;
}}
QScrollBar::handle:vertical:hover {{
    background: {MUTED};
}}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{ height: 0; }}
QScrollBar:horizontal {{
    background: {BG_MID};
    height: 8px;
    border: none;
    margin: 0;
}}
QScrollBar::handle:horizontal {{
    background: {BG_HOVER};
    border-radius: 4px;
    min-width: 24px;
    margin: 2px;
}}
QScrollBar::handle:horizontal:hover {{
    background: {MUTED};
}}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{ width: 0; }}

/* ── Tab widget ── */
QTabWidget::pane {{
    border: 1px solid {BORDER};
    border-radius: {CARD_RADIUS}px;
    background: {BG_MID};
    top: -1px;
}}
QTabBar::tab {{
    background: {BG_DARK};
    color: {TEXT_SECONDARY};
    padding: 7px 18px;
    border: 1px solid {BORDER};
    border-bottom: none;
    border-top-left-radius: {BUTTON_RADIUS}px;
    border-top-right-radius: {BUTTON_RADIUS}px;
    margin-right: 2px;
    font-family: "{FONT_UI}", sans-serif;
    font-size: {_s}pt;
    font-weight: 500;
}}
QTabBar::tab:selected {{
    background: {BG_MID};
    color: {TEXT_PRIMARY};
    border-bottom: 2px solid {ACCENT};
}}
QTabBar::tab:hover:!selected {{
    background: {BG_LIGHT};
    color: {TEXT_PRIMARY};
}}

/* ── Checkbox ── */
QCheckBox {{
    color: {TEXT_PRIMARY};
    spacing: 8px;
    font-family: "{FONT_UI}", sans-serif;
    font-size: {_n}pt;
}}
QCheckBox:checked {{
    color: {ACCENT_BRIGHT};
    font-weight: 600;
}}
QCheckBox::indicator {{
    width: 18px;
    height: 18px;
    border: 2px solid {MUTED};
    border-radius: 4px;
    background: {BG_LIGHT};
}}
QCheckBox::indicator:checked {{
    background-color: {ACCENT};
    border-color: {ACCENT_BRIGHT};
    image: none;
}}
QCheckBox::indicator:hover {{
    border-color: {ACCENT};
}}

/* ── Status bar ── */
QStatusBar {{
    background-color: {BG_MID};
    border-top: 1px solid {BORDER};
    color: {TEXT_SECONDARY};
    font-family: "{FONT_MONO}", monospace;
    font-size: {_s}pt;
    padding: 3px 10px;
}}

/* ── Tooltip ── */
QToolTip {{
    background-color: {BG_LIGHT};
    color: {TEXT_PRIMARY};
    border: 1px solid {BORDER};
    padding: 6px 10px;
    border-radius: {BUTTON_RADIUS}px;
    font-family: "{FONT_UI}", sans-serif;
    font-size: {_s}pt;
}}

/* ── Scroll area ── */
QScrollArea {{
    border: none;
    background: transparent;
}}
QScrollArea > QWidget > QWidget {{
    background: transparent;
}}

/* ── Message box ── */
QMessageBox {{
    background-color: {BG_MID};
    color: {TEXT_PRIMARY};
}}
QMessageBox QPushButton {{
    min-width: 80px;
}}

/* ── File dialog ── */
QFileDialog {{
    background-color: {BG_MID};
    color: {TEXT_PRIMARY};
}}

/* ── Semantic text-colour roles ──
   Used via widget.setObjectName("role_X") instead of baking a colour into
   an inline per-widget stylesheet — inline f-string colours are frozen at
   construction time and never see a later theme change, since Qt gives a
   widget's own stylesheet priority over the global one. These roles use
   the exact same object-name QSS mechanism already relied on elsewhere
   (e.g. QPushButton#primary), so they correctly re-theme immediately
   whenever the global stylesheet is rebuilt and reapplied. */
QLabel#role_title {{
    color: {ACCENT};
}}
QLabel#role_secondary {{
    color: {TEXT_SECONDARY};
}}
QLabel#role_primary_text {{
    color: {TEXT_PRIMARY};
}}
QLabel#role_success {{
    color: {SUCCESS};
}}
QLabel#role_warning {{
    color: {WARNING};
}}
QLabel#role_error {{
    color: {ERROR};
}}
QLabel#role_muted {{
    color: {MUTED};
}}
QPushButton#role_error_btn {{
    color: {ERROR};
    border-color: {ERROR};
}}

/* ── Stat cards and banners ──
   Same rationale as the sidebar sections below: these use background/
   border combinations that a color:-only sweep missed, so they stayed
   frozen at construction time despite looking like they should react
   to the theme toggle. */
QWidget#stat_card {{
    background: {BG_MID};
    border: 1px solid {BORDER};
    border-radius: 6px;
    padding: 10px;
}}
QLabel#role_muted_bg {{
    background: {TEXT_SECONDARY};
    border-radius: 3px;
}}
QLabel#warning_banner {{
    background: {WARNING}22;
    color: {WARNING};
    padding: 8px;
    border-radius: 4px;
}}

/* ── Sidebar structural sections ──
   Same rationale as the text-colour roles above: these were previously
   set via inline f-string stylesheets baked in at construction time,
   which is why the sidebar stayed dark-themed even after the toggle. */
QWidget#sidebar_container {{
    background-color: {BG_MID};
    border-right: 1px solid {BORDER};
}}
QWidget#sidebar_header {{
    background-color: {BG_MID};
    border-bottom: 1px solid {BORDER};
    padding: 0;
}}
QWidget#sidebar_footer_section {{
    background: transparent;
    border-top: 1px solid {BORDER};
    padding: 8px 0 4px 0;
}}
QLabel#role_section_label {{
    color: {MUTED};
    font-size: 7pt;
    font-weight: 600;
    letter-spacing: 0.5px;
    border: none;
}}
QRadioButton#sidebar_radio {{
    color: {TEXT_SECONDARY};
    font-size: 8pt;
    spacing: 4px;
    border: none;
    background: transparent;
}}
QRadioButton#sidebar_radio:checked {{
    color: {ACCENT};
}}
QRadioButton#sidebar_radio::indicator {{
    width: 10px;
    height: 10px;
    border: 1px solid {MUTED};
    border-radius: 5px;
    background: transparent;
}}
QRadioButton#sidebar_radio::indicator:checked {{
    background: {ACCENT};
    border-color: {ACCENT};
}}

/* ── Splitter ── */
QSplitter::handle {{
    background: {BORDER};
}}
QSplitter::handle:hover {{
    background: {ACCENT_DIM};
}}
"""


# Default stylesheet at import time — used by main.py initial setStyleSheet call
STYLESHEET = build_stylesheet(FONT_SIZE_PRESET_DEFAULT, DEFAULT_THEME)