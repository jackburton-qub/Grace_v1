"""
table_utils.py
--------------
Shared table helpers used across result panels (SSR, Primer, BLAST Results).

Previously this class/constant/helper was copy-pasted independently into
ssr_panel.py, primer_panel.py, and results_panel.py (and _lbl into
primer_panel.py + specificity_panel.py). Consolidated here so there's a
single place to change display limits or sort behaviour.
"""

from PyQt6.QtWidgets import QTableWidgetItem, QLabel
from PyQt6.QtCore import Qt
from ui.style import TEXT_SECONDARY, FONT_SIZE_SMALL


# Cap table rows to prevent main-thread freeze on very large datasets.
# (Raw BLAST hits use their own lower RAW_TABLE_LIMIT in results_panel.py,
# since that table can be far larger than any other in the app.)
TABLE_DISPLAY_LIMIT = 10_000


class NumericTableWidgetItem(QTableWidgetItem):
    """QTableWidgetItem that sorts numerically instead of alphabetically."""

    def __init__(self, text, numeric_value=None):
        super().__init__(text)
        self.numeric_value = numeric_value

    def __lt__(self, other):
        if (
            self.numeric_value is not None
            and isinstance(other, NumericTableWidgetItem)
            and other.numeric_value is not None
        ):
            return bool(self.numeric_value < other.numeric_value)
        return super().__lt__(other)


def secondary_label(text, tip=None):
    """
    Small muted QLabel, used for form field captions.

    When `tip` is provided, a small "ⓘ" indicator is appended to the text so
    the presence of help is visible rather than something the user has to
    randomly discover by hovering. The tooltip fires anywhere over the label,
    not just over the icon glyph.

    Deliberately plain-text (not rich text/HTML) so the label correctly
    dims via QSS when .setEnabled(False) is called on it — HTML-embedded
    colors would silently override the stylesheet's :disabled rule.
    """
    display_text = f"{text} \u24d8" if tip else text  # ⓘ circled small letter i
    w = QLabel(display_text)
    w.setObjectName("role_secondary")
    w.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}pt;")
    if tip:
        w.setToolTip(tip)
        w.setCursor(Qt.CursorShape.WhatsThisCursor)
    return w