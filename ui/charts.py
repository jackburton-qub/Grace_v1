"""
charts.py
---------
Shared QPainter-based chart drawing, extracted from the on-screen
PieChart/VerticalBarChart widgets in ssr_summary_panel.py so the exact
same drawing logic can render either to a widget (at its own origin) or
directly onto a PDF page at an arbitrary position — guaranteeing the PDF
report's charts look identical to what's shown on-screen, rather than
maintaining two separate implementations that could silently drift apart.

Every function takes an explicit (x, y, w, h) region rather than assuming
a widget's own (0, 0) origin and size, which is the only real change from
the original widget code — the drawing logic itself is unchanged.
"""

from PyQt6.QtGui import QPainter, QPen, QBrush, QColor, QFont
from PyQt6.QtCore import Qt

CHART_COLORS = [
    "#4C9BE8", "#56C596", "#F4A261", "#E76F51",
    "#9B89C4", "#48CAE4", "#F6C90E", "#E07BA5",
]


def draw_table(painter, x: int, y: int, page_bottom: int, col_widths: list, headers: list,
                rows: list, header_bg: str, header_text: str, text_color: str,
                alt_row_bg: str, font_family: str = "Consolas", font_size: int = 8,
                row_height: int = 20):
    """
    Draw a table starting at (x, y), stopping before page_bottom. Returns
    the number of data rows actually drawn — if that's less than
    len(rows), the caller should start a new page and call this again
    with rows[returned_count:] to continue, redrawing the header there
    too (this function always draws its own header when called).

    col_widths: list of column widths in device units, same length as headers.
    rows: list of row-tuples, each the same length as headers.
    """
    total_w = sum(col_widths)

    # Header row
    header_h = row_height + 6
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(QBrush(QColor(header_bg)))
    painter.drawRect(x, y, total_w, header_h)
    painter.setFont(QFont(font_family, font_size, QFont.Weight.Bold))
    painter.setPen(QPen(QColor(header_text)))
    cx = x
    for header, w in zip(headers, col_widths):
        painter.drawText(cx + 6, y + header_h - 8, str(header))
        cx += w
    y += header_h

    painter.setFont(QFont(font_family, font_size))
    n_drawn = 0
    for i, row_data in enumerate(rows):
        if y + row_height > page_bottom:
            break
        if i % 2 == 1:
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QBrush(QColor(alt_row_bg)))
            painter.drawRect(x, y, total_w, row_height)
        painter.setPen(QPen(QColor(text_color)))
        cx = x
        for val, w in zip(row_data, col_widths):
            painter.drawText(cx + 6, y + row_height - 6, str(val))
            cx += w
        y += row_height
        n_drawn += 1

    return n_drawn


def draw_pie_chart(painter: QPainter, x: int, y: int, w: int, h: int, data: dict,
                    text_color: str, border_color: str, font_family: str = "Arial", font_size: int = 9):
    """Draw a pie chart with a legend into the region (x, y, w, h). `data` is {label: count}."""
    if not data:
        return
    total = sum(data.values())
    if total == 0:
        return
    size = min(w, h) - 40
    if size <= 0:
        return
    cx, cy = x + w // 2, y + h // 2
    start_angle = 0
    legend_x = cx + size // 2 + 20
    legend_y = cy - (len(data) * 20) // 2
    for i, (label, count) in enumerate(data.items()):
        angle = int(360 * count / total)
        color = QColor(CHART_COLORS[i % len(CHART_COLORS)])
        painter.setBrush(QBrush(color))
        painter.setPen(QPen(QColor(border_color), 2))
        painter.drawPie(cx - size // 2, cy - size // 2, size, size, start_angle * 16, angle * 16)
        painter.setPen(QPen(QColor(text_color), 1))
        painter.setFont(QFont(font_family, font_size))
        pct = (count / total) * 100
        painter.drawText(legend_x, legend_y + i * 20, f"{label}: {count:,} ({pct:.1f}%)")
        start_angle += angle


def draw_bar_chart(painter: QPainter, x: int, y: int, w: int, h: int, data: list,
                    text_color: str, axis_color: str,
                    label_font_family: str = "Arial", axis_font_family: str = "Consolas"):
    """Draw a vertical bar chart with axis labels into the region (x, y, w, h).
    `data` is a list of (label, count) tuples."""
    if not data:
        return
    margin_left, margin_right, margin_top, margin_bottom = 50, 30, 30, 50
    chart_w = w - margin_left - margin_right
    chart_h = h - margin_top - margin_bottom
    if chart_w <= 0 or chart_h <= 0:
        return

    max_count = max(c for _, c in data) if data else 1
    max_count = max(max_count, 1)
    n = len(data)

    bar_width = max(10, chart_w // n - 4)

    axis_left   = x + margin_left
    axis_right  = x + w - margin_right
    axis_top    = y + margin_top
    axis_bottom = y + h - margin_bottom

    painter.setPen(QPen(QColor(axis_color), 1))
    painter.drawLine(axis_left, axis_bottom, axis_right, axis_bottom)
    painter.drawLine(axis_left, axis_top, axis_left, axis_bottom)

    painter.setFont(QFont(axis_font_family, 8))
    painter.setPen(QPen(QColor(text_color), 1))
    painter.drawText(axis_left - 25, axis_bottom + 4, "0")
    mid_y = axis_top + chart_h // 2
    painter.drawText(axis_left - 40, mid_y + 4, f"{max_count // 2}")
    painter.drawText(axis_left - 40, axis_top + 10, f"{max_count}")

    for i, (label, count) in enumerate(data):
        bar_x = axis_left + i * (bar_width + 4) + 4
        bar_h = int((count / max_count) * chart_h)
        bar_y = axis_bottom - bar_h
        color = QColor(CHART_COLORS[i % len(CHART_COLORS)])
        painter.setBrush(QBrush(color))
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawRect(bar_x, bar_y, bar_width, bar_h)

        painter.setPen(QColor(text_color))
        painter.setFont(QFont(label_font_family, 8))
        # Truncate to fit the ACTUAL space this bar's slot has, measured
        # with real font metrics — a fixed character-count guess (the
        # previous approach) doesn't account for how narrow bars get once
        # there are many of them; a 10-character label like a contig
        # accession can still overflow into its neighbour even though
        # it's under a generic "12 characters is probably fine" cutoff.
        slot_width = bar_width + 4
        display_label = _fit_label(painter, str(label), slot_width)
        painter.drawText(bar_x, axis_bottom + 18, display_label)


def _fit_label(painter, text, max_width):
    """Shrink `text` (adding an ellipsis) until it fits within max_width,
    measured using the painter's currently-set font. Falls back to just
    the ellipsis if nothing else fits."""
    if painter.fontMetrics().horizontalAdvance(text) <= max_width:
        return text
    for n in range(len(text) - 1, 0, -1):
        candidate = text[:n] + "\u2026"
        if painter.fontMetrics().horizontalAdvance(candidate) <= max_width:
            return candidate
    return "\u2026"