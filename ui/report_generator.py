"""
report_generator.py
--------------------
Generates a comprehensive PDF summary report of the work carried out in
the current session, using core.report_data for all the underlying
numbers (fully unit-tested independently of Qt) and ui.charts for the
same chart-drawing logic used by the on-screen SSR Summary widgets.

Available at any stage of the pipeline — each section is only included
if that stage's data actually exists in state, since GRACE is explicitly
meant to be useful even when a user hasn't run the full pipeline through
to the end.

Always renders in light-theme colours regardless of the app's current
on-screen theme: a PDF meant for printing or sharing should be readable
on paper and in an email client, independent of whatever dark/light
preference happens to be set in the GUI at the time.
"""

from PyQt6.QtGui import QPdfWriter, QPainter, QPageSize, QFont, QColor, QPen, QBrush
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QFileDialog, QMessageBox

from ui.style import FONT_UI, FONT_MONO, THEMES
from ui.charts import draw_pie_chart, draw_bar_chart, draw_table
from core.report_data import compute_full_report_data

import datetime

# Always light theme for the PDF — see module docstring.
_C = THEMES["light"]
BG          = _C["BG_DARK"]
CARD_BG     = _C["BG_MID"]
TEXT        = _C["TEXT_PRIMARY"]
TEXT_MUTED  = _C["TEXT_SECONDARY"]
ACCENT      = _C["ACCENT"]
BORDER      = _C["BORDER"]
SUCCESS     = _C["SUCCESS"]
ERROR       = _C["ERROR"]
WARNING     = _C["WARNING"]

RESOLUTION = 150  # DPI
MARGIN = 60        # device units
CONTENT_TOP = 90   # y-position where section content starts, below the header


def _italic_font(family, size):
    f = QFont(family, size)
    f.setItalic(True)
    return f


def _page_size(writer):
    return writer.width(), writer.height()


def _draw_header_footer(painter, writer):
    """Consistent, minimal running header — just the app identity and a
    rule. The section name itself is only ever shown once, as the large
    title drawn by _section_title(), rather than repeated here too."""
    w, h = _page_size(writer)
    painter.setPen(QPen(QColor(TEXT_MUTED)))
    painter.setFont(_italic_font(FONT_UI, 8))
    painter.drawText(MARGIN, 30, "GRACE \u2014 Genomic Repeat Analysis and Characterisation Engine")
    painter.setPen(QPen(QColor(BORDER)))
    painter.drawLine(MARGIN, 42, w - MARGIN, 42)


def _section_title(painter, x, y, text):
    painter.setPen(QPen(QColor(ACCENT)))
    painter.setFont(QFont(FONT_UI, 16, QFont.Weight.Bold))
    painter.drawText(x, y, text)
    # A short accent underline for a bit of visual polish, distinguishing
    # section titles from ordinary bold text at a glance.
    metrics = painter.fontMetrics()
    title_w = metrics.horizontalAdvance(text)
    painter.setPen(QPen(QColor(ACCENT), 2))
    painter.drawLine(x, y + 8, x + min(title_w, 220), y + 8)
    return y + 34


def _draw_stats_block(painter, x, y, stats, card_width=None):
    """
    Draw a set of label/value pairs as a lightly-shaded card, with the
    value column positioned according to the WIDEST label actually used
    — measured with Qt's own font metrics rather than a guessed fixed
    offset, which is what previously let a long label (e.g. "Markers
    achieved / requested") run straight into its own value with no gap
    at all.

    stats: list of (label, value, value_color_or_None) tuples.
    Returns the y position just below the drawn card.
    """
    label_font = QFont(FONT_UI, 10)
    value_font = QFont(FONT_MONO, 11, QFont.Weight.Bold)

    painter.setFont(label_font)
    metrics = painter.fontMetrics()
    max_label_w = max(metrics.horizontalAdvance(label) for label, _, _ in stats)
    value_x = x + max_label_w + 30

    painter.setFont(value_font)
    v_metrics = painter.fontMetrics()
    max_value_w = max(v_metrics.horizontalAdvance(str(v)) for _, v, _ in stats)

    card_w = card_width or (max_label_w + max_value_w + 60)
    card_h = len(stats) * 24 + 16
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(QBrush(QColor(CARD_BG)))
    painter.drawRect(x - 12, y - 16, card_w, card_h)

    cy = y
    for label, value, color in stats:
        painter.setFont(label_font)
        painter.setPen(QPen(QColor(TEXT_MUTED)))
        painter.drawText(x, cy, label)
        painter.setFont(value_font)
        painter.setPen(QPen(QColor(color or TEXT)))
        painter.drawText(value_x, cy, str(value))
        cy += 24

    return y - 16 + card_h + 20


def _draw_subheading(painter, x, y, text):
    painter.setFont(QFont(FONT_UI, 11, QFont.Weight.Bold))
    painter.setPen(QPen(QColor(TEXT)))
    painter.drawText(x, y, text)
    return y + 14


def _wrapped_text(painter, x, y, max_width, text, color, italic=True):
    """Simple word-wrap within max_width, used for the free-text
    unplaced_reason explanation. Returns the y position after the text."""
    painter.setFont(_italic_font(FONT_UI, 9) if italic else QFont(FONT_UI, 9))
    painter.setPen(QPen(QColor(color)))
    words = text.split()
    line = ""
    for word in words:
        test = f"{line} {word}".strip()
        if painter.fontMetrics().horizontalAdvance(test) > max_width:
            painter.drawText(x, y, line)
            y += 14
            line = word
        else:
            line = test
    if line:
        painter.drawText(x, y, line)
        y += 14
    return y


def generate_report(state, path: str) -> None:
    """
    Write a multi-page PDF summary report to `path`, covering every
    pipeline stage that has data in `state`. Raises on failure — the
    caller (the GUI button handler) is responsible for catching and
    showing a message box.
    """
    data = compute_full_report_data(state)

    writer = QPdfWriter(path)
    writer.setPageSize(QPageSize(QPageSize.PageSizeId.A4))
    writer.setResolution(RESOLUTION)
    writer.setTitle("GRACE Summary Report")

    painter = QPainter(writer)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    page_w, page_h = _page_size(writer)

    # ---------------- Title page ----------------
    painter.setPen(QPen(QColor(ACCENT)))
    painter.setFont(QFont(FONT_UI, 36, QFont.Weight.Bold))
    painter.drawText(MARGIN, 160, "GRACE")
    painter.setFont(_italic_font(FONT_UI, 14))
    painter.setPen(QPen(QColor(TEXT)))
    painter.drawText(MARGIN, 195, "Genomic Repeat Analysis and Characterisation Engine")
    painter.setFont(QFont(FONT_UI, 18, QFont.Weight.Bold))
    painter.drawText(MARGIN, 250, "Summary Report")

    painter.setFont(QFont(FONT_UI, 10))
    painter.setPen(QPen(QColor(TEXT_MUTED)))
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    painter.drawText(MARGIN, 280, f"Generated {timestamp}")
    if state.genome_filename:
        painter.drawText(MARGIN, 300, f"Genome: {state.genome_filename}")

    checklist_y = 350
    painter.setFont(QFont(FONT_UI, 12, QFont.Weight.Bold))
    painter.setPen(QPen(QColor(TEXT)))
    painter.drawText(MARGIN, checklist_y, "Stages included in this report")
    checklist_y += 30
    stage_labels = [
        ("genome", "Genome loaded"),
        ("ssrs", "SSR detection"),
        ("primers", "Primer design"),
        ("specificity", "BLAST specificity screening"),
        ("panel", "Multiplex panel design"),
    ]
    painter.setFont(QFont(FONT_UI, 11))
    for key, label in stage_labels:
        done = data.get(key) is not None
        mark = "\u2713" if done else "\u2014"
        painter.setPen(QPen(QColor(SUCCESS if done else TEXT_MUTED)))
        painter.drawText(MARGIN, checklist_y, f"{mark}  {label}")
        checklist_y += 24

    # ---------------- Genome section ----------------
    if data["genome"]:
        writer.newPage()
        _draw_header_footer(painter, writer)
        g = data["genome"]
        y = _section_title(painter, MARGIN, CONTENT_TOP, "Genome")
        y = _draw_stats_block(painter, MARGIN, y, [
            ("Contigs", f"{g['n_contigs']:,}", None),
            ("Total length (bp)", f"{g['total_bp']:,}", None),
            ("Largest contig (bp)", f"{g['largest_contig_bp']:,}", None),
        ])
        # A single-contig genome makes a "largest contigs" bar chart
        # pointless (one bar, nothing to compare) — only draw it when
        # there's genuinely more than one contig to compare.
        if len(g["top_contigs_bar"]) > 1:
            y = _draw_subheading(painter, MARGIN, y, "Largest contigs")
            draw_bar_chart(
                painter, MARGIN, y, page_w - 2 * MARGIN, 240,
                g["top_contigs_bar"], TEXT, TEXT_MUTED, FONT_UI, FONT_MONO,
            )

    # ---------------- SSR section ----------------
    if data["ssrs"]:
        writer.newPage()
        _draw_header_footer(painter, writer)
        s = data["ssrs"]
        y = _section_title(painter, MARGIN, CONTENT_TOP, "SSR Detection")
        y = _draw_stats_block(painter, MARGIN, y, [
            ("Total SSRs", f"{s['total_ssrs']:,}", None),
            ("SSR density", f"{s['density_per_mb']:.2f} / Mb", None),
            ("Unique motifs", f"{s['unique_motifs']:,}", None),
            ("Mean repeat count", f"{s['mean_repeats']:.1f}", None),
        ])

        half_w = (page_w - 2 * MARGIN - 20) // 2
        chart_top = y
        if s["motif_type_pie"]:
            _draw_subheading(painter, MARGIN, chart_top, "Motif type distribution")
            draw_pie_chart(
                painter, MARGIN, chart_top + 10, half_w, 210,
                s["motif_type_pie"], TEXT, BG, FONT_UI, 9,
            )
        if s["feature_pie"]:
            _draw_subheading(painter, MARGIN + half_w + 20, chart_top, "Genomic feature distribution")
            draw_pie_chart(
                painter, MARGIN + half_w + 20, chart_top + 10, half_w, 210,
                s["feature_pie"], TEXT, BG, FONT_UI, 9,
            )
        y = chart_top + 240

        if s["top_motifs_bar"]:
            y = _draw_subheading(painter, MARGIN, y, "Most common motifs")
            draw_bar_chart(
                painter, MARGIN, y, page_w - 2 * MARGIN, 210,
                s["top_motifs_bar"], TEXT, TEXT_MUTED, FONT_UI, FONT_MONO,
            )
            y += 230

        if s["chrom_bar"]:
            if y > page_h - 260:
                writer.newPage()
                _draw_header_footer(painter, writer)
                y = CONTENT_TOP
            y = _draw_subheading(painter, MARGIN, y, "SSRs per contig/chromosome")
            draw_bar_chart(
                painter, MARGIN, y, page_w - 2 * MARGIN, 230,
                s["chrom_bar"], TEXT, TEXT_MUTED, FONT_UI, FONT_MONO,
            )

    # ---------------- Primer design section ----------------
    if data["primers"]:
        writer.newPage()
        _draw_header_footer(painter, writer)
        p = data["primers"]
        y = _section_title(painter, MARGIN, CONTENT_TOP, "Primer Design")
        y = _draw_stats_block(painter, MARGIN, y, [
            ("Primer pairs designed", f"{p['n_designed']:,}", None),
            ("Unique SSR loci covered", f"{p['n_unique_ssrs']:,}", None),
        ])

        if p["product_size_hist"]:
            y = _draw_subheading(painter, MARGIN, y, "Product size distribution (bp)")
            draw_bar_chart(
                painter, MARGIN, y, page_w - 2 * MARGIN, 210,
                p["product_size_hist"], TEXT, TEXT_MUTED, FONT_UI, FONT_MONO,
            )
            y += 230

        if p["tm_hist"]:
            y = _draw_subheading(painter, MARGIN, y, "Melting temperature distribution (\u00b0C)")
            draw_bar_chart(
                painter, MARGIN, y, page_w - 2 * MARGIN, 210,
                p["tm_hist"], TEXT, TEXT_MUTED, FONT_UI, FONT_MONO,
            )

    # ---------------- Specificity section ----------------
    if data["specificity"]:
        writer.newPage()
        _draw_header_footer(painter, writer)
        sp = data["specificity"]
        y = _section_title(painter, MARGIN, CONTENT_TOP, "BLAST Specificity Screening")
        y = _draw_stats_block(painter, MARGIN, y, [
            ("PASS", f"{sp['n_pass']:,}", SUCCESS),
            ("FAIL", f"{sp['n_fail']:,}", ERROR),
            ("UNKNOWN", f"{sp['n_unknown']:,}", WARNING),
            ("Pass rate", f"{sp['pass_rate']:.1f}%", None),
            ("Unique SSRs passed", f"{sp['n_unique_pass']:,}", None),
        ])

        if sp["status_pie"]:
            y = _draw_subheading(painter, MARGIN, y, "Specificity outcome")
            draw_pie_chart(
                painter, MARGIN, y, page_w - 2 * MARGIN, 240,
                sp["status_pie"], TEXT, BG, FONT_UI, 10,
            )

    # ---------------- Panel design section ----------------
    if data["panel"]:
        writer.newPage()
        _draw_header_footer(painter, writer)
        pd = data["panel"]
        y = _section_title(painter, MARGIN, CONTENT_TOP, "Multiplex Panel Design")
        y = _draw_stats_block(painter, MARGIN, y, [
            ("Markers achieved / requested", f"{pd['achieved']} / {pd['requested']}", None),
            ("Multiplex pools", f"{pd['n_pools']}", None),
        ])
        if pd["unplaced_reason"]:
            y = _wrapped_text(painter, MARGIN, y, page_w - 2 * MARGIN, pd["unplaced_reason"], WARNING)
            y += 10

        half_w = (page_w - 2 * MARGIN - 20) // 2
        chart_top = y
        if pd["dye_bar"]:
            _draw_subheading(painter, MARGIN, chart_top, "Markers per dye")
            draw_bar_chart(
                painter, MARGIN, chart_top + 10, half_w, 210,
                pd["dye_bar"], TEXT, TEXT_MUTED, FONT_UI, FONT_MONO,
            )
        if pd["pool_bar"]:
            _draw_subheading(painter, MARGIN + half_w + 20, chart_top, "Markers per multiplex pool")
            draw_bar_chart(
                painter, MARGIN + half_w + 20, chart_top + 10, half_w, 210,
                pd["pool_bar"], TEXT, TEXT_MUTED, FONT_UI, FONT_MONO,
            )

        # ---- Primer sequence table (paginated) ----
        if pd["markers_table"]:
            writer.newPage()
            _draw_header_footer(painter, writer)
            y = _section_title(painter, MARGIN, CONTENT_TOP, "Panel Design \u2014 Primer Sequences")

            headers  = ["Pool", "Dye", "SSR ID", "Motif", "Size (bp)", "Forward primer (5'\u21923')", "Reverse primer (5'\u21923')"]
            col_w    = [45, 50, 65, 90, 65, 260, 260]
            rows = [
                (m["pool"], m["dye"], m["ssr_id"], m["motif"], m["product_size"],
                 m["forward_primer"], m["reverse_primer"])
                for m in pd["markers_table"]
            ]

            idx = 0
            table_y = y
            while idx < len(rows):
                remaining = rows[idx:]
                n_drawn = draw_table(
                    painter, MARGIN, table_y, page_h - MARGIN, col_w, headers, remaining,
                    header_bg=ACCENT, header_text="#ffffff", text_color=TEXT,
                    alt_row_bg=CARD_BG, font_family=FONT_MONO, font_size=7, row_height=18,
                )
                idx += n_drawn
                if n_drawn == 0:
                    # Not even one row fit (shouldn't normally happen given
                    # the margins used) — avoid an infinite loop.
                    break
                if idx < len(rows):
                    writer.newPage()
                    _draw_header_footer(painter, writer)
                    table_y = CONTENT_TOP

    painter.end()


def export_report_dialog(state, parent_widget) -> None:
    """
    Show a save-file dialog and generate the report, with an error message
    box on failure and a status confirmation on success — the shared
    handler every panel's "Export Report" action calls.
    """
    if not any([state.has_genome, state.has_ssrs, state.has_primers,
                state.has_specificity, state.has_panel_design]):
        QMessageBox.information(
            parent_widget, "Nothing to report",
            "Load a genome and complete at least one pipeline stage before generating a report."
        )
        return

    path, _ = QFileDialog.getSaveFileName(
        parent_widget, "Export Summary Report", "GRACE_report.pdf", "PDF files (*.pdf)"
    )
    if not path:
        return

    try:
        generate_report(state, path)
    except Exception as e:
        QMessageBox.critical(parent_widget, "Report Failed", f"Could not generate report:\n\n{e}")
        return

    QMessageBox.information(parent_widget, "Report Saved", f"Summary report saved to:\n{path}")