"""
ssr_summary_panel.py — SSR Summary Panel
With fallback chromosome name assignment.
"""

from collections import Counter
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QScrollArea, QSizePolicy, QGridLayout,
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QPainter, QPen, QBrush, QColor

from ui.style import (
    ACCENT, WARNING, TEXT_SECONDARY, TEXT_PRIMARY,
    FONT_UI, FONT_MONO, FONT_SIZE_LARGE,
    FONT_SIZE_SMALL, PANEL_PADDING, BG_MID, BORDER,
    THEMES, DEFAULT_THEME,
)
from ui.charts import draw_pie_chart, draw_bar_chart


class PieChart(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._data = {}
        self.setMinimumHeight(200)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.set_theme(DEFAULT_THEME)

    def set_theme(self, theme_name):
        """
        QPainter draws pixels directly and never consults the Qt
        stylesheet system at all, so this widget can't pick up a theme
        change the way QSS-styled widgets do — it must be told explicitly.
        Stores the colours it needs as instance state, read in paintEvent().
        """
        c = THEMES.get(theme_name, THEMES[DEFAULT_THEME])
        self._text_color = c["TEXT_PRIMARY"]
        self._slice_border_color = c["BG_MID"]
        self.update()

    def set_data(self, data):
        self._data = data
        self.update()

    def paintEvent(self, event):
        if not self._data:
            return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        draw_pie_chart(
            painter, 0, 0, self.width(), self.height(), self._data,
            self._text_color, self._slice_border_color, FONT_UI, FONT_SIZE_SMALL,
        )


class VerticalBarChart(QWidget):
    """Vertical bar chart with X‑axis labels and Y‑axis counts."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self._data = []  # list of (label, count)
        self.setMinimumHeight(300)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.set_theme(DEFAULT_THEME)

    def set_theme(self, theme_name):
        c = THEMES.get(theme_name, THEMES[DEFAULT_THEME])
        self._text_color = c["TEXT_PRIMARY"]
        self._axis_color = c["TEXT_SECONDARY"]
        self.update()

    def set_data(self, data):
        self._data = data
        self.update()

    def paintEvent(self, event):
        if not self._data:
            return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        draw_bar_chart(
            painter, 0, 0, self.width(), self.height(), self._data,
            self._text_color, self._axis_color, FONT_UI, FONT_MONO,
        )


class SSRSummaryPanel(QWidget):
    def __init__(self, state, main_window):
        super().__init__()
        self.state = state
        self.mw = main_window
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

        title = QLabel("Genomic distribution and statistics of detected SSRs")
        title.setFont(QFont(FONT_UI, FONT_SIZE_LARGE + 2, QFont.Weight.Bold))
        title.setObjectName("role_title")
        L.addWidget(title)

        # Overview cards
        self._overview_grid = QGridLayout()
        L.addLayout(self._overview_grid)

        # Warning banner for fragmented genomes
        self._warning_banner = QLabel(
            "Genome is fragmented. Showing contigs in arbitrary order. "
            "Load a GFF with chromosome annotations for proper grouping."
        )
        self._warning_banner.setObjectName("warning_banner")
        self._warning_banner.setWordWrap(True)
        self._warning_banner.setVisible(False)
        L.addWidget(self._warning_banner)

        # Chromosome/Contig distribution histogram
        L.addWidget(QLabel("SSR Count by Chromosome / Contig (Top 15)"))
        self._chrom_bar = VerticalBarChart()
        L.addWidget(self._chrom_bar)

        # SSR Count by Sequence (vertical bars)
        L.addWidget(QLabel("SSR Count by Sequence — Top 15 of total sequences"))
        self._contig_bar = VerticalBarChart()
        L.addWidget(self._contig_bar)

        # Pie charts row (side‑by‑side)
        pie_row = QHBoxLayout()
        pie_row.setSpacing(20)

        # Motif type pie chart (left)
        motif_widget = QWidget()
        motif_layout = QVBoxLayout(motif_widget)
        motif_layout.addWidget(QLabel("Motif Type Distribution"))
        self._motif_pie = PieChart()
        motif_layout.addWidget(self._motif_pie)
        pie_row.addWidget(motif_widget)

        # Genomic feature pie chart (right)
        feature_widget = QWidget()
        feature_layout = QVBoxLayout(feature_widget)
        feature_layout.addWidget(QLabel("Genomic Feature Distribution"))
        self._feature_pie = PieChart()
        feature_layout.addWidget(self._feature_pie)
        pie_row.addWidget(feature_widget)

        L.addLayout(pie_row)

        # Store reference to feature widget to hide/show later
        self._feature_widget = feature_widget

        # Top motifs vertical bar chart
        L.addWidget(QLabel("Top 15 Canonical Motifs"))
        self._top_motifs_bar = VerticalBarChart()
        L.addWidget(self._top_motifs_bar)

    def _refresh(self):
        if not self.state.has_ssrs:
            return

        ssrs = self.state.ssrs
        genome = self.state.genome
        has_gff = self.state.gff_features is not None

        # --- FALLBACK: Ensure chrom_name is populated ---
        gff_index = self.state.gff_features
        if gff_index and hasattr(gff_index, 'chrom_names'):
            chrom_names = gff_index.chrom_names
            for ssr in ssrs:
                if "chrom_name" not in ssr:
                    contig = ssr.get("contig", "")
                    ssr["chrom_name"] = chrom_names.get(contig, contig)

        total_ssrs = len(ssrs)
        genome_size = sum(len(s) for s in genome.values()) if genome else 0
        density = (total_ssrs / (genome_size / 1_000_000)) if genome_size > 0 else 0

        motifs = [ssr.get("canonical_motif", ssr.get("motif", "")) for ssr in ssrs]
        motif_counts = Counter(motifs)
        unique_motifs = len(motif_counts)

        repeat_counts = [ssr.get("repeat_count", 0) for ssr in ssrs]
        mean_repeats = sum(repeat_counts) / total_ssrs if total_ssrs else 0

        # Overview cards
        for i in reversed(range(self._overview_grid.count())):
            w = self._overview_grid.itemAt(i).widget()
            if w:
                w.deleteLater()

        cards = [
            ("Total SSRs", f"{total_ssrs:,}"),
            ("SSR density", f"{density:.2f} / Mb"),
            ("Unique motifs", f"{unique_motifs:,}"),
            ("Mean repeat count", f"{mean_repeats:.1f}"),
        ]
        for col, (label, value) in enumerate(cards):
            card = QWidget()
            card.setObjectName("stat_card")
            cl = QVBoxLayout(card)
            vl = QLabel(value)
            vl.setFont(QFont(FONT_MONO, FONT_SIZE_LARGE, QFont.Weight.Bold))
            vl.setObjectName("role_title")
            cl.addWidget(vl)
            cl.addWidget(QLabel(label))
            self._overview_grid.addWidget(card, 0, col)

        # SSR Count by Sequence
        contig_counts = Counter(ssr.get("contig", "?") for ssr in ssrs)
        top_contigs = contig_counts.most_common(15)
        self._contig_bar.set_data(top_contigs)

        # Motif type distribution
        motif_types = {
            "Mononucleotide": 0, "Dinucleotide": 0, "Trinucleotide": 0,
            "Tetranucleotide": 0, "Pentanucleotide": 0, "Hexanucleotide": 0,
            "Compound": 0, "Imperfect": 0,
        }
        for ssr in ssrs:
            if ssr.get("is_compound"):
                motif_types["Compound"] += 1
                continue
            if ssr.get("is_imperfect"):
                motif_types["Imperfect"] += 1
                continue
            length = len(ssr.get("canonical_motif", ssr.get("motif", "")))
            if length == 1:
                motif_types["Mononucleotide"] += 1
            elif length == 2:
                motif_types["Dinucleotide"] += 1
            elif length == 3:
                motif_types["Trinucleotide"] += 1
            elif length == 4:
                motif_types["Tetranucleotide"] += 1
            elif length == 5:
                motif_types["Pentanucleotide"] += 1
            elif length >= 6:
                motif_types["Hexanucleotide"] += 1
        self._motif_pie.set_data(motif_types)

        # Remove any motif type that has zero count (not scanned)
        motif_types = {k: v for k, v in motif_types.items() if v > 0}
        self._motif_pie.set_data(motif_types)

        # Genomic feature distribution (only if GFF loaded)
        if has_gff:
            features = [ssr.get("genomic_feature", "intergenic") for ssr in ssrs]
            feature_counts = Counter(features)
            feature_data = {
                "Exon": feature_counts.get("exon", 0),
                "Intron": feature_counts.get("intron", 0),
                "Intergenic": feature_counts.get("intergenic", 0),
            }
            self._feature_pie.set_data(feature_data)
            self._feature_widget.setVisible(True)
        else:
            self._feature_widget.setVisible(False)

        # Top motifs
        top_motifs = motif_counts.most_common(15)
        self._top_motifs_bar.set_data(top_motifs)

        # Chromosome/Contig distribution histogram
        counts = {}
        for ssr in ssrs:
            display = ssr.get("chrom_name", ssr.get("contig", "?"))
            counts[display] = counts.get(display, 0) + 1

        has_chrom_names = any(ssr.get("chrom_name") for ssr in ssrs)
        many_contigs = len(counts) > 50
        self._warning_banner.setVisible(not has_chrom_names and many_contigs)

        top_items = sorted(counts.items(), key=lambda x: x[1], reverse=True)[:15]
        self._chrom_bar.set_data(top_items)

    def on_show(self):
        self._refresh()

    def set_theme(self, theme_name):
        """Propagate a theme change down to every custom-painted chart —
        QPainter widgets don't respond to the global QSS stylesheet at all,
        so each one must be told explicitly."""
        for chart in (self._chrom_bar, self._contig_bar, self._motif_pie,
                      self._feature_pie, self._top_motifs_bar):
            chart.set_theme(theme_name)