# GRACE
**G**enomic **R**epeat **A**nalysis and **C**haracterisation **E**ngine

A desktop application for genome-wide microsatellite (SSR) discovery, primer design, BLAST-based specificity screening, and multiplex genotyping panel design — from a raw reference genome to a laboratory-ready panel, in one pipeline.

---

## What GRACE does

GRACE takes a reference genome through five stages, each usable independently or as a full pipeline:

1. **Genome loading** — FASTA, FASTQ, or GenBank flat file input, normalised internally to a single representation regardless of source format.
2. **SSR detection** — perfect, imperfect (interrupted), and compound microsatellite detection, with genomic-feature annotation (exon/intron/intergenic) when a GFF is supplied.
3. **Primer design** — Primer3-backed primer design per locus, with multiple ranked candidate pairs per SSR and support for both capillary electrophoresis and amplicon-sequencing workflows.
4. **BLAST specificity screening** — genome-wide off-target checking using NCBI BLAST+, including 3′-end-aware mismatch evaluation.
5. **Multiplex panel design** — automatic assembly of dye-balanced, annealing-temperature-compatible multiplex PCR panels from validated primer pairs, with even pool-size distribution and configurable per-multiplex constraints.

At any stage, GRACE can export a **PDF summary report** (genome stats, SSR distribution charts, primer/Tm histograms, specificity outcomes, panel composition, and full primer sequence tables) — a report can be generated after any single stage, not only after a complete pipeline run.

## Features

- **Perfect, imperfect, and compound SSR detection**
- **Memory- and CPU-aware scanning** — automatically scales worker processes to available system memory, so large genomes degrade gracefully on constrained hardware instead of crashing
- **Coarse-to-fine scanning** using vectorised pre-filtering to skip regions with no repeat potential before running the exact detection pass
- **GFF/GTF genomic feature annotation** — flag which SSRs fall in exons, introns, or intergenic regions
- **Adaptive flanking** for primer design — automatically extends the search region per locus to make requested product sizes genuinely reachable
- **BLAST specificity screening** using short-sequence BLAST settings (word_size=7, configurable 3′-end mismatch threshold) with BLAST's native default e-value of 10, and memory-aware BLAST threading that scales `-num_threads` to available memory relative to genome size
- **Multiplex panel design** — greedy constructive optimiser balancing dye assignment, annealing-temperature compatibility, product-size separation, and contig spread
- **Session save/resume** — save and reload a full working session (genome path, SSRs, primers, specificity results, panel design) as a JSON file
- **PDF report export** — available at any pipeline stage, with charts matching what's shown on-screen
- **Light/dark theme toggle**

## Requirements

- Python 3.10+
- [PyQt6](https://pypi.org/project/PyQt6/)
- [primer3-py](https://github.com/libnano/primer3-py)
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) (`blastn`, `makeblastdb` executables) — Windows binaries (BLAST+ 2.17.0) are bundled in `blast/windows/`, so Windows users don't need to install it. On macOS/Linux, place the binaries under `blast/macos/` or `blast/linux/`, or ensure they're on `PATH`
- `pandas`, `numpy`, `psutil`

Install dependencies:

```bash
pip install -r requirements.txt
```

## Running GRACE

```bash
python main.py
```

## Building a standalone executable

```bash
pip install pyinstaller
python build.py
```

Output is written to `dist/GRACE/` (plus `dist/GRACE.app` on macOS). Distribute the whole `dist/GRACE/` folder; the BLAST+ binaries for the build platform are bundled automatically.

## GUI Workflow

The application is organised as seven sequential steps in the sidebar:

| Step | Purpose |
|---|---|
| Home | Load a genome (and optionally a GFF/GTF annotation file); save/load sessions |
| SSR Detection | Configure motif lengths, minimum repeat counts, and compound/imperfect detection; run the scan |
| SSR Summary | Charts and statistics for detected SSRs |
| Primer Design | Configure Primer3 parameters and design primers for detected loci |
| Specificity Check | Run BLAST specificity screening against the loaded genome |
| BLAST Results | Review PASS/FAIL results, export primer lists |
| Panel Design | Assemble a dye-balanced multiplex panel from validated primers |

An **Export Report (PDF)** button in the sidebar is available at every step, generating a report covering whatever stages have been completed so far.

## Project Structure

```
GRACE/
├── main.py                     # Application entry point
├── app_state.py                # Shared application state (AppState)
├── build.py / GRACE.spec       # PyInstaller build
│
├── core/                       # Business logic — no GUI dependency
│   ├── fasta_loader.py          # FASTA/FASTQ/GenBank parsing
│   ├── gff_parser.py            # GFF/GTF annotation indexing
│   ├── ssr_detection.py         # SSR detection, extension, compound merging
│   ├── primer_design.py         # Primer3 integration
│   ├── primer_import.py         # Import existing primer pairs from file
│   ├── primer_specificity_blast.py  # BLAST invocation and specificity logic
│   ├── panel_design.py          # Multiplex panel optimiser
│   └── report_data.py           # Pure-Python data layer for the PDF report
│
├── ui/                          # PyQt6 presentation layer
│   ├── main_window.py             # Main window, sidebar, theme/session/report controls
│   ├── style.py                  # Colour palette, stylesheet, theme definitions
│   ├── charts.py                  # Shared chart-drawing (on-screen + PDF report)
│   ├── table_utils.py             # Shared table/label helpers
│   ├── report_generator.py        # PDF report generation
│   └── panels/                    # One panel per pipeline step
│       ├── home_panel.py
│       ├── ssr_panel.py
│       ├── ssr_summary_panel.py
│       ├── primer_panel.py
│       ├── specificity_panel.py
│       ├── results_panel.py
│       └── panel_design_panel.py
│
└── blast/                       # Bundled BLAST+ binaries (windows/; add macos/ or linux/ as needed)
```

## Known Limitations

- **In-memory genome architecture**: the whole genome is loaded into memory at once. Memory-aware worker scaling degrades gracefully on constrained hardware, but a genome larger than available RAM will still be slow even in the safest (sequential) mode. A disk-backed, memory-mapped architecture would remove this ceiling but represents a substantial rewrite, not yet undertaken.
- **Panel design is a greedy heuristic**, not a proven-optimal solver — a genuinely different assignment of loci to dyes/pools could in principle do marginally better in some cases.
- **Imperfect-repeat detection is heuristic**, not a formal alignment — in rare cases an unrelated adjacent repeat could be incorrectly absorbed as an interruption of a different locus.
- **BLAST specificity screening can be resource-intensive on large genomes.** BLAST's thread count is scaled to available memory relative to genome size, but BLAST's own search memory still grows with how permissive the search is (the default e-value of 10 favours sensitivity over speed), so large, repeat-dense genomes can take considerably longer and use substantial memory.

## Third-party software

GRACE relies on the following third-party software, which retains its own licence:

- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/) — redistributed in `blast/` as a United States Government work in the public domain (see the NCBI public domain notice).
- [Primer3](https://primer3.org/) via [primer3-py](https://github.com/libnano/primer3-py) — GPL-2.0 (primer3-py's bundled Primer3 library).
- [PyQt6](https://www.riverbankcomputing.com/software/pyqt/) — GPL-3.0 / commercial.
- pandas, NumPy, psutil — BSD-3-Clause.

## Licence

GRACE is released under the MIT Licence — see [LICENSE](LICENSE).
