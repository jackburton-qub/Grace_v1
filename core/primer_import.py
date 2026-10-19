"""
primer_import.py
----------------
Reads existing primer pairs from CSV/TSV or paired FASTA so they can be
screened for specificity (and used in panel design) without running SSR
detection or primer design first.

Accepted formats:

- CSV / TSV with a header row containing name, forward, reverse and
  optionally expected_size (case-insensitive, any column order).
- FASTA with one entry per primer, paired by name: >name_F and >name_R.

Every imported pair becomes a record shaped like a designed primer pair,
so the existing BLAST screening, results page and panel design can use it
unchanged. Imported records are marked with source="imported" and an
ssr_id of "IMP-<name>", so they stay identifiable in every output,
including ones that only show the ID. The ssr_id never parses as an
integer, which keeps it distinct from designed loci when BLAST query
names are mapped back to (ssr_id, pair_rank).
"""

import csv
import io
import os
import re
from dataclasses import dataclass, field
from typing import List, Optional

import primer3


IMPORTED_SOURCE = "imported"
IMPORTED_ID_PREFIX = "IMP-"

MIN_PRIMER_LEN = 15
MAX_PRIMER_LEN = 60

FASTA_EXTENSIONS = {".fa", ".fasta", ".fas", ".fna", ".ffn"}
_VALID_BASES = re.compile(r"^[ACGT]+$")
# Characters kept in the ssr_id. BLAST query IDs are cut at whitespace and
# "|" separates the primer side, so anything else is replaced with "_".
_UNSAFE_ID_CHARS = re.compile(r"[^A-Za-z0-9._-]")

# Same Tm conditions primer design uses (core/primer_design.py presets:
# 50 mM monovalent, 1.5 mM divalent, 0.6 mM dNTP, 50 nM DNA, SantaLucia
# Tm and salt correction), so imported and designed pairs cluster by Ta
# on the same footing in panel design.
_TM_ARGS = {
    "mv_conc": 50.0, "dv_conc": 1.5, "dntp_conc": 0.6, "dna_conc": 50.0,
    "tm_method": "santalucia", "salt_corrections_method": "santalucia",
}


@dataclass
class ImportResult:
    pairs: List[dict] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)    # file is rejected if any
    warnings: List[str] = field(default_factory=list)  # informational only

    @property
    def ok(self) -> bool:
        return not self.errors


def _gc_percent(seq: str) -> float:
    return round(100.0 * sum(1 for b in seq if b in "GC") / len(seq), 2)


def _check_primer(seq: str, where: str, label: str, errors: list) -> Optional[str]:
    """Normalise and validate one primer sequence; append problems to errors."""
    clean = re.sub(r"\s+", "", seq or "").upper()
    if not clean:
        errors.append(f"{where}: {label} primer is empty.")
        return None
    bad = sorted(set(clean) - set("ACGT"))
    if bad:
        errors.append(
            f"{where}: {label} primer contains invalid character(s) "
            f"{', '.join(repr(c) for c in bad)} — only A, C, G and T are allowed."
        )
        return None
    if not MIN_PRIMER_LEN <= len(clean) <= MAX_PRIMER_LEN:
        errors.append(
            f"{where}: {label} primer is {len(clean)} bp — expected "
            f"{MIN_PRIMER_LEN}-{MAX_PRIMER_LEN} bp."
        )
        return None
    return clean


def make_imported_pair(name: str, forward: str, reverse: str,
                       expected_size: Optional[int] = None) -> dict:
    """Build a primer record for an already-validated imported pair."""
    return {
        "ssr_id": IMPORTED_ID_PREFIX + _UNSAFE_ID_CHARS.sub("_", name),
        "pair_rank": 0,
        "source": IMPORTED_SOURCE,
        "primer_name": name,
        "contig": None,           # filled from the on-target BLAST amplicon
        "motif": "",
        "canonical_motif": "",
        "left_primer": forward,
        "right_primer": reverse,
        "left_tm": round(primer3.calc_tm(forward, **_TM_ARGS), 2),
        "right_tm": round(primer3.calc_tm(reverse, **_TM_ARGS), 2),
        "left_gc": _gc_percent(forward),
        "right_gc": _gc_percent(reverse),
        "expected_size": expected_size,
        "product_size": expected_size,  # None until BLAST predicts it
        "product_size_source": "expected_size" if expected_size else None,
    }


def _finish(raw_pairs: list, errors: list, warnings: list) -> ImportResult:
    """
    raw_pairs: [(where, name, forward, reverse, expected_size)]. Validates
    sequences and name uniqueness, then builds records — only if the whole
    file is clean, so a partly-broken file never half-imports.
    """
    if not raw_pairs and not errors:
        errors.append("No primer pairs found in the file.")

    built = []
    ids_seen = {}
    for where, name, fwd, rev, size in raw_pairs:
        fwd_ok = _check_primer(fwd, where, "Forward", errors)
        rev_ok = _check_primer(rev, where, "Reverse", errors)
        ssr_id = IMPORTED_ID_PREFIX + _UNSAFE_ID_CHARS.sub("_", name)
        if ssr_id in ids_seen:
            prev_where, prev_name = ids_seen[ssr_id]
            same = "the same name" if prev_name == name else f"a name that matches '{prev_name}' once special characters are replaced"
            errors.append(f"{where}: '{name}' has {same} as {prev_where}. Pair names must be unique.")
            continue
        ids_seen[ssr_id] = (where, name)
        if fwd_ok and rev_ok:
            if fwd_ok == rev_ok:
                warnings.append(f"{where}: '{name}' has identical forward and reverse primers.")
            built.append(make_imported_pair(name, fwd_ok, rev_ok, size))

    if errors:
        return ImportResult(errors=errors, warnings=warnings)
    return ImportResult(pairs=built, warnings=warnings)


def parse_delimited(text: str, delimiter: Optional[str] = None) -> ImportResult:
    """Parse CSV/TSV text with columns name, forward, reverse[, expected_size]."""
    errors, warnings = [], []
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        return ImportResult(errors=["The file is empty."])

    if delimiter is None:
        header_line = lines[0]
        delimiter = "\t" if header_line.count("\t") > header_line.count(",") else ","

    reader = csv.reader(io.StringIO("\n".join(lines)), delimiter=delimiter)
    header = [h.strip().lower().lstrip("﻿") for h in next(reader)]
    missing = [c for c in ("name", "forward", "reverse") if c not in header]
    if missing:
        return ImportResult(errors=[
            f"Header row is missing required column(s): {', '.join(missing)}. "
            f"Expected columns: name, forward, reverse, and optionally expected_size. "
            f"Found: {', '.join(h for h in header if h) or '(nothing)'}."
        ])
    col = {c: header.index(c) for c in ("name", "forward", "reverse", "expected_size") if c in header}

    raw_pairs = []
    for line_no, row in enumerate(reader, start=2):
        where = f"Row {line_no}"
        cells = [c.strip() for c in row]

        def cell(key):
            i = col.get(key)
            return cells[i] if i is not None and i < len(cells) else ""

        name = cell("name")
        if not name:
            errors.append(f"{where}: name is empty.")
            continue
        size = None
        size_text = cell("expected_size")
        if size_text:
            try:
                size = int(float(size_text))
                if size <= 0 or size != float(size_text):
                    raise ValueError
            except ValueError:
                errors.append(f"{where}: expected_size '{size_text}' is not a positive whole number of bp.")
                continue
        raw_pairs.append((f"{where} ('{name}')", name, cell("forward"), cell("reverse"), size))

    return _finish(raw_pairs, errors, warnings)


_PAIR_SUFFIX = re.compile(r"^(.+)_([FR])$", re.IGNORECASE)


def parse_fasta(text: str) -> ImportResult:
    """Parse FASTA where each pair is two entries named <name>_F and <name>_R."""
    errors, warnings = [], []
    entries = []  # (line_no, header_id, sequence)
    current = None
    for line_no, line in enumerate(text.splitlines(), start=1):
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            header_id = line[1:].split()[0] if line[1:].split() else ""
            current = [line_no, header_id, []]
            entries.append(current)
        elif current is None:
            errors.append(f"Line {line_no}: sequence found before any '>' header line.")
        else:
            current[2].append(line)

    by_name = {}
    order = []
    for line_no, header_id, seq_lines in entries:
        where = f"Line {line_no} ('>{header_id}')"
        m = _PAIR_SUFFIX.match(header_id)
        if not m:
            errors.append(f"{where}: name must end in _F (forward) or _R (reverse) so the pair can be matched.")
            continue
        name, side = m.group(1), m.group(2).upper()
        slot = by_name.setdefault(name, {})
        if name not in order:
            order.append(name)
        if side in slot:
            errors.append(f"{where}: duplicate {name}_{side} (first seen on line {slot[side][0]}).")
            continue
        slot[side] = (line_no, "".join(seq_lines))

    raw_pairs = []
    for name in order:
        slot = by_name[name]
        if "F" not in slot or "R" not in slot:
            have, need = ("F", "R") if "F" in slot else ("R", "F")
            errors.append(f"Line {slot[have][0]}: '{name}_{have}' has no matching '{name}_{need}'.")
            continue
        where = f"Pair '{name}' (lines {slot['F'][0]}/{slot['R'][0]})"
        raw_pairs.append((where, name, slot["F"][1], slot["R"][1], None))

    return _finish(raw_pairs, errors, warnings)


def import_primer_file(path: str) -> ImportResult:
    """Read and validate a primer file, choosing the parser from its extension/content."""
    try:
        with open(path, encoding="utf-8-sig") as fh:
            text = fh.read()
    except UnicodeDecodeError:
        return ImportResult(errors=["The file isn't plain text (UTF-8). Save it as CSV, TSV or FASTA and try again."])
    except OSError as e:
        return ImportResult(errors=[f"Could not read the file: {e}"])

    ext = os.path.splitext(path)[1].lower()
    if ext in FASTA_EXTENSIONS or text.lstrip().startswith(">"):
        return parse_fasta(text)
    return parse_delimited(text, "\t" if ext in (".tsv", ".tab") else ("," if ext == ".csv" else None))


def apply_blast_amplicon_sizes(results: list) -> list:
    """
    After BLAST screening, give each imported pair the contig of its
    predicted on-target amplicon, and its size where no expected_size was
    supplied. Only done when BLAST found exactly one in-range amplicon
    (always the case for a PASS); otherwise the pair has no single
    predicted product and its size stays unknown. Designed pairs are
    returned untouched. Mutates and returns `results`.
    """
    for rec in results:
        if rec.get("source") != IMPORTED_SOURCE:
            continue
        amplicons = rec.get("amplicons") or []
        if len(amplicons) != 1:
            rec["blast_amplicon_size"] = None
            continue
        amp = amplicons[0]
        rec["blast_amplicon_size"] = amp["product_size"]
        rec["contig"] = amp["contig"]
        if not rec.get("expected_size"):
            rec["product_size"] = amp["product_size"]
            rec["product_size_source"] = "BLAST amplicon"
    return results
