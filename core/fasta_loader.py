"""
FASTA / FASTQ / GBFF Loader Module (backend-safe)
--------------------------------------------------
Provides clean functions for loading genome sequence files into a dictionary.

Supported formats
-----------------
- FASTA  (.fa, .fasta, .fna)   — standard multi-FASTA
- FASTQ  (.fastq, .fq)         — sequence only, quality discarded
- GBFF   (.gb, .gbff, .gbk)    — GenBank flat file; extracts both sequence
                                  AND annotation features in one pass

Returns
-------
All loaders return a dict:  {contig_id: sequence_string}

GBFF loader additionally returns a GFFIndex via load_gbff() when annotation
is requested — callers that want both use load_gbff() directly.

Performance notes
-----------------
- FASTA/FASTQ: streamed line-by-line, no full-file buffer
- GBFF: streamed record-by-record; feature table parsed lazily
- IUPAC filter uses str.translate() (C-level, ~10× faster than char loop)
"""

import io
import re
import string
from typing import Optional, Tuple


# ---------------------------------------------------------------------------
# IUPAC filter
# ---------------------------------------------------------------------------

_IUPAC = set("ACGTNRYSWKMBDHV-")
_DELETE_TABLE = str.maketrans("", "", "".join(
    c for c in (string.ascii_letters + string.digits + string.punctuation + " ")
    if c.upper() not in _IUPAC
))


def _clean(line: str) -> str:
    return line.upper().translate(_DELETE_TABLE)


# ---------------------------------------------------------------------------
# FASTA loader
# ---------------------------------------------------------------------------

def load_fasta(handle) -> dict:
    """
    Load a FASTA file from a file path or file-like object.
    Returns {sequence_id: sequence}.
    """
    if isinstance(handle, str):
        file_obj    = open(handle, "r", encoding="utf-8", errors="replace")
        close_after = True
    else:
        content = handle.read()
        if isinstance(content, bytes):
            content = content.decode("utf-8", errors="replace")
        file_obj    = io.StringIO(content)
        close_after = True

    sequences  = {}
    current_id = None
    chunks     = []

    try:
        for raw_line in file_obj:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = _clean("".join(chunks))
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if current_id is not None:
            sequences[current_id] = _clean("".join(chunks))
    finally:
        if close_after:
            file_obj.close()

    return sequences


# ---------------------------------------------------------------------------
# FASTQ loader
# ---------------------------------------------------------------------------

def load_fastq(handle) -> dict:
    """
    Load a FASTQ file from a file path or file-like object.
    Returns {sequence_id: sequence}. Quality scores are discarded.
    """
    if isinstance(handle, str):
        file_obj    = open(handle, "r", encoding="utf-8", errors="replace")
        close_after = True
    else:
        content = handle.read()
        if isinstance(content, bytes):
            content = content.decode("utf-8", errors="replace")
        file_obj    = io.StringIO(content)
        close_after = True

    sequences = {}

    try:
        lines = iter(file_obj)
        for header_line in lines:
            header_line = header_line.strip()
            if not header_line:
                continue
            if not header_line.startswith("@"):
                continue
            seq_id   = header_line[1:].split()[0]
            seq_line = next(lines, "").strip()
            next(lines, None)   # "+" separator
            next(lines, None)   # quality
            if seq_id and seq_line:
                sequences[seq_id] = _clean(seq_line)
    finally:
        if close_after:
            file_obj.close()

    return sequences


# ---------------------------------------------------------------------------
# GBFF loader
# ---------------------------------------------------------------------------

# Feature types we care about for annotation
_GBFF_FEATURES = {"gene", "exon", "CDS", "mRNA"}

# Regex to parse a /qualifier="value" line
_QUALIFIER_RE = re.compile(r'^\s*/(\w+)="?([^"]*)"?\s*$')
_LOCUS_RE     = re.compile(r'^LOCUS\s+(\S+)')
_ACCESSION_RE = re.compile(r'^ACCESSION\s+(\S+)')
_VERSION_RE   = re.compile(r'^VERSION\s+(\S+)')


def _parse_location(loc_str: str):
    """
    Parse a GenBank location string into (start, end, strand).
    Handles complement(), join(), and simple N..M forms.
    Returns (start_0based, end_1based, strand) or None on failure.
    """
    loc_str = loc_str.strip()
    strand  = -1 if "complement" in loc_str else 1
    # Strip all wrappers
    clean = re.sub(r'complement\(|join\(|\)', '', loc_str)
    # Take first range (for join, we take the outer extent)
    ranges = re.findall(r'[<>]?(\d+)\.\.[<>]?(\d+)', clean)
    if not ranges:
        # Single position
        m = re.search(r'(\d+)', clean)
        if m:
            pos = int(m.group(1))
            return (pos - 1, pos, strand)
        return None
    starts = [int(r[0]) for r in ranges]
    ends   = [int(r[1]) for r in ranges]
    return (min(starts) - 1, max(ends), strand)


def load_gbff(
    path: str,
    build_annotation: bool = True,
    progress_callback=None,
) -> Tuple[dict, Optional[object]]:
    """
    Parse a GenBank flat file (.gb / .gbff / .gbk).

    Extracts:
      - Sequence for every LOCUS record → genome dict
      - Feature annotations (gene/exon/CDS) → GFFIndex (if build_annotation=True)
      - Chromosome name mappings from /chromosome qualifiers

    Args:
        path:               path to .gbff file
        build_annotation:   if True, build and return a GFFIndex
        progress_callback:  optional callable(n_records_parsed)

    Returns:
        (genome_dict, gff_index_or_None)
        genome_dict: {accession: sequence}
        gff_index:   GFFIndex with chrom_names populated, or None
    """
    from core.gff_parser import GFFIndex

    genome    = {}          # {accession: sequence}
    gff_index = GFFIndex() if build_annotation else None

    # State per LOCUS record
    current_id    = None    # accession used as contig ID
    current_locus = None    # LOCUS name (fallback if no accession)
    in_features   = False
    in_origin     = False
    seq_chunks    = []

    # Feature parsing state
    current_feat_type = None
    current_feat_loc  = None
    current_qualifiers = {}
    feat_loc_buf      = ""   # multi-line location accumulator

    n_records = 0

    def _flush_feature():
        """Commit current feature to gff_index."""
        if not build_annotation or not current_feat_type:
            return
        if current_feat_type not in _GBFF_FEATURES:
            return
        if not current_feat_loc or not current_id:
            return
        parsed = _parse_location(current_feat_loc)
        if not parsed:
            return
        start, end, strand = parsed

        feat = current_feat_type
        # Normalise to our feature vocabulary
        if feat == "mRNA":
            feat = "exon"

        gff_index._add(current_id, start, end, feat)

        # Extract chromosome name from /chromosome qualifier
        chrom = current_qualifiers.get("chromosome", "")
        if chrom and current_id and chrom != current_id:
            gff_index.chrom_names[current_id] = chrom

    def _flush_record():
        """Commit current LOCUS record to genome dict."""
        nonlocal current_id, current_locus, in_features, in_origin
        nonlocal seq_chunks, current_feat_type, current_feat_loc
        nonlocal current_qualifiers, feat_loc_buf, n_records

        _flush_feature()  # commit last feature of record

        if current_id and seq_chunks:
            seq = "".join(seq_chunks)
            if seq:
                genome[current_id] = _clean(seq)

        # Reset state
        current_id          = None
        current_locus       = None
        in_features         = False
        in_origin           = False
        seq_chunks          = []
        current_feat_type   = None
        current_feat_loc    = None
        current_qualifiers  = {}
        feat_loc_buf        = ""
        n_records          += 1

        if progress_callback and n_records % 100 == 0:
            progress_callback(n_records)

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")

            # ── Record separator ──────────────────────────
            if line.startswith("//"):
                _flush_record()
                continue

            # ── LOCUS ─────────────────────────────────────
            if line.startswith("LOCUS"):
                m = _LOCUS_RE.match(line)
                if m:
                    current_locus = m.group(1)
                    current_id    = current_locus   # default; overridden by VERSION
                continue

            # ── ACCESSION / VERSION ───────────────────────
            if line.startswith("ACCESSION") and not in_features and not in_origin:
                m = _ACCESSION_RE.match(line)
                if m:
                    current_id = m.group(1)
                continue

            if line.startswith("VERSION") and not in_features and not in_origin:
                m = _VERSION_RE.match(line)
                if m:
                    current_id = m.group(1)   # e.g. CP007106.1 — most specific
                continue

            # ── FEATURES block ────────────────────────────
            if line.startswith("FEATURES"):
                in_features = True
                in_origin   = False
                continue

            # ── ORIGIN (sequence) block ───────────────────
            if line.startswith("ORIGIN"):
                if in_features:
                    _flush_feature()
                    current_feat_type  = None
                    current_feat_loc   = None
                    current_qualifiers = {}
                    feat_loc_buf       = ""
                in_features = False
                in_origin   = True
                continue

            # ── Feature table parsing ─────────────────────
            if in_features and build_annotation:
                # Feature key line: 5-space indent, then feature type
                if len(line) > 5 and line[5] != " ":
                    # New feature — flush previous
                    _flush_feature()
                    current_feat_type  = line[5:].split()[0]
                    feat_loc_buf       = line[5 + len(current_feat_type):].strip()
                    current_feat_loc   = None
                    current_qualifiers = {}

                elif line.startswith("                     /"):
                    # Qualifier line — first resolve any pending location
                    if feat_loc_buf and not current_feat_loc:
                        current_feat_loc = feat_loc_buf
                        feat_loc_buf     = ""
                    m = _QUALIFIER_RE.match(line)
                    if m:
                        current_qualifiers[m.group(1)] = m.group(2)

                elif line.startswith("                     ") and feat_loc_buf:
                    # Continuation of a multi-line location
                    feat_loc_buf += line.strip()

                # Finalise location once we have it
                if feat_loc_buf and ")" not in feat_loc_buf and ".." not in feat_loc_buf:
                    pass  # still accumulating
                elif feat_loc_buf and not current_feat_loc:
                    current_feat_loc = feat_loc_buf
                    feat_loc_buf     = ""

            # ── Sequence lines ────────────────────────────
            elif in_origin:
                # ORIGIN lines: "      60 acgt acgt..."
                # Strip leading digits and spaces, keep only nucleotides
                seq_part = re.sub(r'[\d\s]', '', line)
                if seq_part:
                    seq_chunks.append(seq_part)

    # Flush final record if file doesn't end with "//"
    if current_id or seq_chunks:
        _flush_record()

    if build_annotation and gff_index:
        gff_index.build()

    return genome, gff_index


# ---------------------------------------------------------------------------
# Auto-dispatch
# ---------------------------------------------------------------------------

def write_fasta(genome_dict: dict, path: str, line_width: int = 70) -> None:
    """
    Write a {contig: sequence} dict out to a standard multi-line FASTA file.

    Used to reconstruct a genuinely valid FASTA file from an already-parsed
    genome dict — needed anywhere downstream (e.g. BLAST) that requires
    real FASTA input, since the ORIGINAL file the user loaded might not be
    FASTA at all (GBFF and FASTQ are both supported load formats, but
    neither is FASTA — passing them straight through to a tool that
    expects FASTA fails, since their format looks completely different).
    """
    with open(path, "w") as f:
        for contig, seq in genome_dict.items():
            f.write(f">{contig}\n")
            for i in range(0, len(seq), line_width):
                f.write(seq[i:i + line_width] + "\n")


def load_sequence_file(path_or_handle) -> dict:
    """
    Auto-detect format by file extension and dispatch to the correct loader.
    Returns {sequence_id: sequence}.

    For GBFF files, returns the genome dict only.
    Use load_gbff() directly if you also want the annotation index.
    """
    if isinstance(path_or_handle, str):
        lower = path_or_handle.lower()
        if lower.endswith((".fastq", ".fq")):
            return load_fastq(path_or_handle)
        if lower.endswith((".gb", ".gbff", ".gbk")):
            genome, _ = load_gbff(path_or_handle, build_annotation=False)
            return genome
        return load_fasta(path_or_handle)

    # File-like object — peek at first character
    content = path_or_handle.read()
    if isinstance(content, bytes):
        content = content.decode("utf-8", errors="replace")
    first = content.lstrip()[:6].upper()
    buf   = io.StringIO(content)
    if first.startswith("@"):
        return load_fastq(buf)
    if first.startswith("LOCUS"):
        genome, _ = load_gbff(buf, build_annotation=False)
        return genome
    return load_fasta(buf)