"""
dsap_toolkit.py

A small toolkit that automates the repetitive parts of the WSSP/DSAP
practice-clone workflow:

  1. find_insert_boundaries()  - locate the cDNA insert start (after the
                                  GAATTC / ATTAC / G-string adapter) and
                                  the end of the poly-A tail
  2. clean_ns()                - replace ambiguous 'N' base calls with 'A'
                                  (per the "cDNAs end in poly-A" rule)
  3. find_best_orf()           - translate all 3 forward reading frames
                                  and return the longest ORF (start/end
                                  base numbers, protein sequence)
  4. classify_clone()          - apply the "DSAP commandments" to decide
                                  clone type (full ORF / back-end partial /
                                  front-end partial / middle partial /
                                  non-coding) and what UTRs should exist


  5. analyze_clone()           - run ALL of the above automatically, in
                                  the right order, on one raw sequence
  6. analyze_files()           - run analyze_clone() on every sequence in
                                  a list of files / folders (.ab1, FASTA,
                                  plain text) and optionally write a CSV

Command line usage:

  python dsap_toolkit.py                       # demo on a toy sequence
  python dsap_toolkit.py clone1.ab1            # one chromatogram
  python dsap_toolkit.py my_clones/            # every sequence in a folder
  python dsap_toolkit.py my_clones/ -o results.csv
  python dsap_toolkit.py -s "GAATTC...AAAA"    # paste a raw sequence
"""

from __future__ import annotations
import argparse
import csv
import re
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# 1. INSERT BOUNDARY DETECTION
# ---------------------------------------------------------------------------
# The plasmid always carries the same three landmarks right before the
# cDNA insert: GAATTC, then ATTAC, then a run of G's. The insert itself
# starts on the base immediately after that G-string, and (being an
# mRNA/cDNA) ends in a poly-A tail.

ADAPTER_PATTERN = re.compile(r"GAATTC.{0,30}?ATTAC(G+)", re.IGNORECASE)
POLYA_PATTERN = re.compile(r"(A{10,})")  # a run of 10+ A's = poly-A tail


@dataclass
class InsertBoundaries:
    insert_start_base: int      # 1-indexed base where insert begins
    insert_end_base: int        # 1-indexed base of last A in poly-A tail
    adapter_found: bool
    polya_found: bool


def find_insert_boundaries(raw_sequence: str) -> InsertBoundaries:
    """
    Scan a raw (uncropped) waveform sequence for the GAATTC/ATTAC/G-string
    adapter and the poly-A tail, and report the crop points DSAP asks for.

    Bases are reported 1-indexed to match how DSAP/FinchTV label them
    (e.g. "A62").
    """
    seq = raw_sequence.strip().upper()

    adapter_match = ADAPTER_PATTERN.search(seq)
    if adapter_match:
        insert_start = adapter_match.end() + 1  # 1-indexed, base AFTER the G-string
        adapter_found = True
    else:
        insert_start = 1
        adapter_found = False

    # Look for poly-A tail starting from the insert onward (avoids matching
    # any stray A-runs upstream of the actual insert). Use the LAST run of
    # 10+ A's: the real tail is at the 3' end, and a 3' UTR can contain an
    # A-rich stretch of its own that would otherwise cut the insert short.
    search_region = seq[insert_start - 1:]
    polya_runs = list(POLYA_PATTERN.finditer(search_region))
    polya_match = polya_runs[-1] if polya_runs else None
    if polya_match:
        insert_end = (insert_start - 1) + polya_match.end()
        polya_found = True
    else:
        insert_end = len(seq)
        polya_found = False

    return InsertBoundaries(
        insert_start_base=insert_start,
        insert_end_base=insert_end,
        adapter_found=adapter_found,
        polya_found=polya_found,
    )


# ---------------------------------------------------------------------------
# 2. N -> A CLEANUP
# ---------------------------------------------------------------------------

def clean_ns(sequence: str) -> tuple[str, int]:
    """
    Replace ambiguous 'N' calls with 'A', per the DSAP rule that cDNAs end
    in a poly-A tail so unreadable trailing bases are usually A's.
    Returns (cleaned_sequence, number_of_replacements).
    """
    n_count = sequence.upper().count("N")
    cleaned = re.sub("N", "A", sequence, flags=re.IGNORECASE)
    return cleaned, n_count


# ---------------------------------------------------------------------------
# 3. ORF FINDER (translate 3 forward frames, find the longest ORF)
# ---------------------------------------------------------------------------

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate_frame(sequence: str, frame: int) -> str:
    """Translate sequence starting at 0-indexed offset `frame` (0, 1, or 2)."""
    seq = sequence.upper().replace("U", "T")
    protein = []
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        protein.append(CODON_TABLE.get(codon, "X"))
    return "".join(protein)


@dataclass
class ORFResult:
    frame: int              # 1, 2, or 3 (matches DSAP's +1/+2/+3 labels)
    protein: str
    start_base: int         # 1-indexed DNA base where ORF starts
    end_base: int            # 1-indexed DNA base of the LAST base of the stop codon
    has_start_codon: bool   # True if the ORF begins with M (a real ATG)
    has_stop_codon: bool    # True if it ends in a stop codon


MIN_ORF_AA = 30  # shorter "ORFs" are what random non-coding DNA produces


def find_best_orf(sequence: str, min_aa: int = MIN_ORF_AA) -> Optional[ORFResult]:
    """
    Translate all 3 forward reading frames and return the single longest
    ORF found (start codon M ... stop codon *, or the longest M...end /
    start...* run if no full M-to-* ORF exists, matching the "back end
    partial" / "front end partial" cases DSAP describes).

    ORFs shorter than `min_aa` amino acids are ignored. Without this, any
    stray stop codon near the start of a non-coding clone would be called a
    back-end partial (a stop codon turns up roughly every 21 codons in
    random sequence).
    """
    seq = sequence.upper().strip()
    candidates: list[ORFResult] = []

    for frame in range(3):  # 0-indexed internally
        protein = translate_frame(seq, frame)

        # Find every run that starts at 'M' and ends at '*' (full ORF)
        for m in re.finditer(r"M[^*]*\*", protein):
            aa_start, aa_end = m.start(), m.end()  # aa_end exclusive, includes '*'
            dna_start = frame + aa_start * 3 + 1          # 1-indexed
            dna_end = frame + aa_end * 3                   # last base of stop codon
            candidates.append(ORFResult(
                frame=frame + 1, protein=m.group().rstrip("*"),
                start_base=dna_start, end_base=dna_end,
                has_start_codon=True, has_stop_codon=True,
            ))

        # Back-end partial: no M at the front, but ends in a stop codon.
        # Take from base 1 of this frame to the first stop codon.
        first_stop = protein.find("*")
        if first_stop != -1:
            dna_start = frame + 1
            dna_end = frame + (first_stop + 1) * 3
            candidates.append(ORFResult(
                frame=frame + 1, protein=protein[:first_stop],
                start_base=dna_start, end_base=dna_end,
                has_start_codon=False, has_stop_codon=True,
            ))

        # Front-end partial: starts at M but runs off the end with no stop.
        first_m = protein.find("M")
        if first_m != -1 and "*" not in protein[first_m:]:
            dna_start = frame + first_m * 3 + 1
            dna_end = frame + len(protein) * 3
            candidates.append(ORFResult(
                frame=frame + 1, protein=protein[first_m:],
                start_base=dna_start, end_base=dna_end,
                has_start_codon=True, has_stop_codon=False,
            ))

    candidates = [c for c in candidates if len(c.protein) >= min_aa]
    if not candidates:
        return None

    # A full ORF (both start and stop codon present) is always the
    # biologically correct call over a partial, even if a partial happens
    # to translate to a longer string in some other frame. Among ORFs of
    # the same "completeness", the longest protein wins -- matching the
    # "which reading frame is most likely to code for a protein" logic in
    # the DSAP Toolbox.
    def rank(c: ORFResult):
        is_full = c.has_start_codon and c.has_stop_codon
        return (is_full, len(c.protein))

    return max(candidates, key=rank)


# ---------------------------------------------------------------------------
# 4. CLONE CLASSIFIER  ("the DSAP commandments" as logic)
# ---------------------------------------------------------------------------

@dataclass
class ClassificationResult:
    clone_type: str
    has_5utr: bool
    has_3utr: bool
    notes: str


def classify_clone(orf: Optional[ORFResult], sequence_length: int) -> ClassificationResult:
    """
    Apply the DSAP commandments to an ORFResult and report clone type +
    which UTRs should exist. Mirrors:
      I.   stop codon is part of the ORF
      II.  5' UTR always starts at base 1 (when it exists)
      III. back-end partial -> no 5' UTR
      IV.  back-end partial -> ORF starts at base 1
      V.   non-coding -> entire sequence is 3' UTR
      VI.  front-end partial -> no 3' UTR
      VII. middle-ORF partial -> no 5' UTR or 3' UTR
    """
    if orf is None:
        return ClassificationResult(
            clone_type="non-coding",
            has_5utr=False, has_3utr=True,
            notes="No ORF found (BLASTx too weak). Entire sequence is 3' UTR. "
                  "Set ORF and 5' UTR to N/A in DSAP.",
        )

    if orf.has_start_codon and orf.has_stop_codon:
        return ClassificationResult(
            clone_type="full ORF",
            has_5utr=orf.start_base > 1,
            has_3utr=orf.end_base < sequence_length,
            notes=f"Full ORF: {orf.start_base}-{orf.end_base} (frame +{orf.frame}).",
        )

    if orf.has_stop_codon and not orf.has_start_codon:
        return ClassificationResult(
            clone_type="back-end partial",
            has_5utr=False, has_3utr=orf.end_base < sequence_length,
            notes=(f"No start codon (broken off). Per commandment IV, treat ORF as "
                   f"starting at base 1 and ending at base {orf.end_base}. "
                   f"5' UTR = N/A."),
        )

    if orf.has_start_codon and not orf.has_stop_codon:
        return ClassificationResult(
            clone_type="front-end partial",
            has_5utr=orf.start_base > 1, has_3utr=False,
            notes=(f"Has start codon at {orf.start_base} but runs off the end "
                   f"with no stop codon -> no 3' UTR."),
        )

    return ClassificationResult(
        clone_type="middle ORF partial",
        has_5utr=False, has_3utr=False,
        notes="No start and no stop codon found -> no 5' UTR or 3' UTR.",
    )


# ---------------------------------------------------------------------------
# 5. E-VALUE SANITY CHECK
# ---------------------------------------------------------------------------

def passes_evalue_threshold(evalue: float, threshold: float = 1e-6) -> bool:
    """DSAP guidance: an E-value of e-06 or below is suitable for most
    genetic studies; 0.0 is the best possible match."""
    return evalue <= threshold


# ---------------------------------------------------------------------------
# 5b. BLASTx PROTEIN IDENTIFICATION (online, via NCBI)
# ---------------------------------------------------------------------------
# Sends the insert to NCBI BLASTx (the same search DSAP runs) and returns
# the best hit's protein name, e.g. "patatin-like protein 2".
# Results are cached in dsap_blast_cache.json so re-running a clone is
# instant instead of waiting on NCBI again.

BLAST_CACHE = Path(__file__).with_name("dsap_blast_cache.json")
UNINFORMATIVE = re.compile(
    r"hypothetical|uncharacteri[sz]ed|unnamed|unknown|predicted protein|"
    r"putative protein$", re.IGNORECASE)


@dataclass
class BlastHit:
    protein_name: str
    organism: str
    accession: str
    evalue: float
    identity_pct: float
    frame: str
    top_hits: str           # short summary of the first few hits


def _parse_hit_def(hit_def: str) -> tuple[str, str]:
    """'patatin-like protein 2 [Arabidopsis thaliana] >ref|...' ->
    ('patatin-like protein 2', 'Arabidopsis thaliana')"""
    first = hit_def.split(" >")[0].strip()
    m = re.match(r"(.*?)\s*\[([^\[\]]+)\]\s*$", first)
    if m:
        name = re.sub(r"^(PREDICTED|LOW QUALITY PROTEIN):\s*", "", m.group(1),
                      flags=re.IGNORECASE)
        return name.strip(), m.group(2).strip()
    return first, ""


NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def _ncbi_blastx_xml(sequence: str, database: str, hitlist_size: int,
                     timeout_min: float):
    """Submit a BLASTx job to NCBI and poll it, printing progress.
    Unlike Bio.Blast.NCBIWWW.qblast (which polls forever, silently), this
    gives up after `timeout_min` minutes and tells you the RID so you can
    look the job up in a browser."""
    import io
    import time
    import urllib.parse
    import urllib.request

    def call(params: dict, timeout: float = 60) -> str:
        data = urllib.parse.urlencode(params).encode()
        req = urllib.request.Request(NCBI_BLAST_URL, data=data,
                                     headers={"User-Agent": "dsap_toolkit"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read().decode("utf-8", errors="replace")

    page = call({"CMD": "Put", "PROGRAM": "blastx", "DATABASE": database,
                 "QUERY": sequence, "HITLIST_SIZE": hitlist_size,
                 "FORMAT_TYPE": "XML", "TOOL": "dsap_toolkit"})
    rid = re.search(r"RID = (\S+)", page)
    rtoe = re.search(r"RTOE = (\d+)", page)
    if not rid:
        raise RuntimeError("NCBI did not accept the BLAST request "
                           "(no RID returned) - try again in a few minutes")
    rid = rid.group(1)
    est = int(rtoe.group(1)) if rtoe else 30
    print(f"    submitted, RID={rid} (NCBI estimate ~{est}s). "
          f"Check in a browser: {NCBI_BLAST_URL}?CMD=Get&RID={rid}",
          file=sys.stderr)

    start = time.time()
    time.sleep(min(max(est, 10), 60))
    while True:
        elapsed = time.time() - start
        try:
            info = call({"CMD": "Get", "RID": rid,
                         "FORMAT_OBJECT": "SearchInfo"})
        except OSError as e:  # transient network hiccup: keep waiting
            print(f"    poll error ({e}), retrying...", file=sys.stderr)
            info = "Status=WAITING"
        status = re.search(r"Status=(\w+)", info)
        status = status.group(1) if status else "UNKNOWN"
        if status == "READY":
            if "ThereAreHits=yes" not in info:
                return None
            xml = call({"CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML"},
                       timeout=120)
            return io.StringIO(xml)
        if status in ("FAILED", "UNKNOWN"):
            raise RuntimeError(f"NCBI reports BLAST job {rid} as {status}")
        if elapsed > timeout_min * 60:
            raise TimeoutError(
                f"BLAST still running after {timeout_min:g} min (RID {rid}). "
                f"NCBI is busy - try --db swissprot, a longer --blast-timeout, "
                f"or open the RID link above later")
        print(f"    still running... {elapsed / 60:.1f} min elapsed",
              file=sys.stderr)
        time.sleep(30)


def run_blastx(sequence: str, database: str = "nr",
               hitlist_size: int = 10,
               timeout_min: float = 15) -> Optional[BlastHit]:
    """BLASTx the sequence at NCBI. Returns None when there are no hits.
    Prefers the best hit with an informative name ("patatin-like protein 2")
    over "hypothetical protein"-style hits with a similar E-value."""
    import hashlib
    import json

    key = f"{database}:{hashlib.sha1(sequence.encode()).hexdigest()}"
    cache = {}
    if BLAST_CACHE.exists():
        try:
            cache = json.loads(BLAST_CACHE.read_text())
        except ValueError:
            cache = {}
    if key in cache:
        return BlastHit(**cache[key]) if cache[key] else None

    try:
        from Bio.Blast import NCBIXML
    except ImportError:
        raise RuntimeError("BLAST needs Biopython: pip install biopython")

    handle = _ncbi_blastx_xml(sequence, database, hitlist_size, timeout_min)
    record = NCBIXML.read(handle) if handle is not None else None

    hits = []
    for aln in (record.alignments if record else []):
        hsp = aln.hsps[0]
        name, organism = _parse_hit_def(aln.hit_def)
        hits.append(BlastHit(
            protein_name=name, organism=organism, accession=aln.accession,
            evalue=hsp.expect,
            identity_pct=round(100 * hsp.identities / hsp.align_length, 1),
            frame=f"{hsp.frame[0]:+d}" if hsp.frame else "",
            top_hits="",
        ))

    best = None
    if hits:
        best = hits[0]
        for h in hits:  # first informative hit that's still a good match
            if not UNINFORMATIVE.search(h.protein_name) and \
                    h.evalue <= max(best.evalue * 1e5, 1e-6):
                best = h
                break
        best.top_hits = " | ".join(
            f"{h.protein_name} [{h.organism}] {h.evalue:.0e}" for h in hits[:3])

    cache[key] = asdict(best) if best else None
    BLAST_CACHE.write_text(json.dumps(cache, indent=1))
    return best


# ---------------------------------------------------------------------------
# 6. FULL AUTOMATED PIPELINE (one clone)
# ---------------------------------------------------------------------------
# Order matters:
#   clean N's FIRST  -> an N inside the poly-A tail would otherwise break
#                       the A{10,} run and hide the tail
#   find boundaries  -> crop to the insert
#   find ORF         -> on the cropped insert only (DSAP numbers bases
#                       from the start of the cropped sequence)
#   classify         -> clone type + UTR ranges ready to type into DSAP

def _fmt_range(start: int, end: int) -> str:
    return f"{start}-{end}" if start <= end else "N/A"


@dataclass
class CloneReport:
    name: str
    raw_length: int
    n_replaced: int
    adapter_found: bool
    polya_found: bool
    insert_start_base: int      # in raw-sequence coordinates
    insert_end_base: int        # in raw-sequence coordinates
    insert_length: int
    clone_type: str
    frame: str                  # "+1" / "+2" / "+3" or "N/A"
    utr5: str                   # ranges in CROPPED-insert coordinates
    orf: str
    utr3: str
    protein: str
    protein_length: int
    blast_protein: str          # e.g. "patatin-like protein 2"
    blast_organism: str
    blast_accession: str
    blast_identity_pct: Optional[float]
    blast_frame: str
    blast_top_hits: str
    evalue: Optional[float]
    evalue_ok: Optional[bool]
    warnings: str
    notes: str
    insert_sequence: str


def analyze_clone(raw_sequence: str, name: str = "clone",
                  evalue: Optional[float] = None, blast: bool = False,
                  blast_db: str = "nr",
                  blast_timeout: float = 15,
                  min_orf_aa: int = MIN_ORF_AA) -> CloneReport:
    """Run every step of the toolkit on one raw sequence.
    Accepts GenBank-style text too: line numbers, spaces and lowercase
    are stripped out automatically."""
    raw = re.sub(r"[^A-Za-z]", "", raw_sequence).upper()

    cleaned, n_count = clean_ns(raw)
    bounds = find_insert_boundaries(cleaned)
    insert = cleaned[bounds.insert_start_base - 1: bounds.insert_end_base]
    L = len(insert)

    orf = find_best_orf(insert, min_aa=min_orf_aa)
    cls = classify_clone(orf, L)

    # Turn the classification into the three ranges DSAP asks for,
    # applying the commandments to the coordinates.
    if orf is None:                                   # V: all 3' UTR
        utr5, orf_rng, utr3 = "N/A", "N/A", _fmt_range(1, L)
    elif cls.clone_type == "back-end partial":        # III + IV
        utr5 = "N/A"
        orf_rng = _fmt_range(1, orf.end_base)
        utr3 = _fmt_range(orf.end_base + 1, L)
    elif cls.clone_type == "front-end partial":       # VI: runs to the end
        utr5 = _fmt_range(1, orf.start_base - 1)
        orf_rng = _fmt_range(orf.start_base, L)
        utr3 = "N/A"
    elif cls.clone_type == "middle ORF partial":      # VII
        utr5, orf_rng, utr3 = "N/A", _fmt_range(1, L), "N/A"
    else:                                             # full ORF (I + II)
        utr5 = _fmt_range(1, orf.start_base - 1)
        orf_rng = _fmt_range(orf.start_base, orf.end_base)
        utr3 = _fmt_range(orf.end_base + 1, L)

    warnings = []
    if not bounds.adapter_found:
        warnings.append("adapter (GAATTC/ATTAC/G-string) NOT found - insert start "
                        "set to base 1 (fine if the sequence is already cropped)")
    if not bounds.polya_found:
        warnings.append("poly-A tail NOT found - insert end set to last base")
    if n_count:
        warnings.append(f"{n_count} N's replaced with A")
    if L < 30:
        warnings.append(f"insert is very short ({L} bp)")

    hit: Optional[BlastHit] = None
    if blast:
        print(f"[{name}] running BLASTx at NCBI ({blast_db}) - against 'nr' "
              f"this often takes 5-15 min (use --db swissprot for speed)...",
              file=sys.stderr)
        try:
            hit = run_blastx(insert, database=blast_db,
                             timeout_min=blast_timeout)
        except Exception as e:
            warnings.append(f"BLAST failed: {e}")
        else:
            if hit is None:
                warnings.append("BLASTx found no hits - likely non-coding")
            else:
                if evalue is None:
                    evalue = hit.evalue
                if orf and hit.frame and hit.frame != f"+{orf.frame}":
                    warnings.append(f"BLASTx best frame is {hit.frame} but the "
                                    f"longest ORF is in +{orf.frame} - check by hand")

    return CloneReport(
        name=name,
        raw_length=len(raw),
        n_replaced=n_count,
        adapter_found=bounds.adapter_found,
        polya_found=bounds.polya_found,
        insert_start_base=bounds.insert_start_base,
        insert_end_base=bounds.insert_end_base,
        insert_length=L,
        clone_type=cls.clone_type,
        frame=f"+{orf.frame}" if orf else "N/A",
        utr5=utr5, orf=orf_rng, utr3=utr3,
        protein=orf.protein if orf else "",
        protein_length=len(orf.protein) if orf else 0,
        blast_protein=hit.protein_name if hit else "",
        blast_organism=hit.organism if hit else "",
        blast_accession=hit.accession if hit else "",
        blast_identity_pct=hit.identity_pct if hit else None,
        blast_frame=hit.frame if hit else "",
        blast_top_hits=hit.top_hits if hit else "",
        evalue=evalue,
        evalue_ok=passes_evalue_threshold(evalue) if evalue is not None else None,
        warnings="; ".join(warnings),
        notes=cls.notes,
        insert_sequence=insert,
    )


def format_report(r: CloneReport) -> str:
    """Human-readable summary of one CloneReport."""
    lines = [
        f"=== {r.name} ===",
        f"  Raw length        : {r.raw_length} bp",
        f"  Crop raw sequence : base {r.insert_start_base} to {r.insert_end_base}"
        f"  ({r.insert_length} bp insert)",
        f"  Clone type        : {r.clone_type}   (frame {r.frame})",
        f"  5' UTR            : {r.utr5}",
        f"  ORF               : {r.orf}",
        f"  3' UTR            : {r.utr3}",
        f"  Protein ({r.protein_length} aa) : {r.protein or 'N/A'}",
    ]
    if r.blast_protein:
        lines += [
            f"  BLASTx protein    : {r.blast_protein}",
            f"  Organism          : {r.blast_organism or '?'}   "
            f"(accession {r.blast_accession}, {r.blast_identity_pct}% identity, "
            f"frame {r.blast_frame})",
            f"  Top hits          : {r.blast_top_hits}",
        ]
    if r.evalue is not None:
        lines.append(f"  E-value           : {r.evalue:g} -> "
                     f"{'PASS' if r.evalue_ok else 'FAIL'} (<= 1e-06)")
    lines.append(f"  Notes             : {r.notes}")
    if r.warnings:
        lines.append(f"  WARNINGS          : {r.warnings}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 7. BATCH MODE (files / folders)
# ---------------------------------------------------------------------------

AB1_EXTS = {".ab1", ".abi"}
TEXT_EXTS = {".fasta", ".fa", ".fas", ".fna", ".seq", ".txt"}
SUPPORTED_EXTS = AB1_EXTS | TEXT_EXTS


def read_sequences(path: Path) -> list[tuple[str, str]]:
    """Return [(name, sequence), ...] from one file.
    .ab1/.abi chromatograms need Biopython; FASTA may hold many records;
    any other text file is treated as a single bare sequence."""
    ext = path.suffix.lower()
    if ext in AB1_EXTS:
        try:
            from Bio import SeqIO
        except ImportError:
            raise RuntimeError("Reading .ab1 files needs Biopython: pip install biopython")
        rec = SeqIO.read(str(path), "abi")
        return [(path.stem, str(rec.seq))]

    # utf-8-sig strips the BOM that Notepad/PowerShell add, which would
    # otherwise hide the leading '>' of a FASTA file.
    text = path.read_text(encoding="utf-8-sig", errors="ignore")
    if text.lstrip().startswith(">"):
        records, name, chunks = [], None, []
        for line in text.splitlines():
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:].strip() or path.stem, []
            else:
                chunks.append(line.strip())
        if name is not None:
            records.append((name, "".join(chunks)))
        return records
    return [(path.stem, text)]


def collect_files(inputs: list[str]) -> list[Path]:
    files: list[Path] = []
    for item in inputs:
        p = Path(item)
        if p.is_dir():
            files.extend(sorted(f for f in p.iterdir()
                                if f.is_file() and f.suffix.lower() in SUPPORTED_EXTS))
        elif p.is_file():
            files.append(p)
        else:
            print(f"[skip] not found: {item}", file=sys.stderr)
    return files


def analyze_files(inputs: list[str], csv_path: Optional[str] = None,
                  evalue: Optional[float] = None, blast: bool = False,
                  blast_db: str = "nr",
                  blast_timeout: float = 15,
                  min_orf_aa: int = MIN_ORF_AA) -> list[CloneReport]:
    """Run the full pipeline on every sequence in the given files/folders."""
    reports: list[CloneReport] = []
    for f in collect_files(inputs):
        try:
            for name, seq in read_sequences(f):
                reports.append(analyze_clone(seq, name=name, evalue=evalue,
                                             blast=blast, blast_db=blast_db,
                                             blast_timeout=blast_timeout,
                                             min_orf_aa=min_orf_aa))
        except Exception as e:  # keep going if one file is bad
            print(f"[error] {f.name}: {e}", file=sys.stderr)
    if csv_path:
        write_csv(reports, csv_path)
    return reports


def write_csv(reports: list[CloneReport], csv_path: str) -> None:
    if not reports:
        return
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(asdict(reports[0]).keys()))
        writer.writeheader()
        for r in reports:
            writer.writerow(asdict(r))


# ---------------------------------------------------------------------------
# 8. COMMAND LINE
# ---------------------------------------------------------------------------

DEMO_RAW = (
    "NNTTGCTGACGT"                     # junk before adapter
    "GAATTCTACGGCCGGGG" "ATTAC" "GGG"  # GAATTC ... ATTAC ... G-string
    "CCATGGCTAAGGAGCTTTGA"             # insert: 5' UTR + ATG...TGA (full ORF)
    "CTCGT"                            # 3' UTR
    "AAAAAAAAAAAAAAAA"                 # poly-A tail
    "NNTT"                             # trailing junk/N's
)


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Automated DSAP clone analysis: N cleanup -> crop -> ORF -> classify.")
    ap.add_argument("inputs", nargs="*",
                    help="sequence files or folders (.ab1, .fasta, .seq, .txt)")
    ap.add_argument("-s", "--sequence", help="analyze a raw sequence given on the command line")
    ap.add_argument("-o", "--csv", help="write all results to this CSV file")
    ap.add_argument("-e", "--evalue", type=float,
                    help="BLASTx E-value to check against the 1e-06 threshold "
                         "(filled in automatically when BLAST runs)")
    ap.add_argument("--no-blast", action="store_true",
                    help="skip the online NCBI BLASTx protein identification")
    ap.add_argument("--db", default="nr",
                    help="NCBI protein database for BLASTx (default: nr; "
                         "'swissprot' is faster)")
    ap.add_argument("--blast-timeout", type=float, default=15,
                    help="minutes to wait for each NCBI BLAST job before "
                         "giving up (default: 15)")
    ap.add_argument("--min-orf", type=int, default=MIN_ORF_AA,
                    help="ignore ORFs shorter than this many amino acids "
                         f"(default: {MIN_ORF_AA})")
    ap.add_argument("--demo", action="store_true", help="run on a built-in toy sequence")
    args = ap.parse_args(argv)
    opts = dict(evalue=args.evalue, blast=not args.no_blast, blast_db=args.db,
                blast_timeout=args.blast_timeout, min_orf_aa=args.min_orf)

    if args.demo:
        reports = [analyze_clone(DEMO_RAW, name="demo", **{**opts, "blast": False,
                                                            "min_orf_aa": 1})]
    elif args.sequence:
        reports = [analyze_clone(args.sequence, name="pasted sequence", **opts)]
    elif args.inputs:
        reports = analyze_files(args.inputs, **opts)
    else:
        # Paste mode: accepts GenBank-style text with line numbers/spaces.
        if sys.stdin.isatty():
            print("Paste your sequence (line numbers and spaces are fine),\n"
                  "then press Enter on an empty line:")
        lines = []
        for line in sys.stdin:
            if not line.strip():
                if lines:
                    break
                continue
            lines.append(line)
        reports = [analyze_clone("".join(lines), name="pasted sequence", **opts)]

    if args.csv:
        write_csv(reports, args.csv)

    if not reports:
        print("No sequences found.")
        return 1
    for r in reports:
        print(format_report(r), end="\n\n")
    if args.csv:
        print(f"Wrote {len(reports)} result(s) to {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())