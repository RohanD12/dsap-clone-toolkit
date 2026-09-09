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
  5. run_blastx() / summarize_blast_hits()
                                - submit a sequence to NCBI BLAST and pull
                                  out the top hit per distinct organism,
                                  formatted ready to paste into a DSAP table
  6. passes_evalue_threshold() - quick E-value sanity check

Requires: biopython  (pip install biopython --break-system-packages)

Everything is written as plain functions so you can import just the
pieces you want, or run this file directly to see a demo on a toy
sequence.
"""

from __future__ import annotations
import re
import time
from dataclasses import dataclass
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
    # any stray A-runs upstream of the actual insert).
    search_region = seq[insert_start - 1:]
    polya_match = POLYA_PATTERN.search(search_region)
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


def find_best_orf(sequence: str) -> Optional[ORFResult]:
    """
    Translate all 3 forward reading frames and return the single longest
    ORF found (start codon M ... stop codon *, or the longest M...end /
    start...* run if no full M-to-* ORF exists, matching the "back end
    partial" / "front end partial" cases DSAP describes).
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
# 6. BLAST SUBMISSION + TABLE SUMMARY  (requires biopython + internet)
# ---------------------------------------------------------------------------

def run_blastx(sequence: str, database: str = "nr", hitlist_size: int = 50):
    """
    Submit `sequence` to NCBI's BLASTx (protein) search and return the
    parsed Biopython BLAST record. This is a network call and can take
    a minute or two, mirroring the "Searching..." wait in the DSAP UI.

    Requires: from Bio.Blast import NCBIWWW, NCBIXML
    """
    from Bio.Blast import NCBIWWW, NCBIXML

    result_handle = NCBIWWW.qblast(
        "blastx", database, sequence,
        hitlist_size=hitlist_size,
    )
    record = NCBIXML.read(result_handle)
    return record


def summarize_blast_hits(record, max_organisms: int = 3) -> list[dict]:
    """
    Walk a parsed BLAST record and pull the best hit from each of up to
    `max_organisms` DISTINCT organisms -- exactly the "list the best
    matches from three different organisms" step DSAP asks for.

    Returns a list of dicts ready to drop into a DSAP table:
    accession, definition, organism, start, end, evalue
    """
    seen_organisms = set()
    rows = []

    for alignment in record.alignments:
        title = alignment.title
        # Titles look like: "gi|123|ref|XM_1.1| PREDICTED: ... [Genus species]"
        organism_match = re.search(r"\[([^\]]+)\]", title)
        organism = organism_match.group(1) if organism_match else "unknown"

        if organism in seen_organisms:
            continue

        hsp = alignment.hsps[0]  # best HSP for this alignment
        rows.append({
            "accession": alignment.accession,
            "definition": title.split("|")[-1].split("[")[0].strip()[:60],
            "organism": organism,
            "start": hsp.query_start,
            "end": hsp.query_end,
            "evalue": hsp.expect,
        })
        seen_organisms.add(organism)

        if len(rows) >= max_organisms:
            break

    return rows


def print_blast_table(rows: list[dict]) -> None:
    """Pretty-print a summarize_blast_hits() result as a simple table."""
    if not rows:
        print("No hits found.")
        return
    print(f"{'Accession':<15}{'Organism':<25}{'Start':<8}{'End':<8}{'E-value':<10}")
    for r in rows:
        print(f"{r['accession']:<15}{r['organism']:<25}{r['start']:<8}{r['end']:<8}{r['evalue']:<10}")


# ---------------------------------------------------------------------------
# DEMO
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # A toy raw sequence: junk + adapter (GAATTC...ATTAC...GGG) + insert
    # with an ORF (ATG...TAA) + poly-A tail.
    demo_raw = (
        "NNTTGCTGACGT"                     # junk before adapter
        "GAATTCTACGGCCGGGG" "ATTAC" "GGG"  # GAATTC ... ATTAC ... G-string
        "ATGGCTAAGGAGCTTTGA"                # insert: ATG...TAA (full ORF)
        "AAAAAAAAAAAAAAAA"                  # poly-A tail
        "NNTT"                              # trailing junk/N's
    )

    print("=== 1. Insert boundary detection ===")
    boundaries = find_insert_boundaries(demo_raw)
    print(boundaries)

    print("\n=== 2. N -> A cleanup ===")
    cleaned, n_count = clean_ns(demo_raw)
    print(f"Replaced {n_count} N's.")
    print(cleaned)

    print("\n=== 3. ORF finder (on the cropped insert) ===")
    insert_seq = demo_raw[boundaries.insert_start_base - 1: boundaries.insert_end_base]
    orf = find_best_orf(insert_seq)
    print(orf)

    print("\n=== 4. Clone classification ===")
    classification = classify_clone(orf, len(insert_seq))
    print(classification)

    print("\n=== 5. E-value check ===")
    print("2e-40 passes threshold:", passes_evalue_threshold(2e-40))
    print("0.5 passes threshold:", passes_evalue_threshold(0.5))
    print("Blast Results (this may take a minute or two)...")
    print(record = run_blastx(insert_seq))
    print(rows = summarize_blast_hits(record))
    print(print_blast_table(rows))