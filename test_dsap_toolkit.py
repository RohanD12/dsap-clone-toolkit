"""
test_dsap_toolkit.py - 20 synthetic clones with known answers, plus
regression cases for bugs found along the way.

Each clone is built from labelled parts (adapter, 5' UTR, coding region,
3' UTR, poly-A tail), so the correct DSAP answer is known exactly from
how the sequence was built. BLAST is turned off (blast=False) so the tests
run instantly and never touch NCBI.

Run:   python test_dsap_toolkit.py        (or: python -m pytest test_dsap_toolkit.py)

Design notes so the "known answer" really is the only answer:
  * UTR filler uses only C and T, so it can never contain ATG or a stop
    codon (TAA/TAG/TGA all need an A or G), in any reading frame.
  * Coding regions are ATG + GCT/CCT repeats + stop. Shifted by one or two
    bases these read CTG/TGC or CTC/TCC: no extra ATGs or stops.
  * 3' UTRs end in C, so the UTR -> poly-A junction can't form TAA.
"""

import re
import unittest

from dsap_toolkit import analyze_clone

JUNK5 = "NNTTGCTGACGT"                          # junk before the adapter
ADAPTER = "GAATTCTACGGCCGGGG" "ATTAC" "GGG"      # GAATTC ... ATTAC ... G-string
PREFIX = JUNK5 + ADAPTER                         # insert starts right after this
POLYA = "A" * 18
JUNK3 = "CTGTCG"                                 # junk after the poly-A tail


def utr(n: int) -> str:
    """n bases of C/T filler that ends in C (no ATG, no stops, any frame)."""
    s = ("CTTCC" * 100)[:n]
    return s[:-1] + "C" if n else ""


def full_orf(codons: int, stop: str = "TAA", aa: str = "GCT") -> str:
    return "ATG" + aa * codons + stop


def rng(a: int, b: int) -> str:
    return f"{a}-{b}"


def build(utr5="", cds="", utr3="", adapter=True, polya=True, junk3=True):
    """Assemble a raw read and return (raw, expected-insert-info)."""
    insert = utr5 + cds + utr3 + (POLYA if polya else "")
    raw = (PREFIX if adapter else "") + insert + (JUNK3 if polya and junk3 else "")
    start = len(PREFIX) + 1 if adapter else 1
    return raw, dict(insert_start_base=start,
                     insert_end_base=start + len(insert) - 1,
                     insert_length=len(insert))


def expect_full(utr5, cds, utr3, **kw):
    raw, exp = build(utr5, cds, utr3, **kw)
    a, b = len(utr5) + 1, len(utr5) + len(cds)
    exp.update(clone_type="full ORF",
               frame=f"+{len(utr5) % 3 + 1}",
               utr5=rng(1, len(utr5)) if utr5 else "N/A",
               orf=rng(a, b),
               utr3=rng(b + 1, exp["insert_length"]))
    return raw, exp


# ---------------------------------------------------------------------------
# The 20 clones
# ---------------------------------------------------------------------------
CASES = {}

# 1-3: full ORF in each reading frame
CASES["01 full ORF, frame +1"] = expect_full(utr(9), full_orf(40), utr(25))
CASES["02 full ORF, frame +2"] = expect_full(utr(10), full_orf(40), utr(25))
CASES["03 full ORF, frame +3"] = expect_full(utr(11), full_orf(40), utr(25))

# 4: ORF starts at base 1 -> no 5' UTR
CASES["04 full ORF, no 5' UTR"] = expect_full("", full_orf(50), utr(20))

# 5-6: other stop codons (commandment I: stop is part of the ORF)
CASES["05 full ORF, TAG stop"] = expect_full(utr(14), full_orf(35, "TAG"), utr(30))
CASES["06 full ORF, TGA stop"] = expect_full(utr(22), full_orf(60, "TGA"), utr(12))

# 7: 3' UTR is only the poly-A tail (stop codon right before the A's)
CASES["07 full ORF, 3' UTR = poly-A only"] = expect_full(utr(8), full_orf(30), "")

# 8: long clone (~1 kb)
CASES["08 full ORF, ~1 kb clone"] = expect_full(utr(60), full_orf(300), utr(80))

# 9: two full ORFs in different frames -> the longer one wins
_utr5 = utr(6) + full_orf(10, aa="CCT") + "C"      # short decoy ORF, frame +1
_raw, _exp = expect_full(_utr5, full_orf(70), utr(20))   # real ORF, frame +2
CASES["09 two ORFs, longest wins"] = (_raw, _exp)

# 10: full ORF beats a LONGER back-end partial in another frame
#     ORF = bases 2-112 (frame +2, 36 aa). A TAG placed in frame +1 at
#     base 205 gives frame +1 a 68-aa back-end partial.
_raw, _exp = expect_full("C", full_orf(35), "CC" + "CTC" * 30 + "TAG" + "CC")
CASES["10 full ORF beats longer partial"] = (_raw, _exp)

# 11-12: back-end partials (commandments III + IV: ORF starts at base 1)
for name, cds, u3 in [("11 back-end partial, TAA", "CCT" * 45 + "TAA", utr(30)),
                      ("12 back-end partial, TGA", "CCT" * 35 + "TGA", utr(50))]:
    raw, exp = build("", cds, u3)
    exp.update(clone_type="back-end partial", frame="+1", utr5="N/A",
               orf=rng(1, len(cds)), utr3=rng(len(cds) + 1, exp["insert_length"]))
    CASES[name] = (raw, exp)

# 13-14: front-end partials (commandment VI: no 3' UTR) - read runs out
#        before the stop codon, so there's no poly-A tail either
for name, u5, codons in [("13 front-end partial, frame +1", utr(12), 80),
                         ("14 front-end partial, frame +3", utr(17), 55)]:
    cds = "ATG" + "GCT" * codons
    raw, exp = build(u5, cds, "", polya=False)
    exp.update(clone_type="front-end partial", frame=f"+{len(u5) % 3 + 1}",
               utr5=rng(1, len(u5)), orf=rng(len(u5) + 1, exp["insert_length"]),
               utr3="N/A")
    CASES[name] = (raw, exp)

# 15: non-coding clone with no ORF signals (commandment V: all 3' UTR)
raw, exp = build("", "", utr(120))
exp.update(clone_type="non-coding", frame="N/A", utr5="N/A", orf="N/A",
           utr3=rng(1, exp["insert_length"]))
CASES["15 non-coding, no ATG/stops"] = (raw, exp)

# 16: N's inside the insert (bad base calls in the 3' UTR) -> replaced with A
_u3 = "CCTCNCCTTC" + utr(20)
CASES["16 N's inside 3' UTR"] = expect_full(utr(9), full_orf(40), _u3)

# 17: an N in the middle of the poly-A tail must not split the tail
raw, exp = expect_full(utr(9), full_orf(40), utr(25))
raw = raw.replace(POLYA, "A" * 8 + "N" + "A" * 9)
CASES["17 N inside poly-A tail"] = (raw, exp)

# 18: already-cropped sequence (no adapter) -> insert starts at base 1
CASES["18 no adapter (pre-cropped)"] = expect_full(utr(15), full_orf(40), utr(20),
                                                   adapter=False)

# 19: GenBank-style paste: lowercase, line numbers, spaces
raw, exp = expect_full(utr(13), full_orf(45), utr(25))
_low = raw.lower()
raw = "\n".join(f"{i + 1:>9} " + " ".join(_low[i + j:i + j + 10]
                                         for j in range(0, 60, 10))
                for i in range(0, len(_low), 60))
CASES["19 GenBank-formatted paste"] = (raw, exp)

# 20: no poly-A tail at all (insert runs to the last base)
raw, exp = expect_full(utr(10), full_orf(40), utr(30), polya=False)
CASES["20 no poly-A tail"] = (raw, exp)


# ---------------------------------------------------------------------------
# Regression cases: realistic sequences that used to give wrong answers
# ---------------------------------------------------------------------------

def expect_noncoding(insert_body):
    raw, exp = build("", insert_body, "")
    exp.update(clone_type="non-coding", frame="N/A", utr5="N/A", orf="N/A",
               utr3=rng(1, exp["insert_length"]))
    return raw, exp


# 21: real non-coding DNA has stop codons everywhere. This used to be
#     called a back-end partial with a 2-aa "ORF" (bases 1-9).
CASES["21 non-coding with stray stops"] = expect_noncoding(
    "CTAGCTTGACTAACCTTAGCCTTGATCCTAACTTCC" * 3)

# 22: a tiny stop-terminated stretch at the start (15 aa) is noise, not a
#     back-end partial - it's below the 30-aa minimum ORF length.
CASES["22 15-aa back-end partial ignored"] = expect_noncoding(
    "CCT" * 15 + "TGA" + utr(60))

# 23: a run of 10+ A's inside the 3' UTR. The crop used to stop at that
#     run instead of the real poly-A tail at the end (154 instead of 180).
CASES["23 A-run inside 3' UTR"] = expect_full(
    utr(5), full_orf(30), "CTTC" + "A" * 12 + "CTTCCTTC")


FIELDS = ["insert_start_base", "insert_end_base", "insert_length",
          "clone_type", "frame", "utr5", "orf", "utr3"]


class TestKnownClones(unittest.TestCase):
    pass


def _make_test(raw, exp):
    def test(self):
        r = analyze_clone(raw, name="test", blast=False)
        got = {f: getattr(r, f) for f in FIELDS}
        want = {f: exp[f] for f in FIELDS}
        self.assertEqual(want, got)
    return test


for _name, (_raw, _exp) in CASES.items():
    _slug = "test_" + re.sub(r"\W+", "_", _name).strip("_")
    setattr(TestKnownClones, _slug, _make_test(_raw, _exp))


if __name__ == "__main__":
    unittest.main(verbosity=2)
