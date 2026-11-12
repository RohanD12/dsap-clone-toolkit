"""
test_dsap_random.py - 100 randomly generated clones with known answers.

Unlike test_dsap_toolkit.py (hand-built C/T filler), these use random
A/C/G/T sequence and random protein-coding codons, so they look much more
like real reads. Each clone is still assembled from labelled parts, so the
correct DSAP answer is known from how it was built.

A random sequence can accidentally contain a second ORF that makes the
"right" answer ambiguous, even for a human using DSAP. The generator throws
those away and re-rolls, so every clone kept has exactly one correct answer:
  * full ORF / front-end: the start ATG is the only ATG in the insert
  * back-end / non-coding: no ATG anywhere in the insert
  * back-end / front-end: the real ORF is longer than any stop-terminated
    run in the other two frames
  * non-coding: every frame hits a stop within 29 codons (below the 30-aa
    minimum ORF length)

Mix: 40 full ORF, 20 back-end partial, 20 front-end partial, 20 non-coding,
with random extras: N's in UTRs and in the poly-A tail, no adapter,
lowercase, GenBank-style formatting, junk after the tail.

Seeded, so the same 100 clones are generated every run.
Run:   python test_dsap_random.py
"""

import random
import re
import unittest

from dsap_toolkit import analyze_clone

SEED = 20261110
JUNK5 = "NNTTGCTGACGT"
ADAPTER = "GAATTCTACGGCCGGGG" "ATTAC" "GGG"
STOPS = ("TAA", "TAG", "TGA")
SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
         if a + b + c not in STOPS and a + b + c != "ATG"]


def first_stop_aa(seq: str, frame: int) -> int:
    """Codons read in `frame` before the first stop (len of protein if none)."""
    n = 0
    for i in range(frame, len(seq) - 2, 3):
        if seq[i:i + 3] in STOPS:
            return n
        n += 1
    return 10 ** 9   # no stop at all in this frame


def rand_dna(rng, n, alphabet="ACGT"):
    return "".join(rng.choice(alphabet) for _ in range(n))


def rand_utr5(rng, n):
    """Random 5' UTR that doesn't start with G (it would be read as part of
    the adapter's G-string)."""
    return (rng.choice("ACT") + rand_dna(rng, n - 1)) if n else ""


def rand_codons(rng, n):
    return "".join(rng.choice(SENSE) for _ in range(n))


def sprinkle_ns(rng, s, k):
    s = list(s)
    for _ in range(k):
        if s:
            s[rng.randrange(len(s))] = "N"
    return "".join(s)


def make_clone(rng, kind):
    """Return (raw_read, expected_fields) or None if it should be re-rolled."""
    if kind == "full ORF":
        u5 = rand_utr5(rng, rng.randint(0, 80))
        cds = "ATG" + rand_codons(rng, rng.randint(100, 300)) + rng.choice(STOPS)
        u3 = rand_dna(rng, rng.randint(0, 150))
        if rng.random() < 0.3:
            u3 = sprinkle_ns(rng, u3, rng.randint(1, 3))
        polya = True
    elif kind == "back-end partial":
        u5 = ""
        cds = rand_codons(rng, rng.randint(40, 250)) + rng.choice(STOPS)
        u3 = rand_dna(rng, rng.randint(10, 150))
        polya = True
    elif kind == "front-end partial":
        u5 = rand_utr5(rng, rng.randint(0, 80))
        cds = "ATG" + rand_codons(rng, rng.randint(60, 300))
        u3 = ""
        polya = False          # read ran out before the stop codon
    else:  # non-coding
        u5, cds = "", ""
        u3 = rand_utr5(rng, rng.randint(60, 400))
        polya = True

    insert_body = u5 + cds + u3
    cleaned = insert_body.replace("N", "A")
    L_body = len(insert_body)

    # --- re-roll anything ambiguous --------------------------------------
    if re.search("A{10,}", cleaned):
        return None
    atgs = cleaned.count("ATG")
    if kind in ("full ORF", "front-end partial"):
        if atgs != 1:
            return None
    elif atgs != 0:
        return None
    if kind in ("back-end partial", "non-coding") and cleaned[:1] == "G":
        return None
    if kind == "back-end partial":
        real = len(cds) // 3 - 1
        if first_stop_aa(cleaned, 0) != real or \
                any(first_stop_aa(cleaned, f) >= real for f in (1, 2)):
            return None
    if kind == "front-end partial":
        real = len(cds) // 3
        if any(first_stop_aa(cleaned, f) >= real for f in range(3)
               if first_stop_aa(cleaned, f) < 10 ** 9):
            return None
    if kind == "non-coding":
        if any(first_stop_aa(cleaned, f) >= 30 for f in range(3)):
            return None

    # --- assemble the read -------------------------------------------------
    tail = ""
    if polya:
        tail = "A" * rng.randint(12, 30)
        if rng.random() < 0.3:      # an N miscall inside the tail
            i = rng.randrange(3, len(tail) - 3)
            tail = tail[:i] + "N" + tail[i + 1:]
    insert = insert_body + tail
    L = len(insert)

    adapter = rng.random() < 0.9
    prefix = JUNK5 + ADAPTER if adapter else ""
    junk3 = rand_dna(rng, rng.randint(0, 20), "CGT") if polya else ""
    raw = prefix + insert + junk3

    start = len(prefix) + 1 if adapter else 1
    exp = dict(insert_start_base=start, insert_end_base=start + L - 1,
               insert_length=L, clone_type=kind)

    a, b = len(u5) + 1, len(u5) + len(cds)
    if kind == "full ORF":
        exp.update(frame=f"+{len(u5) % 3 + 1}",
                   utr5=f"1-{len(u5)}" if u5 else "N/A",
                   orf=f"{a}-{b}", utr3=f"{b + 1}-{L}")
    elif kind == "back-end partial":
        exp.update(frame="+1", utr5="N/A", orf=f"1-{b}", utr3=f"{b + 1}-{L}")
    elif kind == "front-end partial":
        exp.update(frame=f"+{len(u5) % 3 + 1}",
                   utr5=f"1-{len(u5)}" if u5 else "N/A",
                   orf=f"{a}-{L}", utr3="N/A")
    else:
        exp.update(frame="N/A", utr5="N/A", orf="N/A", utr3=f"1-{L}")

    # --- cosmetic variations that must not change the answer --------------
    style = rng.random()
    if style < 0.2:
        raw = raw.lower()
    elif style < 0.35:
        low = raw.lower()
        raw = "\n".join(f"{i + 1:>9} " + " ".join(low[i + j:i + j + 10]
                                                 for j in range(0, 60, 10))
                        for i in range(0, len(low), 60))
    return raw, exp


def generate(seed=SEED):
    rng = random.Random(seed)
    plan = (["full ORF"] * 40 + ["back-end partial"] * 20 +
            ["front-end partial"] * 20 + ["non-coding"] * 20)
    cases = []
    for i, kind in enumerate(plan, 1):
        while True:
            clone = make_clone(rng, kind)
            if clone:
                cases.append((f"{i:03d}_{kind.replace(' ', '_').replace('-', '_')}",
                              *clone))
                break
    return cases


FIELDS = ["insert_start_base", "insert_end_base", "insert_length",
          "clone_type", "frame", "utr5", "orf", "utr3"]


class TestRandomClones(unittest.TestCase):
    pass


def _make_test(raw, exp):
    def test(self):
        r = analyze_clone(raw, name="test", blast=False)
        self.assertEqual({f: exp[f] for f in FIELDS},
                         {f: getattr(r, f) for f in FIELDS})
    return test


for _name, _raw, _exp in generate():
    setattr(TestRandomClones, "test_" + _name, _make_test(_raw, _exp))


if __name__ == "__main__":
    unittest.main(verbosity=1)
