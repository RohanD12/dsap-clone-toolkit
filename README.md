# dsap-clone-toolkit
Python toolkit that automates manual cDNA clone analysis — adapter/insert boundary detection, sequence cleanup, ORF discovery across reading frames, and rule-based clone classification (full ORF, partial, non-coding). Built for the WSSP/DSAP bioinformatics workflow.

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
  6. passes_evalue_threshold() - quick E-value sanity check

Everything is written as plain functions so you can import just the
pieces you want, or run this file directly to see a demo on a toy
sequence.
"""
