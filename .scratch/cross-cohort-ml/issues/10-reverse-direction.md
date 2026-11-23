# 10: Reverse direction (PDBP → PPMI)

**What to build:** The whole pipeline trained and frozen in PDBP and scored in PPMI, including the baseline ladder and the cross-sex transfer-gap verdict, reported as a secondary analysis. The gene list and scaling for this direction come from PDBP alone.

**Blocked by:** 06

**Status:** ready-for-agent

- [ ] The leakage guards hold with the cohort roles swapped (no PPMI sample influences the PDBP-trained frozen models).
- [ ] Ladder and transfer-gap tables for this direction, labelled secondary.
