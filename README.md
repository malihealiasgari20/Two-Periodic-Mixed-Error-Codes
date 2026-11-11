# Reproducibility code for two-periodic mixed-error codes

These verification programs were developed by **Malihe Aliasgari** to support the manuscript **Adjacent-Transposition Syndromes and Sharp Length Bounds for Two-Periodic Mixed-Error Codes over Finite Fields**. The software authorship presented here is separate from the manuscript author list.

**Requirements:** Python **3.8 or later**, standard library only; no external packages, datasets, or network access. The commands below use `python3` (use `python` if that is how Python 3 is installed on your computer). The indicated runtimes were observed with Python 3.13.5 in one reference environment; other computers may be slower or faster.

| Program | Manuscript results verified | Run | Typical reference runtime |
| --- | --- | --- | --- |
| `GF64_Sharp_Bound_Verification.py` | Theorem 5.2: sharp `m_max(64)=29`, 13 formal `(b,D)` cases and cyclotomic polynomial GCD certificates, 28-symbol even-length example. | `python3 GF64_Sharp_Bound_Verification.py` | ~1 second |
| `Generalized_Two_Periodic_Bounds_Verification.py` | Table 1, `k=2,...,18`, using the short-exponent reduction in Lemma 5.3 to certify the maximum lengths. | `python3 Generalized_Two_Periodic_Bounds_Verification.py 18` | ~19 seconds |
| `Two_Periodic_Exhaustive_Verification.py` | Tables 1 and 2 (default `k=2,...,10`; can extend to 12), Proposition 5.9, Lemma 4.1/Corollary 4.2 checks, `GF(64)` encoder/decoder/DNA examples, the `k>=29` numerical threshold, and a numerical check of Equation (30). | `python3 Two_Periodic_Exhaustive_Verification.py 10` | ~10 seconds |

## Saved results

Text outputs from all three commands are in [`outputs/`](outputs/). Each ends with an `Elapsed:` line measured in the reference environment. To rerun a program and save new output, for example:

```bash
python3 GF64_Sharp_Bound_Verification.py > outputs/GF64_Sharp_Bound_Verification.txt 2>&1
python3 Generalized_Two_Periodic_Bounds_Verification.py 18 > outputs/Generalized_Two_Periodic_Bounds_Verification.txt 2>&1
python3 Two_Periodic_Exhaustive_Verification.py 10 > outputs/Two_Periodic_Exhaustive_Verification.txt 2>&1
```

An additional `outputs/Two_Periodic_Exhaustive_Verification_k12.txt` records `python3 Two_Periodic_Exhaustive_Verification.py 12`, which passed and took about 96 seconds in the reference environment. The original saved outputs include `/usr/bin/time` wall-clock durations; the above shell redirections do not automatically add timing lines. All code aborts with an assertion or error if an explicitly asserted verification fails.

## Scope and precision

- The **GF(64)** certificate and finite-field searches use exact binary-polynomial arithmetic. No randomization or floating-point approximation enters these code-length results.
- The **generalized** script relies on the *proved* short-exponent restriction of Lemma 5.3 for the large-field search; it is not an unrestricted search through all exponents or codewords at `k=18`.
- The **exhaustive** script enumerates primitive-field elements up to Frobenius conjugacy, all admissible alternating weights, and all candidate lengths (using the monotonicity of direction separation). It tests the equality in Proposition 5.9 for finite degrees, rather than replacing its mathematical proof.
- The **Jacobi-sum check** (Part 7) uses deterministic seeded pseudorandom sets and floating-point complex values. This is a numerical spot-check of the corrected conjugated character identity, **not** an exact verification of the general character-sum lower-bound theorem.
- The **threshold check** (Part 6) evaluates the explicit real-valued inequality numerically for `k=26,...,32`; its general mathematical justification belongs to the manuscript.
- References to **Theorem 5.2**, **Lemma 5.3**, **Proposition 5.9**, **Tables 1 and 2**, and **Equation (30)** correspond to the 26-page preview manuscript dated November 2026. Recheck numbering after the final MDPI-template compilation.

## Citation and archiving

The files in this folder are software and output logs only; the manuscript's `.tex` and PDF, editing scripts, and Python cache files are intentionally excluded.

To obtain a permanent software DOI, publish the folder as a public GitHub repository, enable Zenodo's GitHub integration for that repository, and create a GitHub release (for example, `v1.0.0`). Add the GitHub URL and minted Zenodo DOI to the manuscript's Data Availability Statement **after they exist**. `CITATION.cff` can then be updated with the published software DOI.

## License

MIT License; see [LICENSE](LICENSE).
