# Reproducibility Code for Two-Periodic Mixed-Error Codes

**Associated manuscript:** *Adjacent-Transposition Syndromes and Sharp Length Bounds for Two-Periodic Mixed-Error Codes over Finite Fields*

**Verification-software author:** Malihe Aliasgari

**Repository presented by:** Malihe Aliasgari and Yousef Nejatbakhsh

This public repository contains Python verification programs, saved computational results, and scientific figures associated with the manuscript. The verification software was developed by **Malihe Aliasgari**; repository presentation and software authorship are separate from the research manuscript's author list. The results concerning optimal lengths refer to the specified **two-periodic family**, not to every possible error-correcting code.

## Reproducing Section 6.5

The computations were run using **Python 3.13** (reference environment: Python **3.13.5**) and **only the Python standard library**. No third-party packages, datasets, or network connection are needed.

From the repository's root directory, run these commands in the order used for the paper's Section 6.5:

```bash
python3 GF64_Sharp_Bound_Verification.py
python3 Generalized_Two_Periodic_Bounds_Verification.py 18
python3 Two_Periodic_Exhaustive_Verification.py 12
```

**Use `12` for the third command** to reproduce the full exhaustive table through extension degree `k = 12`. Omitting the argument makes that program use its shorter default of `k = 10` and **does not** reproduce the complete Section 6.5 run. The extended computation can take substantially longer than the other two.

### What each program checks

| Program | Main checks |
| --- | --- |
| [`GF64_Sharp_Bound_Verification.py`](GF64_Sharp_Bound_Verification.py) | Theorem 5.2: the sharp `m_max(64) = 29` result, 13 candidate residue pairs, cyclotomic-polynomial GCD certificates, and an even-length construction. |
| [`Generalized_Two_Periodic_Bounds_Verification.py`](Generalized_Two_Periodic_Bounds_Verification.py) | Table 1 for `k = 2, ..., 18`, using the proved short-exponent reduction of Lemma 5.3 for larger fields. |
| [`Two_Periodic_Exhaustive_Verification.py`](Two_Periodic_Exhaustive_Verification.py) | The GF(7) example (Appendix A); the GF(64) sharp result; parameter comparisons `6`, `22`, and `29`; encoder, decoder, and DNA examples; reject-and-re-check decoder logic; Tables 1 and 2 through `k = 12`; finite checks related to Proposition 5.9 and Lemma 4.1/Corollary 4.2; the `k >= 29` threshold; and a numerical check of Equation (30). |

### Additional checks in the updated exhaustive program

- **GF(7), Appendix A:** Checks odd- and even-position adjacent swaps, the general locator equation, failure of the characteristic-two shortcut, and the explicit collision illustrating the necessity of the length condition.
- **Fixed-alpha optimization over GF(64):** With `alpha = t`, the fixed choice `beta = alpha^2` yields a maximum separated length of `6`. Optimizing `beta = alpha^b` while retaining `alpha = t` yields a maximum information length of **22**, achieved for `b = 8`. Jointly optimizing `alpha` and `beta` yields **29**.
- **Decoder safety checks:** The demonstration rejects unusable syndrome coordinates, rejects unrecognized projective locators, and **recomputes both syndromes after any proposed correction**, rejecting candidates that fail this re-check. A seeded simulation also runs the decoder on 20,000 arbitrary words and verifies the reporting behavior. These checks do not imply correction guarantees for error patterns outside the model proved in the manuscript.

## Saved Outputs

The [`outputs/`](outputs/) directory contains the saved console results:

| Output file | Corresponding run |
| --- | --- |
| [`GF64_Sharp_Bound_Verification.txt`](outputs/GF64_Sharp_Bound_Verification.txt) | First command above |
| [`Generalized_Two_Periodic_Bounds_Verification.txt`](outputs/Generalized_Two_Periodic_Bounds_Verification.txt) | Second command above, with `18` |
| [`Two_Periodic_Exhaustive_Verification_K12.txt`](outputs/Two_Periodic_Exhaustive_Verification_K12.txt) | Third command above, with `12` — **reference log for Section 6.5** |
| [`Two_Periodic_Exhaustive_Verification.txt`](outputs/Two_Periodic_Exhaustive_Verification.txt) | Optional shorter run, with `10` |

To regenerate these output files, use ordinary shell redirection:

```bash
python3 GF64_Sharp_Bound_Verification.py > outputs/GF64_Sharp_Bound_Verification.txt 2>&1
python3 Generalized_Two_Periodic_Bounds_Verification.py 18 > outputs/Generalized_Two_Periodic_Bounds_Verification.txt 2>&1
python3 Two_Periodic_Exhaustive_Verification.py 12 > outputs/Two_Periodic_Exhaustive_Verification_K12.txt 2>&1
```

Timing information, when present in older reference logs, was measured separately; simple redirection does not add timing lines.

## Research Visualizations

### GF(64) Syndrome-Direction Separation

![GF(64) syndrome-direction separation](figures/GF64_Syndrome_Separation.png)

*For `alpha = t^2 + t + 1`, `beta = alpha^6`, and 29 information symbols, the 57 permitted error directions are distinct: 29 substitution directions (blue), 28 adjacent-transposition directions (orange), and six unused projective directions (gray). The two check symbols are assumed intact.*

### Computed Optimal-Length Gaps, Extension Degrees 2–18

![Computed optimal-length gaps for extension degrees 2 through 18](figures/Computed_Optimal_Length_Gaps.png)

*The plotted quantity is `2^(k-1) - m_max(2^k)` for the specified two-periodic family, using the computed values reported in Table 1.*

## Computational Scope and Precision

- **Exact results:** The finite-field and polynomial computations, including the GF(64) certificates, parameter searches, GF(7) example, and decoder arithmetic, use exact arithmetic.
- **Large-field search:** The generalized script applies the proved short-exponent restriction of Lemma 5.3; its `k = 18` computation is not a brute-force search through every parameter.
- **Finite corroboration:** Exhaustive calculations through `k = 12` support the tabulated results and test Proposition 5.9 for those degrees. They do not replace the general mathematical proof.
- **Numerical checks:** The threshold test uses a real-valued numerical inequality. The Equation (30) character-sum check uses seeded pseudorandom finite sets and floating-point complex arithmetic. It is a numerical spot-check of the corrected identity, not a proof of the general estimate.
- **Numbering:** Theorem, proposition, table, and equation references reflect the manuscript draft associated with this release; check them against the final journal version.

## Repository Structure

```text
Two-Periodic-Mixed-Error-Codes/
├── README.md
├── GF64_Sharp_Bound_Verification.py
├── Generalized_Two_Periodic_Bounds_Verification.py
├── Two_Periodic_Exhaustive_Verification.py
├── CITATION.cff
├── LICENSE
├── SHA256SUMS
├── figures/
│   ├── GF64_Syndrome_Separation.png
│   └── Computed_Optimal_Length_Gaps.png
└── outputs/
    ├── GF64_Sharp_Bound_Verification.txt
    ├── Generalized_Two_Periodic_Bounds_Verification.txt
    ├── Two_Periodic_Exhaustive_Verification.txt
    └── Two_Periodic_Exhaustive_Verification_K12.txt
```

The manuscript PDF/LaTeX source, temporary editing scripts, and cache files are not part of this reproducibility archive.

## Citation and Archiving

The manuscript refers to **release `v1.0`** for the Section 6.5 computations. Publish that tag only **after** updating the executable file, the `K12` output log, and this README; never reassign a published version tag.

Zenodo can archive the published release and assign a DOI. Once a DOI has actually been issued, add it to the software reference and Data Availability Statement.

**Verification-software author:** Malihe Aliasgari.

**Repository presented by:** Malihe Aliasgari and Yousef Nejatbakhsh.

See [`CITATION.cff`](CITATION.cff) for software-citation metadata.


### Archived Software Release

**Version:** v1.0

**Zenodo DOI:** https://doi.org/10.5281/zenodo.24170395

**GitHub Release:** https://github.com/malihealiasgari20/Two-Periodic-Mixed-Error-Codes/releases/tag/v1.0

The archived release contains the Python verification programs, computational results, and supporting scientific figures accompanying the manuscript "Adjacent-Transposition Syndromes and Sharp Length Bounds for Two-Periodic Mixed-Error Codes over Finite Fields."

Verification software developed by Malihe Aliasgari. Repository presented by Malihe Aliasgari and Yousef Nejatbakhsh.


## License

MIT License; see [`LICENSE`](LICENSE).

Copyright (c) 2026 Malihe Aliasgari.
