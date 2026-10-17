# Mathematical Foundations of Revised Spherical Cap Harmonic Analysis

## License

Copyright (c) 2026 R-SCHA contributors. The project-authored handbook PDF and
the accompanying notes in this directory are licensed under
[Creative Commons Attribution 4.0 International](../../LICENSES/CC-BY-4.0.txt).
See [LICENSING.md](../../LICENSING.md) for attribution and third-party exclusions.
Canonical license: <https://creativecommons.org/licenses/by/4.0/>.

This is the corrected pure-mathematics edition of the R-SCHA handbook.  It
develops the closed-shell three-dimensional theory from the Laplace equation
and complementary mixed boundary-value problems, rather than treating the
four coefficient families as preassigned formulas.

## Mathematical scope

The project contains:

- complete derivations of the three-dimensional subpotentials `V1` and `V2`;
- the Legendre, Mehler, zero-eigenvalue constant, and radial
  Sturm--Liouville systems;
- the complete four-family R-SCHA series;
- orthogonality, componentwise completeness, tensor, coordinate, flux, and
  boundary-regularity theory;
- full derivations of the exact shell energy, the exact fixed-surface energy
  with all surviving cross terms, and the simplified Legendre surface
  spectrum;
- exact coefficient projections and the spherical-harmonic-to-R-SCHA map;
- a concise, separate account of Thébault's classical one-surface R-SCHA2D;
- a notation register and explicit statements of unresolved mathematical
  questions.

The handbook intentionally omits algorithms for root finding,
special-function evaluation, iterative inversion, and other numerical
implementation details.

## Correction status

This edition removes all Pan Zhang citations and attribution.  It also fixes
the zero-mode `V2` projection normalisation, the sign of the Mehler bilinear
Wronskian identity, the null-space statement for the continuous fixed-surface
full-gradient Gram form, and the trace-space formulation of the mixed
boundary problem.  The former coordinate-frame illustration was removed;
the analytical coordinate formulas and the separate cap-rotation figure are
retained.

This corrected edition additionally repairs the Legendre cap-norm sign,
reconciles the Wigner-$D$ defining identity with the coefficient-space
transpose, replaces the incorrect non-unique-continuation argument by the
proper existence/unique-continuation/ill-posedness distinction, and separates
$L^2$ boundary projection from the trace regularity required for an
$H^1(\Omega)$ weak solution.

Known open mathematical questions are not presented as established theorems.
In particular, the triangular Mehler index set and the complete functional
status of the combined R-SCHA2D family remain explicitly qualified.

## Compilation

Run:

```bash
pdflatex main.tex
bibtex8 main
pdflatex main.tex
pdflatex main.tex
```

The compiled handbook is `main.pdf`; a ready-to-use `main.bbl` is included.

## Source organisation

The main text is divided into twenty-four section files under `sections/` and
assembled through `main.tex`.  Figures are stored in `figures/`.
