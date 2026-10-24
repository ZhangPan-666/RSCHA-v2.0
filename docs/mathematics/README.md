# Mathematical Foundations of Revised Spherical Cap Harmonic Analysis

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

The complete Mehler expansion uses independent rectangular bounds
$1\le p\le P_{\max}$ and $0\le m\le M_{\mathrm{Meh}}$. The legacy
triangular rule $m\le p$ is retained only as an explicitly restricted model;
its non-density in the full lateral-data and energy spaces is demonstrated by
the missing $(p,m)=(1,2)$ mode. The combined R-SCHA2D continuation family
remains qualified as an open functional-analytic question.

## Current audited build

The reviewed edition was compiled on October 24, 2026 with TeX Live 2026.
It has 167 pages, 24 section modules, 1,222 unique labels, no missing
cross-references or bibliography keys, and no overfull/underfull box warnings.
The mathematical audit includes 450 independent high-precision or quadrature
comparisons; details are in `notes/handbook_audit_20261024.md` and
`reproducibility/handbook_math_validation.json`.

The AI-use disclosure remains intentionally deferred at the author's request
and is recorded as a pending editorial item in the audit report.

## Compilation

Run:

```bash
pdflatex main.tex
bibtex8 main
pdflatex main.tex
pdflatex main.tex
```

The compiled handbook is `main.pdf`; a ready-to-use `main.bbl` is included. The portable `build.ps1` script compiles into `build/` and copies the PDF back to the root.

## Source organisation

The main text is divided into twenty-four section files under `sections/` and
assembled through `main.tex`.  Figures are stored in `figures/`.


