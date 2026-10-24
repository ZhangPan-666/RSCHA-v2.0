# Validation report

## Deliverable

- Title: *Mathematical Foundations of Revised Spherical Cap Harmonic Analysis*
- Edition: audited mathematical-handbook edition
- Date: October 24, 2026
- PDF pages: 167
- Main section modules: 24
- LaTeX labels: 1,222, all unique
- Missing cross-reference labels: 0
- Missing bibliography keys: 0
- Missing referenced images: 0

## Structural and layout checks

- LaTeX compilation with TeX Live 2026: passed.
- BibTeX8 bibliography pass: passed.
- Undefined-reference and citation warnings: 0 after the final rerun.
- Overfull/underfull box warnings: 0.
- PDF rendered page by page with Poppler for visual inspection: passed.
- Page size: US Letter; PDF openable and searchable.
- Latest vector figure links: 17/17 resolved from `figures/vector/`.

## Mathematical audit

The rectangular Mehler expansion now uses independent positive radial index
`p` and azimuthal order `m`. The legacy triangular selection is explicitly
restricted and is proved non-dense in the full lateral-data and energy spaces
by the missing `(p,m)=(1,2)` mode. The two Legendre recurrence identities in
Section 9, Schmidt-normalised counterparts, Wigner-D Euler phase convention,
repeated-root radial wording, polar-axis component statement, and ND
complementary orthogonality statement were corrected.

The independent script `reproducibility/verify_handbook.py` performs 450
representative high-precision, quadrature, field, rotation, and energy checks.
The JSON result reports `PASS`; the largest scaled discrepancy is approximately
`1.52e-13`, from finite quadrature of a volume-energy block.

The associated-Legendre recurrence convention was cross-checked against the
[NIST DLMF definitions](https://dlmf.nist.gov/14.3) and
[recurrence relations](https://dlmf.nist.gov/14.10).

## Editorial item intentionally deferred

The AI-use disclosure remains deferred as previously instructed and is not
present in the manuscript. It remains listed in the audit report as a
submission-stage item.

## Scope and limits

The numerical checks cover representative finite orders and parameter values;
they are not a formal proof for every possible cap aperture, shell thickness,
high-degree truncation, or external C++ implementation path. The research
questions retained in Section 23 concern restricted-index approximation,
edge singularities, fixed-surface spectral diagonalisation, R-SCHA2D
continuation spaces, and non-spherical boundaries.
