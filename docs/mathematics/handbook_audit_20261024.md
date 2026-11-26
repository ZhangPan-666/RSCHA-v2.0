# Mathematical handbook audit (2026-10-24)

## Results

The review covered all 24 sections of *Mathematical Foundations of Revised Spherical Cap Harmonic Analysis*, including captions, cross-references, bibliography, and the final PDF. The manuscript was rebuilt to 167 pages using the vector artwork in `figures/vector/*.pdf`.

- LaTeX compiled successfully with all cross-references and citations resolved.
- Automated checks found 1,222 labels, zero duplicate labels, zero missing references or figures, and balanced environments.
- The log contained no overfull or underfull boxes or citation or reference warnings.
- All 167 pages were rendered and checked in thumbnail form, including the geometry, index, 3D basis, derivative, rotation, and reference pages.

## Mathematical changes

1. **Rectangular Mehler indexing.** The text uses independent ranges `1 ≤ p ≤ P_max` and `0 ≤ m ≤ M_Meh`. The count is `P_max(2M_Meh+1)`, giving 36 real modes at `P_max = M_Meh = 4`. The legacy restricted model with `m ≤ p` has its own count.
2. **Completeness of the triangular truncation.** Section 13 gives the omitted mode `Ψ_* = Ψ^c_{1,2}`. It is energy-orthogonal to the triangular and constant families, giving a strictly positive lower bound on its distance from the triangular family. Its approximation rates can be studied within the closure of that family, a proper subspace.
3. **Legendre recurrence.** Section 9 gives the degree recurrence first and uses `sqrt(1-x^2)` in the second denominator. The Schmidt semi-normalized forms follow the same convention. The order-lowering identity applies at `m ≥ 1`.
4. **Wigner-D convention.** The defining identity fixes the coefficient transpose and the phases on the first and last Euler angles. The mapping section specifies active and passive matrix conventions.
5. **Boundary and radial details.** The repeated root `μ = 1/4` belongs to the degenerate one-surface construction. The statement about the noninteger exterior branch is formulated on the finite shell. At the axis, spherical and NED components depend on the meridian, so the unique limit is evaluated in a fixed Cartesian frame.
6. **Alternative boundary conditions.** A Green identity establishes strict gradient orthogonality for the complementary lateral-Neumann and radial-Dirichlet split. The spectrum, weights, and edge conditions require their own derivation.
7. **Artwork.** The 17 figures use the latest vector PDFs. The dense Legendre and Mehler atlases use landscape pages, and the three constant-family figures share a page.

## Independent numerical verification

`reproducibility/verify_handbook.py` uses 45-digit working precision to check Legendre recurrences and derivatives, root norms, SH boundary projection, Mehler bilinear and parameter norms, constant-family Lerch norms, radial Neumann eigenfunctions and closed-form overlaps, two-boundary coefficient recovery, Cartesian fields and Hessians, NED rotation for six tensor components, Wigner-D rotation, volume energy, and fixed-radius surface Gram blocks. All 450 comparisons passed. The largest scaled error was about `1.52e-13` in a volume-energy quadrature. Results are in `reproducibility/handbook_math_validation.json`.

Branch and recurrence conventions follow [NIST DLMF 14.3](https://dlmf.nist.gov/14.3) and [14.10](https://dlmf.nist.gov/14.10).

## Scope and editorial items

- The author deferred the AI-use statement for this revision. It remains an editorial item for submission.
- Section 23 lists questions involving R-SCHA2D continuation spaces, edge singularities, fixed-surface spectral diagonalization, and nonspherical boundaries.
- The numerical verification covers representative parameters and finite-order identities. A uniform proof over all degrees, cap apertures, or C++ execution paths is outside its scope.
