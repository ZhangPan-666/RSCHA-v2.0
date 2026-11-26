# R-SCHA mathematical handbook review and revision

Review date: 2026-11-06. The baseline was the 165-page edition after the previous figure revision.

## Conclusions

The review confirmed substantial parts of the separation of variables, radial spectrum, special-function identities, and energy formulas. It identified missing assumptions, overly broad statements, function-space issues, bibliographic errors, and repetition. Several mathematical problems had explicit counterexamples.

The revised handbook is 160 pages. It corrects the issues below, supplies proofs and conditions, and removes repeated derivations and unsupported claims about research status. The structural and numerical checks passed for the formulas and cases they covered.

## Method and coverage

The review covered all 24 sections, including boundary-value and spectral theory, energy and tensor formulas, coordinates and flux, completeness and mapping, projection and convergence, and the surrounding text, notation, and references.

- The analytical review checked domains, normals, regular branches, gauge constants, weak data spaces, convergence topologies, and the distinction between finite and infinite expansions. Explicit counterexamples were used where appropriate.
- The numerical review repeated the 45-digit checks and added independent checks of special functions, radial modes, Cartesian harmonic fields, energy, flux, and Wigner rotation.
- Publication records, original papers, and NIST material were used to verify metadata and conventions.
- Editorial checks covered labels, references, LaTeX environments, compilation, and rendering of revised pages, including tables and formulas.

The review covered the handbook and representative numerical checks. The 3D mode figures and boxed axes from the previous figure revision were retained, and the earlier Mehler index diagram remained excluded.

## 1. Mathematical corrections and missing conditions

### 1.1 Radius of the R-SCHA2D edge data

In source chapter 05 (Section 12 of the handbook), the correction potential had a circular-edge derivative condition at `θ = θ0`, while the radial factor `R2(r)` varied with radius. Imposing the same edge data at all radii would be inconsistent. The condition is now stated at `r = a, θ = θ0`, as required by the single-surface construction.

### 1.2 Decay at infinity and uniqueness

Surface normal data and northward data on the circular edge, together with decay at infinity, admit multiple correction potentials. For `a = 1` and any `τ > 0`, take

```text
Rτ(r) = r^(−1/2)[cos(τ log r) + sin(τ log r)/(2τ)].
```

Here `Rτ′(1) = 0` and `Rτ(r) → 0`. With a regular Mehler angular function normalized by its edge derivative, different `τ` give the same surface normal and edge data but different interior potentials. At the tested point, `τ = 0.7` and `τ = 1.3` give a difference of about 2.24267. The revised text presents `n = −1/2` as the classical R-SCHA2D choice. [Original boundary formulation](https://academic.oup.com/gji/article/174/1/118/574974)

### 1.3 Information from a single radius

For normal data alone, the relation `z = −(n+1) gi + n ge` leaves two unknowns. With the surface potential coefficient `c = gi + ge`, the determinant is `2n+1 > 0` and both branches can be recovered. Exact, complete tangential and normal fields give Cauchy data up to a potential gauge. Stability with finite noisy observations is a separate question. The Neumann degree set now specifies `n ≥ 0`, a nonzero angular function, and one count for repeated degrees.

### 1.4 Matching both boundaries under truncation

Sections 20–21 and the summary now require both radial branches to retain the two-boundary projected normal data for each angular mode. In a hemispherical cap with `ra = 1` and `rb = 2`, the potential `r³ P3(cos θ)` has outward normal coefficients `(−3, 12)`. The boundary matrix is

```text
[ 4       −3 ]
[ −0.125  12 ]
```

The solution has internal coefficient 0 and external coefficient 1. Omitting the external branch loses the simultaneous match. Asymmetric truncation remains available as a restriction of the full expansion.

### 1.5 Energy on a fixed sphere

`∇V ∈ L²(Ω)` alone does not bound the restricted gradient in `L²(Σβ)` for every fixed sphere. Boundary spheres require a gradient trace. In the unit ball,

```text
V = Σ(n≥1) rⁿ Pn(cos θ)/n
```

has finite volume Dirichlet energy `Σ 4π/[n(2n+1)]`. The outer-sphere potential behaves as `−log θ + O(1)` near the north pole, and the tangential gradient behaves as `−1/θ`. Its squared surface gradient is nonintegrable on any cap containing that pole. Restriction to `ra < r < 1` still has finite volume energy. The fixed-sphere formula now requires an L² gradient trace, and its infinite Gram expansion requires convergence of partial sums in the surface gradient norm. The finite-mode formula still applies.

### 1.6 Finite and infinite Gram matrices

A finite real Gram matrix is symmetric positive semidefinite. Its infinite counterpart requires an explicit quadratic-form domain. The coefficient kernel is the inverse image of the constant potentials. It is one-dimensional if the representation contains exactly one gauge constant and the other modes are independent. Redundancy may enlarge it. Uniqueness for the continuous gradient does not bound the smallest finite eigenvalue uniformly in truncation order.

### 1.7 Weak Neumann and Dirichlet data

The revision distinguishes the existence of boundary integrals from boundedness of the weak Neumann functional. L² densities give a simple sufficient condition, with the general formulation in dual trace spaces. Dirichlet data require an H¹ lifting. The weak uniqueness proof tests with the solution itself, and the classical facewise Green identity is used under the needed regularity. [Trace and elliptic theory reference](https://assets.cambridge.org/97805216/63328/frontmatter/9780521663328_frontmatter.pdf)

### 1.8 Completeness in harmonic energy

Section 13 now explains the passage from component L² completeness to energy convergence.

- For `V1`, expand in the angular Dirichlet basis and apply Parseval to the angular quadratic form and radial derivative. The projected weak Laplace equation gives the two Euler branches.
- For `V2`, expand in radial Neumann and longitude bases and use the radial quadratic form. The weak angular ODE and finite energy select the regular polar branch.
- Summability of the nonnegative energy terms gives vanishing tails in the energy norm.

A Lax–Milgram construction in the full harmonic H¹ space supplies the direct-sum decomposition and the solvability of the subproblems.

### 1.9 Uniqueness of the projection and its coefficients

The finite-dimensional orthogonal projection is unique as a function. Its coefficient vector is unique when the generating family is linearly independent. For nonhomogeneous essential conditions, the Galerkin discretization uses a lifting. The lifted problem defines the comparison with prescribed boundary data.

### 1.10 Injectivity of restriction and rank after projection

Unique continuation supports injectivity of exact function restriction. A finite R-SCHA projection may change the rank, which must be checked separately. Regularization is useful for stable estimation, while mathematical uniqueness of an exact finite system depends on its rank.

### 1.11 Global and facewise normal traces

Normals jump between faces of the cap shell. For `u = z` and `θ0 = π/3`, the outward derivatives at the edge of the inner sphere and the conical side tend to `−1/2` and `−√3/2`, respectively. This jump prevents global H½ regularity even though the potential is smooth. The revised result uses facewise H½ regularity under H² assumptions and distinguishes it from the global weak H⁻½ normal trace supported by harmonic H¹ convergence.

### 1.12 Interior estimates

Shrinking the angular cap can leave points on the top and bottom spheres. The interior estimates in Section 22 are now stated on compactly contained subdomains. The Poincaré estimate uses a common volume-mean gauge.

### 1.13 Branches and scope of the remaining formulas

- A Cartesian domain definition includes the polar axis and distinguishes it from a coordinate singularity.
- A finite product defines the Pochhammer symbol at negative integers.
- The `τ = 0` normalization limit is distinguished from the positive radial spectrum of a nondegenerate closed shell.
- The regular ODE pole series gives the derivative bounds beyond the function estimate `A = O(θᵐ)`.
- Thin-shell asymptotics are stated for a fixed angular degree.
- The noninteger radial branches have internal/external labels with a radial meaning.
- Active and passive rotations specify which basis changes. The separation constant is `Λ` and longitude is `λ`.

Special-function conventions are based on [NIST DLMF Chapter 14](https://dlmf.nist.gov/14).

## 2. Claims and repetition addressed

1. The introduction states its scope and points to the relevant sections.
2. The conclusion summarizes scope, assumptions, and remaining topics.
3. Topics outside the handbook are identified as such. The finite real symmetric Gram matrix already has an orthogonal diagonalization.
4. Relative convergence claims about NN, DD, and DN require comparable data, truncations, and error norms.
5. The edge-buffer discussion now gives the Cauchy–Schwarz bound `|Eβ[V]−Eβ[V1]| ≤ Eβ[V2] + 2√(Eβ[V1]Eβ[V2])`. The relative bound `2ε + ε²` uses `Eβ[V2]/Eβ[V1] ≤ ε²`.
6. Edge identities are stated for finite expansions. Claims about Gibbs error concentration and infinite H¹ limits require additional assumptions.
7. Repeated energy notation, Green-identity steps, single-mode factor derivations, and repeated summary items were consolidated.
8. The previous removal of the Mehler indexing debate and dedicated diagram was retained.

## 3. Bibliographic checks

- Torta's article is cited as 2020, volume 41(2), pages 201–247, with the 2019 online publication date. The citation key is retained. [Publisher record](https://link.springer.com/article/10.1007/s10712-019-09576-2)
- Thébault and Gaya-Piqué (2008) now has issue 7, article Q07005, and the DOI. [Publisher record](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008GC001953)
- Verified identifiers were added for two De Santis papers and Talarn's paper, alongside the weak-solution and trace-theory references. [1996 rotation paper](https://academic.oup.com/gji/article/126/1/263/607834)

These checks cover the cited metadata and conventions. The handbook's derivations are stated under their own assumptions.

## 4. Verification results

| Check | Scope | Outcome |
|---|---|---|
| Original formulas | 450 comparisons, 45-digit working precision | Passed, maximum scaled error about `1.52e-13` |
| Boundary and special functions | 144 comparisons, including small τ, wider caps, thin and thick shells, and the R-SCHA2D counterexample | Passed, maximum scaled error about `4.78e-13` |
| Energy, tensor, and flux | 250 comparisons at 45 digits, seven Cartesian harmonics and several cap and longitude settings | Passed, maximum scaled error about `2.98e-43` |
| Wigner rotation | Degrees 0–6, three Euler angle sets, several points, and 588 defining-identity comparisons | Maximum absolute error about `2.01e-15`, with unitarity and coefficient transpose checks passing |
| Structure | 1,177 labels | Balanced environments, with all references and citation keys resolved and no duplicate labels |
| Build and layout | 160 pages, 26 key pages rendered for detailed review | Compiled without unresolved references, overfull or underfull boxes, or missing-glyph warnings |

Each check uses its own error metric and scale.

## 5. Scope of the conclusions

The checks support the stated special-function identities, radial Neumann spectrum and normalization, finite-mode volume and surface energy blocks, NED gradient and Hessian, rotation conventions, and closed-shell flux. The energy-completeness proof follows the stated quadratic forms and convergence conditions.

The handbook states its boundary, regularity, gauge, and convergence conditions. Uniform high-degree stability, term-by-term differentiation of arbitrary infinite expansions, finite energy for arbitrary boundary data, high-order regularity at edges, stable source separation from finite noisy data, and unresolved questions in the wider literature require their own evidence.
