# Manuscript editorial and reproducibility review

Revision `manuscript-20261024-r2`, reviewed on 2026-10-24.

**Stable Computation of Revised Spherical Cap Harmonics** contains 38 pages, five figures, 12 tables, and 50 references. The review covered the text, equations, cross-references, figures, statistics, bibliography, reproduction material, and PDF layout. Numerical results refer to the saved runs. The full C++ acceptance, throughput, and regional inversion runs were outside this review.

## Changes

| Item | Resolution |
|---|---|
| Section 2 geometry | Added the approved four-panel Figure 1 on page 3 and explained `O`, `z_g`, the cap axis, `P`, `r1/r2`, `S1/S2`, and `theta/theta0/phi`. |
| Algorithm overview | Linked Section 4's root, norm, and residual-compression descriptions to Figure 1. The caption identifies the output scales and Mehler normalization identity. The runtime does not integrate a Riccati ODE. |
| Figure numbering | Let LaTeX renumber the original figures as 2–5. Retained filenames and linked the manuscript and reproduction documentation to the figure source. |
| Notation | Used `lambda_d` for damping and `tau_stop` for the BiCGSTAB stopping tolerance to distinguish them from spectral parameters. |
| Total-field operator | Described the linear projection of the anomalous field onto the prescribed background-field direction. |
| Spatial-resolution statement | Described representable shorter wavelengths. Independent observations are needed to determine which wavelengths are resolved. |
| Publication metadata | Corrected Zhang et al. (2025), DOI 10.26464/epp2025057, to the final page range 686–696. |
| Availability | Linked the manuscript revision and figure code from `docs/manuscript/`. |
| Geometry package | Included executable drawing code, requirements, caption, manifest, and vector and raster images. The publication figure is a byte-identical PDF copy. |

The existing title, author order, confirmed grant details, first-author affiliations, four data figures, numerical source, and historical acceptance findings were carried forward.

## Verification

- Checked 108 labels and 81 references with no duplicate labels or missing targets. Algorithm references use labels.
- All 50 bibliography records were cited. BibTeX completed without warnings.
- The five figures, 12 tables, captions, and references were checked. Figure 1 has embedded fonts and no raster layers. Regeneration matched the approved PNG pixel for pixel.
- The LaTeX and BibTeX build completed with all references resolved and no overfull or underfull warnings. All 38 pages were rendered and reviewed, with a closer check of Figure 1.
- Eighteen statistics and artwork checks matched historical `acceptance.json` within printed rounding. The earlier result remains tied to the tested source snapshot.
- Six archived high-precision requests were regenerated using Python, and the original 24 controls passed. This checks the reference script. The historical 40,995-request C++ comparison was a separate run.
- All 39 DOI records were queried. Thirty-two succeeded, and seven interface failures were logged. The publisher or institutional record was checked separately for key references. Blakely et al. (2005) retains pages 445–448.
- The five highlights contain 60–71 characters. Each is within the 85-character limit.
- The MHGA implementation's clipping matches the published saturation function. Its pseudo-vertical derivative and `Hz/den−pi/3` agree with the manuscript equation. The text describes the use of `Bz` and the processing differences involving RTP and MNLM.

## Items for submission

1. The author deferred the AI-use declaration. The source ZIP includes the draft for author review.
2. Complete regional inversion logs and execution manifests are needed to establish actual iteration counts and stopping reasons. The manuscript records the 100-iteration limit. The residuals describe fit.
3. Spectral power and geological correspondence support interpretation. Withheld observations are required to test independent spatial resolution.
4. The analytical review checked internal consistency and saved numerical comparisons. Uniform error bounds over continuous parameters and a full review of all 50 references are outside its scope.
5. The authors should verify current JCP form fields, their consent, and final funding and affiliation entries at submission.

## Sources checked

- [Final issue for Zhang et al. (2025)](https://www.eppcgs.org/article/2025/3)
- [MHGA article and equations (10)–(11)](https://journal.geo.sav.sk/cgg/article/view/514)
- [Gil et al. (2012) institutional record](https://ir.cwi.nl/pub/18881) and [author publication list](https://personales.unican.es/segurajj/publirk.html)
- [Elsevier AI-use policy](https://www.elsevier.com/about/policies-and-standards/generative-ai-policies-for-journals)

## Evidence files

- `structure_audit.json`, `statistics_audit.json`, and `bibliography_registry_audit.json` record the structure, 18 statistics comparisons, and 39 DOI queries.
- `reference_regeneration_summary.json` and `reference_self_check.json` record the six requests and 24 controls.
- `build_qa.json` records the document counts and log checks.
- `MANUSCRIPT_SHA256.json` lists the source-package checksums. The earlier release and numerical acceptance records retain their original scope.
