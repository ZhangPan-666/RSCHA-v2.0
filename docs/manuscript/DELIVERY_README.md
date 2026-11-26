# R-SCHA manuscript delivery notes

Updated 2026-10-24.

## Delivery at this revision

- Manuscript revision `manuscript-20261024-r2`. The software source directories were v2.0 and the release was v1.0 at delivery.
- Manuscript `Zhang et al., 2026, Manuscript_RSCHAFunctions`, with 38 pages, five figures, and 50 references.
- Source `Zhang et al., 2026, Manuscript_RSCHAFunctions.tex`, bibliography `references.bib`, and generated `.bbl`.
- The checked PDF is `Zhang et al., 2026, Manuscript_RSCHAFunctions - manuscript-20261024-r2.pdf`. The revised PDF and source ZIP identify the checked copy.
- Submission highlights are in `JCP_Highlights.txt`.

## Directory guide

- `figures/` contains publication figures and high-resolution previews.
- `reproducibility/` contains numerical requests, reference-generation scripts, spectral data, inversion parameters, and the regional residual audit.
- `submission/` contains the JCP cover letter draft, declarations, and preparation notes.
- `Archive/` was the working directory for pre-edit backups.

## Second revision, 2026-10-24

- Inserted the approved geometry and algorithm figure in Section 2 as Figure 1. The original four figures retain their filenames and appear as Figures 2–5.
- Packaged the plotting code, requirements, PDF/PNG/SVG, caption, and checksums in `Geometry_Algorithm_Figure/`. Regenerating the PNG matched the approved artwork pixel for pixel.
- Added local coordinates, shell and side boundaries, and the algorithm overview. The relevant algorithm paragraphs refer to the corresponding Figure 1 panels.
- Used `lambda_d` for inversion damping and `tau_stop` for the stopping tolerance. Specified the linear projection used for the total-field anomaly.
- Changed the short-wavelength statement to refer to representability and kept independent resolution validation as a separate requirement. Corrected the final Zhang et al. (2025) page range to 686–696.
- Checked labels, cross-references, 50 references, figures, tables, and principal statistics. The six Python reference regenerations and 24 original controls passed. The report is in `reproducibility/editorial_audit_20261024/`. The earlier C++ accuracy, timing, and regional inversion runs remain historical results.
- The AI-use declaration was deferred at the author's request. The draft remains in the submission material.
- The revised PDF and standalone source ZIP were published under `docs/manuscript/`.

## First revision, 2026-10-24

- Retained **Stable Computation of Revised Spherical Cap Harmonics** and changed the two relevant “error measures” headings to “error evaluation.”
- Listed National Natural Science Foundation of China grants 42574077, 42250103, and 42674085 and Natural Science Foundation of Hubei Province grant 2024AFE003. The laboratory open-fund application was pending at this revision.
- Listed Pan Zhang's affiliations 1, 2, and 3, including the school and key laboratory at China University of Geosciences.
- Updated the submission declarations. The figures, numerical sources, and historical results were carried forward.
- Rebuilt 37 pages and 50 references. The title page, revised headings, and funding pages were checked. The log had no unresolved references or overfull boxes.
- Published the first revision as `manuscript-20261024`.

## Software and downloads

At delivery, the R-SCHA distribution provided OpenMP and MPI sources, lookup tables, build scripts, validation summaries, and the mathematical handbook. [Project repository](https://github.com/ZhangPan-666/RSCHA-v2.0) · [Release history](https://github.com/ZhangPan-666/RSCHA-v2.0/releases)

## Build

Run these commands from the extracted manuscript source directory:

```powershell
pdflatex "Zhang et al., 2026, Manuscript_RSCHAFunctions.tex"
bibtex "Zhang et al., 2026, Manuscript_RSCHAFunctions"
pdflatex "Zhang et al., 2026, Manuscript_RSCHAFunctions.tex"
pdflatex "Zhang et al., 2026, Manuscript_RSCHAFunctions.tex"
```

## Synchronizing a revision

1. Make the changes locally and check the PDF and relevant numerical outputs.
2. Update the revision identifier, change notes, and source package.
3. Commit the revised source and artifacts to GitHub. Use a new release or tag for substantive software releases.
4. Use consistent links and versions in **Data and code availability**.
5. Retain the manuscript PDF, sources, figures, and reproduction files with the revision.

Before submitting, verify names, affiliations, email addresses, ORCIDs, corresponding authorship, contributions, funding, conflicts of interest, all authors' approval, and the submission system's file requirements.
