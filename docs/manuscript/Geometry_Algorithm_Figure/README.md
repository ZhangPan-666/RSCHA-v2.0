# R-SCHA geometry and algorithm figure

This directory contains the approved four-panel figure, its plotting code, and caption. The geometry is rotated 12 degrees clockwise about the geographic axis. The outer sphere uses a sparse dashed white grid with 1.2 pt lines at 0.80 opacity.

## Files

- `draw_geometry_algorithm.py` reproduces the figure without external observational data.
- `images/RSCHA_Geometry_Algorithm.pdf` is the vector manuscript figure with embedded fonts.
- `images/RSCHA_Geometry_Algorithm.svg` is the editable vector version.
- `images/RSCHA_Geometry_Algorithm.png` is the complete preview. The `_geometry.png` file previews the geometry panel.
- `caption.tex` contains the LaTeX caption.
- `requirements.txt` and `manifest.json` record the tested dependencies and file checksums.

## Run

Use Python 3.12 or newer from this directory:

```bash
python -m pip install -r requirements.txt
python draw_geometry_algorithm.py
```

The default destination is `images/`. To choose another directory, supply it as the first argument. PDF timestamps may cause checksum differences between runs. The geometry and layout parameters are fixed by the script.

## Manuscript reference

The figure appears in Section 2 as Figure 1 with the label `fig:geometry_algorithm`. The manuscript copy is `figures/Fig1_Geometry_Algorithm.pdf` relative to the extracted source root. The other four figures retain their filenames and appear as Figures 2–5.

Panel (b) shows stable evaluation for the three function families. Panel (c) separates offline residual compression from online recovery. Panel (d) connects Legendre degree sensitivity, boundary roots, and cap normalization. The Riccati label on the Mehler panel refers to the normalization identity. The runtime evaluates the normalization without integrating a Riccati ODE. The drawings are schematic.
