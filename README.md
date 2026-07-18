# 3D Molecular Vibration PNG Generator

Fork of [ashendeema/molecular-vibration-png-generator](https://github.com/ashendeema/molecular-vibration-png-generator), maintained at [ishgits/molecular-vibrations](https://github.com/ishgits/molecular-vibrations). This fork adds a `molvib` package for Gaussian 16 frequency-log parsing, IR-intensity mode selection, and CSV manifest output.

A Python tool for generating **publication-quality 3D molecular vibration images** from XYZ files containing atomic coordinates and vibrational displacement vectors.

Developed by **Ashen Deemantha Liyanage**  
**Zayak's Lab**  
Department of Physics and Astronomy  
Bowling Green State University (BGSU)

---

## Overview

Visualizing molecular vibrational modes is an essential part of computational chemistry and materials science. This program automatically converts molecular vibration data stored in XYZ files into high-resolution PNG images suitable for:

- Research publications
- Conference presentations
- Teaching materials
- Molecular vibration analysis
- Computational chemistry visualization

The script automatically detects all XYZ files in the working directory, renders the molecular structure in 3D, draws vibrational displacement vectors, and saves a publication-quality PNG image for each molecule.

---

## Features

- Automatically detects all `.xyz` files in the current folder
- Generates high-resolution PNG images
- Publication-quality 3D rendering using PyVista
- Automatic bond detection using covalent radii
- Element-specific atom colors
- Vibrational displacement vectors displayed as red arrows
- Automatic extraction of vibrational frequency from the filename
- Transparent or colored backgrounds
- Adjustable rendering resolution
- Customizable atom, bond, and arrow sizes
- Suitable for DFT vibrational analysis

---

## Example

Input file:

```
m_165_3150.62.xyz
```

Output:

```
m_165_3150.62.png
```

The generated image contains:

- Colored atoms
- Chemical bonds
- Vibrational displacement vectors
- Vibrational frequency label

---

## Input File Format

Each XYZ file should contain one atom per line in the following format:

```
Element   X   Y   Z   Fx   Fy   Fz
```

Example:

```
C   0.000   0.000   0.000   0.12   0.03   0.00
O   1.210   0.000   0.000  -0.10   0.02   0.01
H  -0.620   0.930   0.000   0.03  -0.05   0.00
```

where

- **X Y Z** = Atomic coordinates (Å)
- **Fx Fy Fz** = Vibrational displacement vectors

### Benzene

![Benzene](images/m_20_1023.39.png)

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/ishgits/molecular-vibrations.git
```

### 2. Move into the project

```bash
cd molecular-vibrations
```

### 3. Install the required Python libraries

```bash
python -m pip install -e '.[test]'
```

---

## Required Libraries

- Python 3.9 or newer
- NumPy
- PyVista
- VTK

---

## How to Use

Place one or more `.xyz` files in the same directory as the script.

Run

```bash
python xyz_to_png.py
```

The program automatically:

1. Finds all XYZ files
2. Reads atomic coordinates
3. Detects chemical bonds
4. Draws atoms as spheres
5. Draws bonds as cylinders
6. Draws vibrational displacement vectors
7. Extracts vibrational frequency from the filename
8. Saves a PNG image

---

## Gaussian 16 log support (`molvib`)

In addition to `.xyz` files, the `molvib` package reads **Gaussian 16 frequency
output** (`.log` / `.out`) directly, selects modes by **IR intensity**, renders a
PNG per selected mode, and writes a **pandas CSV manifest**. The original
`xyz_to_png.py` script is unchanged and continues to work.

### Command line

```bash
molvib --input examples/logs/adenine_pcm_water.log
# equivalent:
python -m molvib --input examples/logs/adenine_pcm_water.log
```

Two example Gaussian logs ship in `examples/logs/`: `adenine_pcm_water.log` (a
converged minimum) and `neg_test_H2SO4_V_freq.log` (a transition state with one
imaginary mode). See `end_user_test.ipynb` for a runnable walkthrough of both.

Common options:

| Flag | Meaning |
|---|---|
| `--input` | file or directory (globs `*.xyz`, `*.log`, `*.out`) |
| `--ir-threshold` | IR intensity cutoff in KM/Mole (default `10.0`) |
| `--include-imaginary` / `--no-include-imaginary` | always keep imaginary modes (default: include) |
| `--modes 1,27,39` | explicit mode indices; overrides the IR filter |
| `--outdir` | output directory (default: alongside input) |
| `--manifest` | manifest CSV path (default: `outdir/vibration_manifest.csv`) |
| `--resolution 2000x2000` | output image size |
| `--background transparent` | `transparent` or a color |
| `--zoom 0.9` | camera zoom; `<1` adds margin so atoms and arrows are never cropped (default `0.9`) |
| `--no-arrows` | hide displacement arrows |
| `--no-render` | parse + select + write manifest without rendering PNGs |

Output PNGs are named `{logstem}_mode{index:03d}_{freq:.1f}cm.png` to avoid
collisions across many modes. XYZ inputs keep the original `foo.xyz → foo.png`
naming. Logs that fail validation are routed to a `vibration_bad_logs.csv`
alongside the manifest rather than aborting the run.

### Gaussian 16 parsing policy

`molvib` supports standard low-precision **Gaussian 16** frequency output only;
HPModes/high-precision coordinates are not supported. For Link1 or concatenated
files, it parses the **last complete frequency job**. An incomplete final job
falls back to the preceding complete job. Every frequency block must contain an
`Atom AN` header and exactly one complete displacement row per atom; malformed
or truncated blocks are rejected and never rendered. Atomic jobs are rejected.
Linear molecules are accepted when their geometry is collinear and they contain
the expected `3N-5` modes; nonlinear molecules must contain `3N-6` modes.

### Importable API

```python
from molvib import read_any, select_modes
from molvib.cli import render_file, render_directory

# One DataFrame per file, or (manifest_df, bad_logs_df) for a directory:
manifest_df = render_file("adenine_pcm_water.log", ir_threshold=10.0)
manifest_df, bad_df = render_directory("logs/", ir_threshold=10.0)
```

### A note on displacement semantics

Gaussian's printed normal coordinates are **normalized Cartesian displacements**,
not mass-weighted eigenvectors and not forces. (The original tool's `forces`
variable name is a misnomer; in `molvib` these are called `displacements`.) The
renderer normalizes them to the largest displacement purely for visualization —
**arrow lengths are not physical amplitudes.** Geometry is taken from the *last
`Standard orientation` table before the frequency section*, the frame in which
the normal modes are reported, so arrows point along the correct bonds.

---

## Adjustable Settings

The following parameters can easily be modified inside the script.

### Rendering

```python
OUTPUT_RESOLUTION = (2000, 2000)
BACKGROUND = "transparent"
```

### Atom Size

```python
ATOM_RADIUS = 0.28
```

### Bond Size

```python
BOND_RADIUS = 0.10
```

### Arrow Size

```python
ARROW_SHAFT_RADIUS = 0.05
ARROW_TIP_RADIUS = 0.10
ARROW_TIP_LENGTH = 0.25
```

### Show or Hide Force Vectors

```python
SHOW_ARROWS = True
```

---

## Scientific Applications

This tool can be used for

- Density Functional Theory (DFT)
- Molecular vibrational analysis
- Raman spectroscopy visualization
- Infrared spectroscopy
- Computational chemistry
- Materials science
- Molecular dynamics visualization

---

## Output

For every input file

```
molecule.xyz
```

the program automatically generates

```
molecule.png
```

at publication quality.

---

## Repository Structure

```
molecular-vibrations/
│
├── molvib/                     # Gaussian-log + XYZ rendering package
│   ├── readers.py              # read_xyz, read_gaussian_log, read_any
│   ├── select.py               # IR-intensity / imaginary-mode selection
│   ├── render.py               # PyVista rendering + RenderSettings
│   ├── manifest.py             # pandas CSV manifest
│   ├── elements.py             # Z→symbol, colors, covalent radii
│   ├── model.py                # Geometry / Mode / RenderJob dataclasses
│   └── cli.py                  # argparse driver + render_file/render_directory
│
├── xyz_to_png.py               # original standalone XYZ script (unchanged)
├── end_user_test.ipynb         # positive + negative walkthrough / smoke test
│
├── examples/
│   ├── logs/                   # example Gaussian 16 frequency logs
│   └── *.xyz
├── images/                     # reference/example PNGs
├── tests/
│   ├── test_readers.py
│   └── data/                   # parser fixtures
│
├── pyproject.toml
├── requirements.txt
├── LICENSE
└── README.md
```

---

## Future Improvements

Planned features include

- Multiple lighting models
- Custom camera angles
- Automatic legends
- Batch rendering options
- Surface visualization
- High-resolution publication presets

---

## Citation

If this program contributes to your research, please cite this GitHub repository.

---

## License

This project is distributed under the MIT License.

---

## Authors

**Original tool** — **Ashen Deemantha Liyanage**  
PhD Student, Department of Physics and Astronomy  
Bowling Green State University (BGSU) · Zayak's Lab  
GitHub: https://github.com/ashendeema

**`molvib` Gaussian-log extension** — **Ish** ([@ishgits](https://github.com/ishgits))  
Gaussian 16 frequency-log parsing, IR-intensity mode selection, CSV manifest, and CLI.

---

If you find this project useful, consider giving it a ⭐ on GitHub.

# molecular-vibration-png-generator
Python tool for generating publication-quality 3D molecular vibration images from XYZ files using PyVista.
