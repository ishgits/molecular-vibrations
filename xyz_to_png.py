#!/usr/bin/env python3
"""
===============================================================================
3D Molecular Vibration PNG Generator

Author:
    Ashen Deemantha Liyanage

Affiliation:
    Zayak's Lab
    Department of Physics and Astronomy
    Bowling Green State University (BGSU)

GitHub:
    https://github.com/ashendeema

Description:
    A Python tool for generating publication-quality 3D molecular vibration
    images from XYZ files containing atomic coordinates and vibrational
    displacement vectors.

Features:
    • Automatic detection of XYZ files
    • High-resolution PNG rendering
    • Element-colored atoms
    • Automatic bond detection using covalent radii
    • Vibrational displacement vectors (red arrows)
    • Automatic vibrational frequency extraction from filenames
    • Transparent or colored backgrounds
    • Publication-ready rendering using PyVista

Input Format:
    Element   X   Y   Z   Fx   Fy   Fz

where

    X, Y, Z  = Atomic coordinates (Å)
    Fx, Fy, Fz = Vibrational displacement vectors

Output:
    molecule.xyz
        ↓
    molecule.png

Requirements:
    numpy
    pyvista
    vtk

Developed for:
    Visualization of molecular vibrational normal modes obtained from
    Density Functional Theory (DFT) calculations.

Version:
    1.0.0

License:
    MIT License
===============================================================================
"""

import glob
import os
import re
import numpy as np
import pyvista as pv


# ==============================================================================
# SETTINGS
# ==============================================================================

OUTPUT_RESOLUTION = (2000, 2000)
BACKGROUND = "transparent"

ATOM_RADIUS = 0.28
BOND_RADIUS = 0.10

ARROW_SHAFT_RADIUS = 0.05
ARROW_TIP_RADIUS = 0.10
ARROW_TIP_LENGTH = 0.25

SHOW_ARROWS = True


# ==============================================================================
# ELEMENT COLORS
# ==============================================================================

def get_element_color(sym):
    colors = {
        "O": "red",
        "Ba": "green",
        "Ti": "gray",
        "C": "black",
        "N": "blue",
        "Cu": "orange",
        "H": "white",
    }
    return colors.get(sym, "lightgray")


def get_covalent_radius(sym):
    radii = {
        "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66,
        "Ti": 1.36, "Ba": 1.61, "Cu": 1.32,
    }
    return radii.get(sym, 0.8)


# ==============================================================================
# READ XYZ
# ==============================================================================

def read_xyz(filename):
    symbols, positions, forces = [], [], []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) < 4:
                continue

            symbols.append(parts[0])
            positions.append(list(map(float, parts[1:4])))

            if len(parts) >= 7:
                forces.append(list(map(float, parts[4:7])))
            else:
                forces.append([0.0, 0.0, 0.0])

    return symbols, np.array(positions), np.array(forces)


# ==============================================================================
# BOND DETECTION
# ==============================================================================

def get_bonds(symbols, positions, scale=1.2):
    bonds = []
    radii = np.array([get_covalent_radius(s) for s in symbols])

    for i in range(len(symbols)):
        for j in range(i + 1, len(symbols)):
            cutoff = scale * (radii[i] + radii[j])
            dist = np.linalg.norm(positions[i] - positions[j])

            if dist <= cutoff:
                bonds.append((i, j))

    return bonds


# ==============================================================================
# EXTRACT FREQUENCY FROM FILENAME
# ==============================================================================

def extract_frequency(filename):
    """
    Example:
        m_165_3150.62.xyz → 3150.62
    """
    name = os.path.splitext(filename)[0]
    nums = re.findall(r"\d+\.\d+|\d+", name)

    if nums:
        return nums[-1]   # last number = frequency
    return None


# ==============================================================================
# MAIN FUNCTION
# ==============================================================================

def make_png(xyz_file):

    symbols, positions, forces = read_xyz(xyz_file)

    # normalize forces
    norms = np.linalg.norm(forces, axis=1)
    maxn = norms.max() if norms.max() != 0 else 1.0
    forces = forces / maxn

    bonds = get_bonds(symbols, positions)

    plot = pv.Plotter(off_screen=True, window_size=OUTPUT_RESOLUTION)

    # background
    if BACKGROUND == "transparent":
        plot.set_background("white")
        transparent = True
    else:
        plot.set_background(BACKGROUND)
        transparent = False

    # ---------- ATOMS ----------
    for sym, pos in zip(symbols, positions):
        sphere = pv.Sphere(radius=ATOM_RADIUS, center=pos,
                           theta_resolution=64, phi_resolution=64)

        plot.add_mesh(
            sphere,
            color=get_element_color(sym),
            smooth_shading=True,
            ambient=0.3,
            diffuse=0.7,
            specular=0.6,
            specular_power=40
        )

    # ---------- BONDS ----------
    for i1, i2 in bonds:
        p1, p2 = positions[i1], positions[i2]
        d = p2 - p1
        L = np.linalg.norm(d)

        if L > 0:
            cyl = pv.Cylinder(center=(p1+p2)/2, direction=d,
                              height=L, radius=BOND_RADIUS, resolution=60)

            plot.add_mesh(cyl, color="gray", smooth_shading=True)

    # ---------- ARROWS ----------
    if SHOW_ARROWS:
        for pos, vec in zip(positions, forces):
            if np.linalg.norm(vec) == 0:
                continue

            arrow = pv.Arrow(start=pos, direction=vec,
                             tip_length=ARROW_TIP_LENGTH,
                             tip_radius=ARROW_TIP_RADIUS,
                             shaft_radius=ARROW_SHAFT_RADIUS)

            plot.add_mesh(arrow, color="red", opacity=0.9)

    # ---------- CAMERA ----------
    plot.reset_camera()
    plot.camera_position = "iso"
    plot.camera.zoom(1.5)

    # ---------- ADD FREQUENCY LABEL ----------
    freq = extract_frequency(xyz_file)
    if freq:
        label = f"{freq} 1/cm"
        plot.add_text(label, position="upper_left",
                      font_size=28, color="black")

    # render before saving
    plot.show(auto_close=False)

    # ---------- SAVE ----------
    out = xyz_file.replace(".xyz", ".png")
    plot.screenshot(out, transparent_background=transparent)

    plot.close()

    print(f"[OK] Saved: {out}")


# ==============================================================================
# RUN
# ==============================================================================

if __name__ == "__main__":
    xyz_files = sorted(glob.glob("*.xyz"))

    if not xyz_files:
        print("No .xyz files found.")
    else:
        for f in xyz_files:
            make_png(f)

