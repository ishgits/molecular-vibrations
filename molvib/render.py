"""Rendering: lifted from the original ``xyz_to_png.make_png()``.

The geometry/mode/camera pipeline (Sphere/Cylinder/Arrow, covalent-radius bond
detection, iso camera) is unchanged from the original tool. Only the inputs
(displacements + label from parsed data) and the output naming differ.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .elements import get_covalent_radius, get_element_color
from .model import RenderJob


@dataclass
class RenderSettings:
    """Rendering parameters, defaults matching the original xyz_to_png constants."""

    resolution: tuple[int, int] = (2000, 2000)
    background: str = "transparent"      # "transparent" | color name/hex
    atom_radius: float = 0.28
    bond_radius: float = 0.10
    arrow_shaft_radius: float = 0.05
    arrow_tip_radius: float = 0.10
    arrow_tip_length: float = 0.25
    show_arrows: bool = True
    bond_scale: float = 1.2
    camera_zoom: float = 0.9   # <1 adds margin so atoms+arrows aren't cropped


def get_bonds(symbols: list[str], positions: np.ndarray, scale: float = 1.2):
    """Bond list by covalent-radius cutoff (unchanged from original)."""
    bonds = []
    radii = np.array([get_covalent_radius(s) for s in symbols])
    for i in range(len(symbols)):
        for j in range(i + 1, len(symbols)):
            cutoff = scale * (radii[i] + radii[j])
            if np.linalg.norm(positions[i] - positions[j]) <= cutoff:
                bonds.append((i, j))
    return bonds


def default_output_name(job: RenderJob) -> str:
    """Output filename for a job.

    - XYZ: preserve original ``foo.xyz -> foo.png`` naming.
    - Log: ``{logstem}_mode{index:03d}_{freq:.1f}cm.png`` to avoid collisions
      across many modes.
    """
    src = job.source_file
    if src.suffix.lower() == ".xyz":
        return src.stem + ".png"
    return f"{src.stem}_mode{job.mode.index:03d}_{job.mode.frequency:.1f}cm.png"


def render_mode(job: RenderJob, settings: RenderSettings | None = None,
                outdir: str | Path | None = None) -> Path:
    """Render one RenderJob to a PNG and return its path."""
    import pyvista as pv  # imported lazily so parsing/tests don't require pyvista

    settings = settings or RenderSettings()

    symbols = job.geometry.symbols
    positions = job.geometry.positions
    disps = np.asarray(job.mode.displacements, dtype=float)

    # Normalize displacements to the largest (visualization only).
    norms = np.linalg.norm(disps, axis=1)
    maxn = norms.max() if norms.size and norms.max() != 0 else 1.0
    disps = disps / maxn

    bonds = get_bonds(symbols, positions, scale=settings.bond_scale)

    plot = pv.Plotter(off_screen=True, window_size=settings.resolution)

    if settings.background == "transparent":
        plot.set_background("white")
        transparent = True
    else:
        plot.set_background(settings.background)
        transparent = False

    # ---------- ATOMS ----------
    for sym, pos in zip(symbols, positions):
        sphere = pv.Sphere(radius=settings.atom_radius, center=pos,
                           theta_resolution=64, phi_resolution=64)
        plot.add_mesh(sphere, color=get_element_color(sym), smooth_shading=True,
                      ambient=0.3, diffuse=0.7, specular=0.6, specular_power=40)

    # ---------- BONDS ----------
    for i1, i2 in bonds:
        p1, p2 = positions[i1], positions[i2]
        d = p2 - p1
        L = float(np.linalg.norm(d))
        if L > 0:
            cyl = pv.Cylinder(center=(p1 + p2) / 2, direction=d,
                              height=L, radius=settings.bond_radius, resolution=60)
            plot.add_mesh(cyl, color="gray", smooth_shading=True)

    # ---------- ARROWS ----------
    if settings.show_arrows:
        for pos, vec in zip(positions, disps):
            if np.linalg.norm(vec) == 0:
                continue
            arrow = pv.Arrow(start=pos, direction=vec,
                             tip_length=settings.arrow_tip_length,
                             tip_radius=settings.arrow_tip_radius,
                             shaft_radius=settings.arrow_shaft_radius)
            plot.add_mesh(arrow, color="red", opacity=0.9)

    # ---------- CAMERA ----------
    # Orient first, THEN fit to all actors (atoms + arrows), then apply zoom.
    # The original zoomed 1.5x *after* framing, which cropped the molecule.
    plot.camera_position = "iso"
    plot.reset_camera()
    plot.camera.zoom(settings.camera_zoom)

    # ---------- LABEL ----------
    if job.label:
        plot.add_text(job.label, position="upper_left", font_size=28, color="black")

    plot.show(auto_close=False)

    out_dir = Path(outdir) if outdir is not None else job.source_file.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / default_output_name(job)
    plot.screenshot(str(out), transparent_background=transparent)
    plot.close()
    return out
