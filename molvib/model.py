"""Core data model shared across readers, selection, and rendering."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class Geometry:
    """A single molecular geometry (one shared frame per file)."""

    symbols: list[str]          # element symbols, len == n_atoms
    positions: np.ndarray       # (n_atoms, 3), Angstrom
    source_frame: str           # "standard" | "input" | "xyz"

    @property
    def n_atoms(self) -> int:
        return len(self.symbols)


@dataclass
class Mode:
    """A single vibrational normal mode."""

    index: int                  # 1-based Gaussian mode number
    frequency: float            # cm^-1 (negative => imaginary)
    ir_intensity: float | None  # KM/Mole; None if not present
    reduced_mass: float | None
    force_const: float | None
    displacements: np.ndarray   # (n_atoms, 3), normalized Cartesian normal coords
    is_imaginary: bool          # frequency < 0

    @property
    def label(self) -> str:
        return f"{self.frequency:.1f} cm^-1 (mode {self.index})"


@dataclass
class RenderJob:
    """One geometry + one mode = one PNG to render."""

    geometry: Geometry
    mode: Mode
    source_file: Path
    label: str
