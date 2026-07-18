"""Manifest builders: good-run CSV (one row per PNG) and bad-logs CSV."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .model import Geometry, Mode

MANIFEST_COLUMNS = [
    "source_log", "png", "mode_index", "frequency_cm-1", "ir_intensity_km/mol",
    "reduced_mass_amu", "is_imaginary", "n_atoms", "geometry_frame", "rendered",
]

BAD_LOG_COLUMNS = ["source_log", "reason"]


def manifest_row(source: Path, geometry: Geometry, mode: Mode,
                 png: Path | None, rendered: bool) -> dict:
    """Build a single manifest row dict."""
    return {
        "source_log": Path(source).name,
        "png": str(png) if png is not None else "",
        "mode_index": mode.index,
        "frequency_cm-1": mode.frequency,
        "ir_intensity_km/mol": mode.ir_intensity,
        "reduced_mass_amu": mode.reduced_mass,
        "is_imaginary": mode.is_imaginary,
        "n_atoms": geometry.n_atoms,
        "geometry_frame": geometry.source_frame,
        "rendered": rendered,
    }


def build_manifest(rows: list[dict]) -> pd.DataFrame:
    """Assemble manifest rows into a DataFrame with a stable column order."""
    return pd.DataFrame(rows, columns=MANIFEST_COLUMNS)


def build_bad_logs(rows: list[dict]) -> pd.DataFrame:
    """Assemble bad-log rows (source_log + reason) into a DataFrame."""
    return pd.DataFrame(rows, columns=BAD_LOG_COLUMNS)
