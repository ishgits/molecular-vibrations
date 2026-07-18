"""Readers for XYZ and Gaussian 16 log/out frequency files.

All readers return ``(Geometry, list[Mode])`` so the rest of the pipeline is
format-agnostic. ``read_any`` dispatches on file suffix.

Displacement semantics
----------------------
Gaussian's printed normal coordinates are *normalized Cartesian displacements*,
not mass-weighted eigenvectors and not forces. We store them as ``displacements``
(the original tool's ``forces`` name was a misnomer). Arrow lengths are for
visualization only — the renderer normalizes them to the largest displacement —
and are not physical amplitudes.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np

from .elements import symbol_from_z
from .model import Geometry, Mode

# --------------------------------------------------------------------------- #
# Regex helpers (mirror the notebook's regex-helper block)
# --------------------------------------------------------------------------- #
RE_NORMAL_TERM = re.compile(r"Normal termination of Gaussian")
RE_HARMONIC = re.compile(r"Harmonic frequencies \(cm\*\*-1\)")
RE_STD_ORIENT = re.compile(r"Standard orientation:")
RE_INPUT_ORIENT = re.compile(r"Input orientation:")
RE_LINK1 = re.compile(r"(?:--)?Link1(?:--)?", re.IGNORECASE)
RE_FREQ_LINE = re.compile(r"Frequencies --\s+(.*)")
RE_RED_MASS = re.compile(r"Red\. masses --\s+(.*)")
RE_FRC_CONST = re.compile(r"Frc consts\s*--\s+(.*)")
RE_IR_INTEN = re.compile(r"IR Inten\s*--\s+(.*)")
RE_ATOM_HEADER = re.compile(r"^\s*Atom\s+AN\s+X\s+Y\s+Z")
# A displacement row: idx  AN  then 3*k floats
RE_DISP_ROW = re.compile(r"^\s*\d+\s+\d+\s+[-\d.]")


class GaussianLogError(Exception):
    """Raised when a Gaussian log cannot be parsed as a valid freq calculation.

    The CLI/driver routes these to the bad-logs manifest instead of crashing.
    """


# --------------------------------------------------------------------------- #
# XYZ reader
# --------------------------------------------------------------------------- #
def extract_frequency(filename: str | os.PathLike) -> str | None:
    """Extract a trailing frequency from a filename, e.g. m_165_3150.62.xyz -> 3150.62."""
    name = os.path.splitext(os.path.basename(os.fspath(filename)))[0]
    nums = re.findall(r"\d+\.\d+|\d+", name)
    return nums[-1] if nums else None


def read_xyz(path: str | os.PathLike) -> tuple[Geometry, list[Mode]]:
    """Read an ``Element X Y Z [dx dy dz]`` XYZ file into one Geometry + one Mode."""
    path = Path(path)
    symbols: list[str] = []
    positions: list[list[float]] = []
    disps: list[list[float]] = []

    for line in path.read_text(errors="ignore").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 4:
            continue
        symbols.append(parts[0])
        positions.append([float(x) for x in parts[1:4]])
        if len(parts) >= 7:
            disps.append([float(x) for x in parts[4:7]])
        else:
            disps.append([0.0, 0.0, 0.0])

    if not symbols:
        raise ValueError(f"No atoms found in XYZ file: {path}")

    geometry = Geometry(symbols=symbols,
                        positions=np.array(positions, dtype=float),
                        source_frame="xyz")

    freq_str = extract_frequency(path.name)
    frequency = float(freq_str) if freq_str is not None else 0.0
    mode = Mode(
        index=1,
        frequency=frequency,
        ir_intensity=None,
        reduced_mass=None,
        force_const=None,
        displacements=np.array(disps, dtype=float),
        is_imaginary=frequency < 0,
    )
    return geometry, [mode]


# --------------------------------------------------------------------------- #
# Gaussian log reader
# --------------------------------------------------------------------------- #
def _parse_orientation_table(lines: list[str], header_idx: int, path: Path) -> tuple[list[int], np.ndarray]:
    """Parse an orientation table starting at the header line index.

    Layout:
        Standard orientation:
        -----
        Center  Atomic  Atomic   Coordinates (Angstroms)
        Number  Number   Type      X       Y       Z
        -----
        1   7   0   x y z
        ...
        -----
    Returns (atomic_numbers, positions (n,3)).
    """
    # Skip to the first data row: past two dashed lines after the header.
    i = header_idx + 1
    dash_count = 0
    while i < len(lines):
        if set(lines[i].strip()) <= {"-"} and lines[i].strip():
            dash_count += 1
            i += 1
            if dash_count == 2:  # after column-name block
                break
            continue
        i += 1

    atomic_numbers: list[int] = []
    positions: list[list[float]] = []
    while i < len(lines):
        stripped = lines[i].strip()
        if set(stripped) <= {"-"} and stripped:  # closing dashed line
            break
        parts = stripped.split()
        if len(parts) >= 6 and parts[0].isdigit():
            try:
                atomic_numbers.append(int(parts[1]))
                positions.append([float(parts[3]), float(parts[4]), float(parts[5])])
            except ValueError as exc:
                raise GaussianLogError(f"Malformed orientation numeric value in {path.name} at line {i + 1}") from exc
        i += 1
    return atomic_numbers, np.array(positions, dtype=float)


def _parse_floats(text: str, path: Path, context: str) -> list[float]:
    try:
        return [float(x) for x in text.split()]
    except ValueError as exc:
        raise GaussianLogError(f"Malformed numeric value in {path.name} ({context})") from exc


def _is_linear(positions: np.ndarray, tolerance: float = 1e-5) -> bool:
    """Return whether atom positions are collinear, using a scale-aware SVD test."""
    if len(positions) <= 2:
        return True
    singular_values = np.linalg.svd(positions - positions.mean(axis=0), compute_uv=False)
    return bool(singular_values[1] <= tolerance * max(float(singular_values[0]), 1.0))


def _block_context(mode_numbers: list[int], freqs: list[float]) -> str:
    return f"modes {mode_numbers}, frequencies {freqs}"


def _parse_frequency_section(lines: list[str], path: Path, start: int, end: int,
                             harmonic_idx: int) -> tuple[Geometry, list[Mode]]:
    """Parse one bounded Gaussian 16 frequency-job section strictly."""
    std_indices = [i for i in range(start, harmonic_idx) if RE_STD_ORIENT.search(lines[i])]
    input_indices = [i for i in range(start, harmonic_idx) if RE_INPUT_ORIENT.search(lines[i])]
    if std_indices:
        orient_idx, source_frame = std_indices[-1], "standard"
    elif input_indices:
        orient_idx, source_frame = input_indices[-1], "input"
    else:
        raise GaussianLogError(f"No orientation table before frequency section in {path.name}")
    atomic_numbers, positions = _parse_orientation_table(lines, orient_idx, path)
    if positions.size == 0:
        raise GaussianLogError(f"Empty orientation table in {path.name}")
    geometry = Geometry([symbol_from_z(z) for z in atomic_numbers], positions, source_frame)
    n_atoms = geometry.n_atoms
    if n_atoms == 1:
        raise GaussianLogError(f"{path.name}: atomic frequency job has no vibrational modes")

    modes: list[Mode] = []
    i = harmonic_idx
    while i < end:
        match = RE_FREQ_LINE.search(lines[i])
        if not match:
            i += 1
            continue
        mode_numbers: list[int] = []
        for back in range(i - 1, max(i - 4, -1), -1):
            tokens = lines[back].split()
            if tokens and all(token.isdigit() for token in tokens):
                mode_numbers = [int(token) for token in tokens]
                break
        freqs = _parse_floats(match.group(1), path, "Frequencies -- row")
        ncols = len(freqs)
        if not ncols:
            raise GaussianLogError(f"Empty Frequencies -- row in {path.name}")
        if not mode_numbers or len(mode_numbers) != ncols:
            first = modes[-1].index + 1 if modes else 1
            mode_numbers = list(range(first, first + ncols))
        context = _block_context(mode_numbers, freqs)
        red_masses: list[float | None] = [None] * ncols
        frc_consts: list[float | None] = [None] * ncols
        ir_inten: list[float | None] = [None] * ncols
        j = i + 1
        while j < end and not RE_ATOM_HEADER.search(lines[j]):
            if RE_FREQ_LINE.search(lines[j]):
                break
            if (metadata := RE_RED_MASS.search(lines[j])):
                red_masses = _pad(_parse_floats(metadata.group(1), path, "Red. masses row"), ncols)
            elif (metadata := RE_FRC_CONST.search(lines[j])):
                frc_consts = _pad(_parse_floats(metadata.group(1), path, "Frc consts row"), ncols)
            elif (metadata := RE_IR_INTEN.search(lines[j])):
                ir_inten = _pad(_parse_floats(metadata.group(1), path, "IR Inten row"), ncols)
            j += 1
        if j >= end or not RE_ATOM_HEADER.search(lines[j]):
            raise GaussianLogError(f"{path.name}: incomplete frequency block ({context}): missing Atom AN header")
        disp_block = np.empty((ncols, n_atoms, 3), dtype=float)
        block_atomic_numbers: list[int] = []
        for expected_index in range(1, n_atoms + 1):
            row_idx = j + expected_index
            if row_idx >= end or not RE_DISP_ROW.search(lines[row_idx]):
                raise GaussianLogError(f"{path.name}: incomplete frequency block ({context}): expected {n_atoms} displacement rows")
            tokens = lines[row_idx].split()
            expected_width = 2 + 3 * ncols
            if len(tokens) != expected_width:
                raise GaussianLogError(f"{path.name}: incomplete frequency block ({context}): displacement row {expected_index} has {len(tokens) - 2} values, expected {3 * ncols}")
            try:
                atom_index, atomic_number = int(tokens[0]), int(tokens[1])
                values = [float(value) for value in tokens[2:]]
            except ValueError as exc:
                raise GaussianLogError(f"Malformed numeric value in {path.name} ({context}, displacement row {expected_index})") from exc
            if atom_index != expected_index:
                raise GaussianLogError(f"{path.name}: frequency block ({context}) has atom index {atom_index}; expected {expected_index}")
            block_atomic_numbers.append(atomic_number)
            for col in range(ncols):
                disp_block[col, expected_index - 1] = values[3 * col:3 * col + 3]
        if block_atomic_numbers != atomic_numbers:
            raise GaussianLogError(f"{path.name}: frequency block ({context}) atomic numbers do not match orientation table")
        for col, frequency in enumerate(freqs):
            modes.append(Mode(mode_numbers[col], frequency, ir_inten[col], red_masses[col],
                              frc_consts[col], disp_block[col], frequency < 0))
        i = j + n_atoms + 1
    if not modes:
        raise GaussianLogError(f"No complete frequency blocks in {path.name}")
    linear = _is_linear(positions)
    expected = 3 * n_atoms - (5 if linear else 6)
    if len(modes) != expected:
        kind = "linear (3N-5)" if linear else "nonlinear (3N-6)"
        raise GaussianLogError(f"{path.name}: parsed {len(modes)} modes; expected {expected} for a {kind} molecule")
    return geometry, modes


def read_gaussian_log(path: str | os.PathLike) -> tuple[Geometry, list[Mode]]:
    """Parse the last complete low-precision Gaussian 16 frequency job.

    Link1/concatenated files are evaluated job-by-job; an incomplete final job
    falls back to the preceding complete frequency job. HPModes is unsupported.
    """
    path = Path(path)
    text = path.read_text(errors="ignore")
    lines = text.splitlines()

    harmonic_indices = [i for i, line in enumerate(lines) if RE_HARMONIC.search(line)]
    if not harmonic_indices:
        raise GaussianLogError(f"No 'Harmonic frequencies' section in {path.name}")
    link_indices = [i for i, line in enumerate(lines) if RE_LINK1.search(line)]
    last_error: GaussianLogError | None = None
    for harmonic_idx in reversed(harmonic_indices):
        section_start = max((i + 1 for i in link_indices if i < harmonic_idx), default=0)
        later_headers = [i for i in harmonic_indices if i > harmonic_idx]
        later_links = [i for i in link_indices if i > harmonic_idx]
        section_end = min(later_headers + later_links + [len(lines)])
        if not any(RE_NORMAL_TERM.search(line) for line in lines[harmonic_idx:section_end]):
            last_error = GaussianLogError(f"No 'Normal termination' for frequency job in {path.name}")
            continue
        try:
            return _parse_frequency_section(lines, path, section_start, section_end, harmonic_idx)
        except GaussianLogError as exc:
            last_error = exc
    raise last_error or GaussianLogError(f"No complete frequency section in {path.name}")


def _pad(values: list[float], n: int) -> list[float | None]:
    """Pad/truncate a value list to length n (None-fill)."""
    out: list[float | None] = list(values[:n])
    out += [None] * (n - len(out))
    return out


# --------------------------------------------------------------------------- #
# Dispatcher
# --------------------------------------------------------------------------- #
def read_any(path: str | os.PathLike) -> tuple[Geometry, list[Mode]]:
    """Dispatch to the right reader based on file suffix."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".xyz":
        return read_xyz(path)
    if suffix in (".log", ".out"):
        return read_gaussian_log(path)
    raise ValueError(f"Unsupported file type: {suffix!r} ({path})")
