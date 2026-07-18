"""Tests for the Gaussian-log and XYZ readers, mode selection, and rendering."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from molvib.model import Mode
from molvib.readers import GaussianLogError, read_any, read_gaussian_log, read_xyz
from molvib.select import select_modes
from molvib.cli import render_directory

DATA = Path(__file__).parent / "data"
ADENINE_LOG = DATA / "adenine_pcm_water.log"


def _frequency_job(coords=((0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), frequency=100.0,
                   rows=None, atom_numbers=(1, 1), header=True, terminate=True,
                   malformed_frequency=False, extra_metadata=False):
    """Small, valid two-atom linear Gaussian-like frequency job."""
    orientation = [" Standard orientation:", " ---------------------------------------------------------------------",
                   " Center     Atomic      Atomic             Coordinates (Angstroms)",
                   " Number     Number       Type             X           Y           Z",
                   " ---------------------------------------------------------------------"]
    orientation += [f" {i:5d} {z:10d}    0 {xyz[0]:12.6f} {xyz[1]:12.6f} {xyz[2]:12.6f}"
                    for i, (z, xyz) in enumerate(zip(atom_numbers, coords), 1)]
    orientation += [" ---------------------------------------------------------------------",
                    " Harmonic frequencies (cm**-1), IR intensities (KM/Mole)",
                    "                      1", " Frequencies --    BAD" if malformed_frequency
                    else f" Frequencies --    {frequency:.4f}",
                    " Red. masses --      1.0000", " Frc consts  --      0.1000",
                    " IR Inten    --      5.0000"]
    if extra_metadata:
        orientation.append(" Raman Activ --      2.0000")
    if header:
        orientation.append("  Atom  AN      X      Y      Z")
    if rows is None:
        rows = [(1, atom_numbers[0], "0.1 0.0 0.0"), (2, atom_numbers[1], "-0.1 0.0 0.0")]
    orientation += [f" {idx:5d} {z:3d} {values}" for idx, z, values in rows]
    if terminate:
        orientation.append(" Normal termination of Gaussian 16")
    return "\n".join(orientation) + "\n"


# --------------------------------------------------------------------------- #
# Gaussian log parsing
# --------------------------------------------------------------------------- #
def test_parse_adenine_basic():
    geometry, modes = read_gaussian_log(ADENINE_LOG)
    assert geometry.n_atoms == 15
    assert len(modes) == 39
    assert geometry.source_frame == "standard"


def test_atom_symbol_order():
    geometry, _ = read_gaussian_log(ADENINE_LOG)
    # adenine C5H5N5: orientation table lists 5 N, then 5 C, then 5 H.
    assert geometry.symbols == (
        ["N"] * 5 + ["C"] * 5 + ["H"] * 5
    )


def test_first_and_last_frequency():
    _, modes = read_gaussian_log(ADENINE_LOG)
    assert modes[0].frequency == pytest.approx(157.5304, abs=1e-3)
    assert modes[-1].frequency == pytest.approx(3708.8798, abs=1e-3)


def test_displacement_shape_per_mode():
    geometry, modes = read_gaussian_log(ADENINE_LOG)
    for m in modes:
        assert m.displacements.shape == (15, 3)


def test_mode_metadata_present():
    _, modes = read_gaussian_log(ADENINE_LOG)
    m0 = modes[0]
    assert m0.index == 1
    assert m0.ir_intensity == pytest.approx(3.9622, abs=1e-3)
    assert m0.reduced_mass == pytest.approx(4.8989, abs=1e-3)
    assert m0.force_const == pytest.approx(0.0716, abs=1e-3)
    assert not m0.is_imaginary


def test_frame_matching_uses_last_standard_before_freq():
    """Regression: positions must come from the standard-orientation table
    immediately preceding the freq section, not input orientation."""
    geometry, _ = read_gaussian_log(ADENINE_LOG)
    # First atom (N, center 1) from the last Standard orientation table.
    np.testing.assert_allclose(
        geometry.positions[0], [-1.278202, -1.324006, -0.004147], atol=1e-5)
    # This differs from the corresponding Input orientation coordinates.


def test_read_any_dispatch():
    geometry, modes = read_any(ADENINE_LOG)
    assert len(modes) == 39
    with pytest.raises(ValueError):
        read_any(Path("foo.pdb"))


def test_bad_log_raises(tmp_path):
    bad = tmp_path / "bad.log"
    bad.write_text("garbage, no normal termination here")
    with pytest.raises(GaussianLogError):
        read_gaussian_log(bad)


@pytest.mark.parametrize("kwargs", [
    {"header": False},
    {"rows": [(1, 1, "0.1 0.0 0.0")]},
    {"rows": [(1, 1, "0.1 0.0"), (2, 1, "-0.1 0.0 0.0")]},
    {"rows": [(1, 8, "0.1 0.0 0.0"), (2, 1, "-0.1 0.0 0.0")]},
])
def test_incomplete_frequency_blocks_raise_and_route_to_bad_logs(tmp_path, kwargs):
    bad = tmp_path / "bad.log"
    bad.write_text(_frequency_job(**kwargs))
    with pytest.raises(GaussianLogError):
        read_gaussian_log(bad)
    manifest, bad_logs = render_directory(tmp_path, outdir=tmp_path / "out", render=False)
    assert manifest.empty
    assert bad_logs["source_log"].tolist() == ["bad.log"]


def test_raman_metadata_and_out_suffix(tmp_path):
    out = tmp_path / "raman.out"
    out.write_text(_frequency_job(extra_metadata=True))
    geometry, modes = read_any(out)
    assert geometry.n_atoms == 2
    assert len(modes) == 1


def test_malformed_numeric_value_has_clear_error(tmp_path):
    bad = tmp_path / "numeric.log"
    bad.write_text(_frequency_job(malformed_frequency=True))
    with pytest.raises(GaussianLogError, match="Malformed numeric value.*numeric.log"):
        read_gaussian_log(bad)


def test_linear_molecule_is_accepted(tmp_path):
    log = tmp_path / "linear.log"
    log.write_text(_frequency_job())
    geometry, modes = read_gaussian_log(log)
    assert geometry.n_atoms == 2
    assert len(modes) == 1


def test_last_complete_link1_frequency_job_wins(tmp_path):
    log = tmp_path / "two_jobs.log"
    log.write_text(
        _frequency_job(coords=((0, 0, 0), (0, 0, 1)), frequency=111.0)
        + " --Link1--\n"
        + _frequency_job(coords=((5, 0, 0), (5, 0, 2)), frequency=222.0))
    geometry, modes = read_gaussian_log(log)
    np.testing.assert_allclose(geometry.positions[0], [5, 0, 0])
    assert [mode.frequency for mode in modes] == [pytest.approx(222.0)]


def test_incomplete_final_link1_falls_back_to_complete_job(tmp_path):
    log = tmp_path / "fallback.log"
    log.write_text(_frequency_job(frequency=111.0) + " --Link1--\n" +
                   _frequency_job(frequency=222.0, header=False, terminate=False))
    _, modes = read_gaussian_log(log)
    assert [mode.frequency for mode in modes] == [pytest.approx(111.0)]


# --------------------------------------------------------------------------- #
# XYZ reading
# --------------------------------------------------------------------------- #
def test_read_xyz(tmp_path):
    f = tmp_path / "m_1_1650.30.xyz"
    f.write_text(
        "O 0.0 0.0 0.0 0.1 0.0 0.0\n"
        "H 0.9 0.0 0.0 -0.1 0.0 0.0\n"
        "H -0.3 0.8 0.0 0.0 -0.1 0.0\n"
    )
    geometry, modes = read_xyz(f)
    assert geometry.symbols == ["O", "H", "H"]
    assert geometry.source_frame == "xyz"
    assert len(modes) == 1
    assert modes[0].frequency == pytest.approx(1650.30)
    assert modes[0].ir_intensity is None
    assert modes[0].displacements.shape == (3, 3)


# --------------------------------------------------------------------------- #
# Mode selection
# --------------------------------------------------------------------------- #
def _mk(index, freq, ir):
    return Mode(index=index, frequency=freq, ir_intensity=ir, reduced_mass=None,
               force_const=None, displacements=np.zeros((1, 3)), is_imaginary=freq < 0)


def test_select_ir_threshold():
    modes = [_mk(1, 100, 5.0), _mk(2, 200, 50.0), _mk(3, 300, 15.0)]
    out = select_modes(modes, ir_threshold=10.0, include_imaginary=True)
    assert [m.index for m in out] == [2, 3]


def test_select_imaginary_always_kept():
    modes = [_mk(1, -50, 0.0), _mk(2, 200, 1.0)]
    out = select_modes(modes, ir_threshold=10.0, include_imaginary=True)
    assert [m.index for m in out] == [1]  # imaginary kept, low-IR real dropped
    out2 = select_modes(modes, ir_threshold=10.0, include_imaginary=False)
    assert out2 == []


def test_select_explicit_indices_overrides():
    modes = [_mk(1, 100, 5.0), _mk(2, 200, 50.0), _mk(3, 300, 15.0)]
    out = select_modes(modes, ir_threshold=10.0, indices=[1, 3])
    assert [m.index for m in out] == [1, 3]


def test_select_none_ir_kept():
    modes = [_mk(1, 100, None)]
    out = select_modes(modes, ir_threshold=10.0)
    assert [m.index for m in out] == [1]


# --------------------------------------------------------------------------- #
# Smoke render (requires pyvista; skipped if unavailable)
# --------------------------------------------------------------------------- #
def test_smoke_render(tmp_path):
    pv = pytest.importorskip("pyvista")
    from molvib.model import RenderJob
    from molvib.render import RenderSettings, render_mode

    geometry, modes = read_gaussian_log(ADENINE_LOG)
    # A high-frequency N-H/C-H stretch to eyeball arrow direction.
    mode = max(modes, key=lambda m: m.frequency)
    job = RenderJob(geometry=geometry, mode=mode,
                   source_file=ADENINE_LOG, label=mode.label)
    settings = RenderSettings(resolution=(300, 300))
    try:
        out = render_mode(job, settings=settings, outdir=tmp_path)
    except Exception as exc:  # headless env without a GL backend
        pytest.skip(f"pyvista render unavailable: {exc}")
    assert out.exists()
    assert out.stat().st_size > 0
