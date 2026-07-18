"""argparse driver + importable render_file / render_directory helpers."""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import pandas as pd

from .manifest import build_bad_logs, build_manifest, manifest_row
from .model import RenderJob
from .readers import GaussianLogError, read_any
from .render import RenderSettings, render_mode
from .select import select_modes

SUFFIXES = (".xyz", ".log", ".out")


def _iter_inputs(input_path: Path) -> list[Path]:
    if input_path.is_dir():
        files: list[Path] = []
        for suf in SUFFIXES:
            files.extend(sorted(input_path.glob(f"*{suf}")))
        return files
    return [input_path]


def _jobs_for_file(path: Path, ir_threshold, include_imaginary, indices):
    """Read + select; return (geometry, [RenderJob])."""
    geometry, modes = read_any(path)
    chosen = select_modes(modes, ir_threshold=ir_threshold,
                          include_imaginary=include_imaginary, indices=indices)
    jobs = [RenderJob(geometry=geometry, mode=m, source_file=path, label=m.label)
            for m in chosen]
    return geometry, jobs


def render_file(
    path: str | Path,
    outdir: str | Path | None = None,
    ir_threshold: float | None = 10.0,
    include_imaginary: bool = True,
    indices=None,
    settings: RenderSettings | None = None,
    render: bool = True,
) -> pd.DataFrame:
    """Render all selected modes of one file. Returns the manifest DataFrame."""
    path = Path(path)
    out_dir = Path(outdir) if outdir is not None else path.parent
    settings = settings or RenderSettings()

    rows: list[dict] = []
    geometry, jobs = _jobs_for_file(path, ir_threshold, include_imaginary, indices)
    for job in jobs:
        png = None
        rendered = False
        if render:
            png = render_mode(job, settings=settings, outdir=out_dir)
            rendered = True
        rows.append(manifest_row(path, geometry, job.mode, png, rendered))
    return build_manifest(rows)


def render_directory(
    directory: str | Path,
    outdir: str | Path | None = None,
    ir_threshold: float | None = 10.0,
    include_imaginary: bool = True,
    indices=None,
    settings: RenderSettings | None = None,
    render: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Render every input file in a directory.

    Returns ``(manifest_df, bad_logs_df)``. Files that fail validation are
    routed to the bad-logs DataFrame instead of aborting the run.
    """
    directory = Path(directory)
    out_dir = Path(outdir) if outdir is not None else directory
    settings = settings or RenderSettings()

    rows: list[dict] = []
    bad_rows: list[dict] = []
    for path in _iter_inputs(directory):
        try:
            df = render_file(path, outdir=out_dir, ir_threshold=ir_threshold,
                             include_imaginary=include_imaginary, indices=indices,
                             settings=settings, render=render)
            rows.extend(df.to_dict("records"))
        except (GaussianLogError, ValueError) as exc:
            bad_rows.append({"source_log": path.name, "reason": str(exc)})
    return build_manifest(rows), build_bad_logs(bad_rows)


# --------------------------------------------------------------------------- #
# argparse
# --------------------------------------------------------------------------- #
def _parse_resolution(s: str) -> tuple[int, int]:
    w, _, h = s.lower().partition("x")
    return int(w), int(h)


def _parse_indices(s: str | None):
    if not s:
        return None
    return [int(x) for x in s.split(",") if x.strip()]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="molvib",
        description="Render molecular vibration PNGs from XYZ or Gaussian 16 log files.")
    p.add_argument("--input", required=True,
                   help="file or directory (globs *.xyz, *.log, *.out)")
    p.add_argument("--ir-threshold", type=float, default=10.0,
                   help="IR intensity threshold in KM/Mole (default: 10.0)")
    imag = p.add_mutually_exclusive_group()
    imag.add_argument("--include-imaginary", dest="include_imaginary",
                      action="store_true", default=True,
                      help="always keep imaginary modes (default)")
    imag.add_argument("--no-include-imaginary", dest="include_imaginary",
                      action="store_false")
    p.add_argument("--modes", default=None,
                   help="explicit mode indices (comma list); overrides IR filter")
    p.add_argument("--outdir", default=None,
                   help="output directory (default: alongside input)")
    p.add_argument("--manifest", default=None,
                   help="manifest CSV path (default: outdir/vibration_manifest.csv)")
    p.add_argument("--resolution", type=_parse_resolution, default=(2000, 2000),
                   help="WxH, default 2000x2000")
    p.add_argument("--background", default="transparent",
                   help="transparent|<color>")
    p.add_argument("--zoom", type=float, default=0.9,
                   help="camera zoom; <1 adds margin so atoms/arrows aren't cropped (default 0.9)")
    p.add_argument("--no-arrows", dest="show_arrows", action="store_false",
                   default=True, help="disable displacement arrows")
    p.add_argument("--no-render", dest="render", action="store_false", default=True,
                   help="parse + select + write manifest without rendering PNGs")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    input_path = Path(args.input)
    out_dir = Path(args.outdir) if args.outdir else (
        input_path if input_path.is_dir() else input_path.parent)

    settings = RenderSettings(
        resolution=args.resolution,
        background=args.background,
        show_arrows=args.show_arrows,
        camera_zoom=args.zoom,
    )
    indices = _parse_indices(args.modes)

    if input_path.is_dir():
        manifest_df, bad_df = render_directory(
            input_path, outdir=out_dir, ir_threshold=args.ir_threshold,
            include_imaginary=args.include_imaginary, indices=indices,
            settings=settings, render=args.render)
    else:
        try:
            manifest_df = render_file(
                input_path, outdir=out_dir, ir_threshold=args.ir_threshold,
                include_imaginary=args.include_imaginary, indices=indices,
                settings=settings, render=args.render)
            bad_df = build_bad_logs([])
        except (GaussianLogError, ValueError) as exc:
            manifest_df = build_manifest([])
            bad_df = build_bad_logs([{"source_log": input_path.name, "reason": str(exc)}])

    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(args.manifest) if args.manifest else out_dir / "vibration_manifest.csv"
    manifest_df.to_csv(manifest_path, index=False)
    print(f"[OK] Wrote manifest: {manifest_path} ({len(manifest_df)} rows)")

    if len(bad_df):
        bad_path = out_dir / "vibration_bad_logs.csv"
        bad_df.to_csv(bad_path, index=False)
        warnings.warn(f"{len(bad_df)} file(s) failed validation; see {bad_path}")
        print(f"[WARN] Wrote bad-logs manifest: {bad_path} ({len(bad_df)} rows)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
