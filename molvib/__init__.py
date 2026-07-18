"""molvib — Gaussian-log + XYZ molecular vibration rendering pipeline.

Public API:
    read_any, read_xyz, read_gaussian_log     (readers)
    select_modes                              (select)
    render_mode                               (render)
    build_manifest                            (manifest)
    render_file, render_directory             (cli)
"""

from .model import Geometry, Mode, RenderJob
from .readers import read_any, read_xyz, read_gaussian_log, GaussianLogError
from .select import select_modes
from .manifest import build_manifest

__all__ = [
    "Geometry",
    "Mode",
    "RenderJob",
    "read_any",
    "read_xyz",
    "read_gaussian_log",
    "GaussianLogError",
    "select_modes",
    "build_manifest",
]

__version__ = "0.1.0"
