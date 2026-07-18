"""Mode selection: IR-intensity threshold + imaginary-mode handling."""

from __future__ import annotations

from collections.abc import Iterable

from .model import Mode


def select_modes(
    modes: list[Mode],
    ir_threshold: float | None = None,
    include_imaginary: bool = True,
    indices: Iterable[int] | None = None,
) -> list[Mode]:
    """Select which modes to render.

    - ``indices``: if given, an explicit set of 1-based Gaussian mode numbers.
      Overrides IR filtering entirely.
    - ``ir_threshold``: keep modes with ``ir_intensity >= ir_threshold``
      (KM/Mole). Modes with no IR data (``ir_intensity is None``, e.g. from XYZ)
      bypass the filter and are always kept.
    - ``include_imaginary``: when True, imaginary modes (frequency < 0) are
      always kept regardless of IR intensity — they are the diagnostically
      important ones.

    Order is preserved from the input list.
    """
    if indices is not None:
        wanted = set(indices)
        return [m for m in modes if m.index in wanted]

    selected: list[Mode] = []
    for m in modes:
        if include_imaginary and m.is_imaginary:
            selected.append(m)
            continue
        if ir_threshold is None:
            selected.append(m)
            continue
        if m.ir_intensity is None:  # no intensity data -> keep
            selected.append(m)
            continue
        if m.ir_intensity >= ir_threshold:
            selected.append(m)
    return selected
