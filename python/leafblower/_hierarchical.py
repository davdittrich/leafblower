"""HierarchicalConfig — typed config for two-stage coarse-then-fine calibration.

Mirrors the R-side ``hierarchical`` list accepted by ``harvest(hierarchical=...)``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True)
class HierarchicalConfig:
    """Configuration for two-stage (coarse → fine) hierarchical calibration.

    Parameters
    ----------
    coarse_margins:
        Names of the margins that define the coarse-stage partition.  Every
        name must appear as a key in the ``targets`` dict passed to
        ``harvest()``.  At least one name is required.
    min_cell_n:
        Minimum number of observations a coarse cell must contain to qualify
        for Stage-2 refinement.  Cells with fewer observations inherit the
        coarse-stage weights unchanged.  Must be >= 1.
    mode:
        ``"refine"`` (default) — iterates Stage-1 / Stage-2 until
        ``outer_tol`` is met or ``outer_iterations`` is exhausted.
        ``"exact"`` — single-pass; requires orthogonal coarse/fine split
        (BADARG if not orthogonal).
    outer_tol:
        Outer-loop convergence tolerance.  Must be > 0.
    outer_iterations:
        Maximum outer-loop iterations.  Must be in [1, 10000].
    """

    coarse_margins: list[str]
    min_cell_n: int = 30
    mode: Literal["refine", "exact"] = "refine"
    outer_tol: float = 1e-4
    outer_iterations: int = 10

    def __post_init__(self) -> None:
        # coarse_margins: non-empty list of strings
        if not isinstance(self.coarse_margins, list) or len(self.coarse_margins) == 0:
            raise ValueError(
                "HierarchicalConfig: coarse_margins must be a non-empty list of margin names"
            )
        if not all(isinstance(m, str) for m in self.coarse_margins):
            raise ValueError(
                "HierarchicalConfig: all coarse_margins entries must be strings"
            )

        # mode enum
        if self.mode not in ("refine", "exact"):
            raise ValueError(
                f"HierarchicalConfig: mode must be 'refine' or 'exact', got {self.mode!r}"
            )

        # min_cell_n >= 1  (mirrors R-side BADARG: min_cell_n < 1)
        if not isinstance(self.min_cell_n, int) or self.min_cell_n < 1:
            raise ValueError(
                f"HierarchicalConfig: min_cell_n must be an integer >= 1, got {self.min_cell_n!r}"
            )

        # outer_iterations in [1, 10000]
        if not isinstance(self.outer_iterations, int) or not (1 <= self.outer_iterations <= 10000):
            raise ValueError(
                f"HierarchicalConfig: outer_iterations must be in [1, 10000], "
                f"got {self.outer_iterations!r}"
            )

        # outer_tol > 0
        if not (self.outer_tol > 0):
            raise ValueError(
                f"HierarchicalConfig: outer_tol must be > 0, got {self.outer_tol!r}"
            )
