from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass(frozen=True)
class MixingResult:
    """
    Mixing result container.

    Attributes
    ----------
    a_mix : float
        The mixture parameter a.
    b_mix : float
        The mixture parameter b.
    a_ij : np.ndarray
        The matrix of interaction parameters a_ij.
    A_mix : float | None
        The mixture parameter A (optional).
    B_mix : float | None
        The mixture parameter B (optional).
    metadata : dict[str, Any]
        Additional metadata related to the mixture calculation.
    """
    a_mix: float
    b_mix: float
    a_ij: np.ndarray
    A_mix: float | None = None
    B_mix: float | None = None
    metadata: dict[str, Any] = field(default_factory=dict)
