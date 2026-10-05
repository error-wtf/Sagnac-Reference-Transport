from __future__ import annotations

import math

from .analytic import SagnacConfig, delta_tau_detector


def detector_phase(cfg: SagnacConfig, omega_detector: float) -> float:
    if not math.isfinite(omega_detector):
        raise ValueError("omega_detector must be finite")
    return omega_detector * delta_tau_detector(cfg)
