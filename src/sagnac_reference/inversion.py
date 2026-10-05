from __future__ import annotations

import math


def velocity_from_times(t_plus: float, t_minus: float, c: float) -> float:
    if c <= 0 or t_plus <= 0 or t_minus <= 0:
        raise ValueError("c, t_plus and t_minus must be > 0")
    return c * (t_plus - t_minus) / (t_plus + t_minus)


def velocity_from_delta_t(delta_t: float, L: float, c: float) -> float:
    """Stable physical-root inversion of Δt = 2Lv/(c²-v²)."""
    if L <= 0 or c <= 0:
        raise ValueError("L and c must be > 0")
    if delta_t == 0:
        return 0.0
    root = math.sqrt(L * L + (c * delta_t) ** 2)
    return (c * c * delta_t) / (L + root)
