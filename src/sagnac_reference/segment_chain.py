from __future__ import annotations

from dataclasses import dataclass
import math

from .analytic import SagnacConfig


@dataclass(frozen=True)
class ChainResult:
    return_time: float
    steps: int
    dt: float
    winding_target: float


def _signed_relative_increment(cfg: SagnacConfig, direction: int, dt: float) -> float:
    # Unwrapped light-minus-detector displacement in the inertial axle frame.
    # +1: light and detector move in same orientation -> (c-v)dt.
    # -1: light moves opposite detector -> -(c+v)dt.
    if direction == +1:
        return (cfg.c - cfg.v) * dt
    if direction == -1:
        return -(cfg.c + cfg.v) * dt
    raise ValueError("direction must be +1 or -1")


def return_time_segment_chain(cfg: SagnacConfig, direction: int, segments: int = 10000) -> ChainResult:
    """Discrete moving-closure propagation with event interpolation.

    The algorithm advances an unwrapped relative coordinate and detects the
    first non-zero winding by a sign-safe crossing. It intentionally does not
    call the analytic return-time formula.
    """
    if segments < 8:
        raise ValueError("segments must be >= 8")
    dx = cfg.L / segments
    dt = dx / cfg.c
    target = cfg.L if direction == +1 else -cfg.L
    rel = 0.0
    t = 0.0
    steps = 0

    while steps < 10 * segments:
        prev_rel = rel
        prev_t = t
        rel += _signed_relative_increment(cfg, direction, dt)
        t += dt
        steps += 1
        crossed = (direction == +1 and rel >= target) or (direction == -1 and rel <= target)
        if crossed:
            denom = rel - prev_rel
            if denom == 0:
                raise RuntimeError("zero relative increment")
            alpha = (target - prev_rel) / denom
            return ChainResult(prev_t + alpha * dt, steps, dt, target)
    raise RuntimeError("return event not found")
