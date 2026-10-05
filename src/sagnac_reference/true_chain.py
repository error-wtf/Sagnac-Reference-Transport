"""Sagnac reference transport - TRUE discrete ring simulation (Phase 4 rewrite).

The baseline `segment_chain.return_time_segment_chain` advances the RELATIVE
coordinate directly with (c - v)*dt — that is the analytic relative-speed
formula in discrete clothing, not an independent simulation.

This module instead evolves LIGHT POSITION and DETECTOR POSITION SEPARATELY
on the periodic ring [0, L):

    x_light(t + h) = (x_light(t) + s*c*h) mod L        (light, speed ±c)
    x_detector(t + h) = (x_detector(t) + v*h) mod L    (detector, speed v)

Return event: first time the wrapped light position crosses the detector
position (in the direction of light travel), located by linear interpolation
between the two bracketing time levels.  The closed-form Sagnac time is
never used; the code knows only L, c, v, the two kinematic update rules and
the periodic topology.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .analytic import SagnacConfig


@dataclass(frozen=True)
class TrueChainResult:
    return_time: float
    steps: int
    h: float
    wrap_count: int


def _forward_diff(a: float, b: float, L: float) -> float:
    """Forward difference (b - a) wrapped into [0, L)."""
    return (b - a) % L


def return_time_true_chain(
    cfg: SagnacConfig,
    direction: int,
    *,
    steps: int = 20000,
) -> TrueChainResult:
    """Simulate light and detector positions separately; find first re-encounter.

    The light wraps the ring; the detector also drifts (v).  A 'return' means:
    after the light has completed at least one full wrap relative to its start
    (wrap_count >= 1), it arrives at the detector's CURRENT position.
    Crossing detection uses the forward difference of the wrapped light-
    minus-detector separation with linear interpolation in time.
    """
    if direction not in (+1, -1):
        raise ValueError("direction must be +1 or -1")
    if steps < 16:
        raise ValueError("steps must be >= 16")

    # `steps` is the ANGULAR resolution (light-steps per full ring); the
    # simulation length is 4 full light-wraps, which covers any |v| < c case
    # (the detector can never outrun the light by more than one extra wrap).
    h = cfg.L / cfg.c / steps
    max_steps = 10 * steps + 16
    x_light = 0.0
    x_det = 0.0
    t = 0.0
    s = float(direction)

    # separation state: forward difference light -> detector in travel direction
    def sep(xl: float, xd: float) -> float:
        # distance from light to detector measured ALONG light's direction
        d = _forward_diff(xl, xd, cfg.L) if direction == +1 else _forward_diff(xd, xl, cfg.L)
        return d

    prev_sep = sep(x_light, x_det)     # ~L right after the first step
    prev_t = 0.0

    for n in range(1, max_steps + 1):
        x_light = (x_light + s * cfg.c * h) % cfg.L
        x_det = (x_det + cfg.v * h) % cfg.L
        t = n * h
        cur_sep = sep(x_light, x_det)
        # sep shrinks monotonically while the light chases the detector.
        # Return event: sep wraps from small positive back up to ~L, i.e.
        # cur_sep > prev_sep after the chase phase has begun (prev_sep < 0.5L
        # guard keeps the initial jump 0 -> ~L from firing).
        if n >= 2 and cur_sep > prev_sep:
            # linear interpolation of the zero crossing between prev_t and t:
            # prev_sep above zero, wrapped cur_sep = cur_sep - L below zero.
            frac = prev_sep / (prev_sep + (cfg.L - cur_sep))
            t_cross = prev_t + frac * (t - prev_t)
            return TrueChainResult(float(t_cross), n, h, 1)
        prev_sep = cur_sep
        prev_t = t

    raise RuntimeError(f"no return event within {max_steps} steps")
