"""Second independent PDE route: first-order upwind finite volume.

Independent of the Lax-Wendroff route (different update rule, different
dispersion/dissipation character).  First-order in space/time: observed
order ~1.  Used as a cross-check that BOTH numerics agree on the same
transport observable (SAG-S9 strengthening).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .analytic import SagnacConfig
from .transport_pde import _periodic_interp, _wrapped_distance


@dataclass(frozen=True)
class UpwindResult:
    return_time: float
    nx: int
    dt: float
    cfl: float


def return_time_pde_upwind(
    cfg: SagnacConfig,
    direction: int,
    *,
    nx: int = 2400,
    cfl: float = 0.45,
    pulse_sigma: float | None = None,
) -> UpwindResult:
    if direction not in (+1, -1):
        raise ValueError("direction must be +1 or -1")
    if nx < 200:
        raise ValueError("nx must be >= 200")
    if not (0.05 <= cfl < 1.0):
        raise ValueError("cfl must be in [0.05, 1)")

    a = direction * cfg.c
    dx = cfg.L / nx
    dt = cfl * dx / abs(a)
    sigma = pulse_sigma or 0.04 * cfg.L
    if sigma < 8.0 * dx:
        raise ValueError(
            f"pulse sigma {sigma:.3e} too narrow for nx={nx} (needs >= 8*dx)")
    x = np.arange(nx) * dx
    u = np.exp(-0.5 * (_wrapped_distance(x, 0.0, cfg.L) / sigma) ** 2)
    lam = a * dt / dx

    rel_speed = abs(a - cfg.v)
    max_steps = int(np.ceil(1.05 * cfg.L / rel_speed / dt)) + 2
    exclude_until = max(8.0 * sigma / abs(a), 8.0 * dt)

    ts: list[float] = []
    ys: list[float] = []
    for n in range(1, max_steps + 1):
        up = np.roll(u, -1)
        um = np.roll(u, 1)
        if a > 0:
            u = u - lam * (u - um)          # upwind from the left
        else:
            u = u - lam * (up - u)          # upwind from the right
        t = n * dt
        if t < exclude_until:
            continue
        det_x = cfg.v * t
        y = _periodic_interp(u, det_x, cfg.L)
        ts.append(t)
        ys.append(y)
        if len(ys) >= 3:
            y0, y1, y2 = ys[-3:]
            if y1 > y0 and y1 >= y2 and y1 > 0.15:
                den = y0 - 2.0 * y1 + y2
                off = 0.0 if abs(den) < 1e-15 else 0.5 * (y0 - y2) / den
                off = float(np.clip(off, -1.0, 1.0))
                return UpwindResult(ts[-2] + off * dt, nx, dt, cfl)
    raise RuntimeError("upwind PDE return peak not detected")
