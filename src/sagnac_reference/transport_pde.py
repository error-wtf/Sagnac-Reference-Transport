from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .analytic import SagnacConfig


@dataclass(frozen=True)
class PDEResult:
    return_time: float
    nx: int
    dt: float
    cfl: float
    peak_value: float


def _periodic_interp(field: np.ndarray, x: float, L: float) -> float:
    nx = field.size
    y = (x % L) / L * nx
    i0 = int(np.floor(y)) % nx
    frac = y - np.floor(y)
    i1 = (i0 + 1) % nx
    return float((1.0 - frac) * field[i0] + frac * field[i1])


def _wrapped_distance(x: np.ndarray, x0: float, L: float) -> np.ndarray:
    return (x - x0 + 0.5 * L) % L - 0.5 * L


def return_time_pde(
    cfg: SagnacConfig,
    direction: int,
    *,
    nx: int = 2400,
    cfl: float = 0.55,
    pulse_sigma: float | None = None,
) -> PDEResult:
    """Finite-difference periodic advection with a moving detector.

    Uses Lax-Wendroff for u_t + a u_x = 0 with a=±c and detects the first
    non-trivial detector peak. The analytic return-time formula is not used by
    the solver.
    """
    if direction not in (+1, -1):
        raise ValueError("direction must be +1 or -1")
    if nx < 200:
        raise ValueError("nx must be >= 200")
    if not (0.05 <= cfl < 1.0):
        raise ValueError("cfl must be in [0.05, 1)")

    dx = cfg.L / nx
    dt = cfl * dx / cfg.c
    # Pulse width MUST be grid-independent (convergence study 2026-10-05):
    # sigma = 6*dx made the parabolic peak-interpolation error dx^2/sigma^2 =
    # const, destroying convergence (observed order ~0.5 instead of 2).
    sigma = pulse_sigma or 0.04 * cfg.L
    if sigma < 8.0 * dx:
        raise ValueError(
            f"pulse sigma {sigma:.3e} too narrow for nx={nx} "
            f"(needs sigma >= 8*dx = {8.0*dx:.3e}); reduce nx or raise sigma")
    x = np.arange(nx) * dx
    u = np.exp(-0.5 * (_wrapped_distance(x, 0.0, cfg.L) / sigma) ** 2)
    a = direction * cfg.c
    lam = a * dt / dx

    # Conservative upper horizon derived only from domain bounds |v|<c.
    # 4 L/c is enough for moderate validation cases; for relativistic v,
    # scale by the minimum guaranteed relative speed c-|v|.
    max_time = 1.20 * cfg.L / max(cfg.c - abs(cfg.v), 1e-15)
    max_steps = int(np.ceil(max_time / dt)) + 2
    exclude_until = max(8.0 * sigma / cfg.c, 8.0 * dt)

    samples_t: list[float] = []
    samples_y: list[float] = []

    for n in range(1, max_steps + 1):
        up = np.roll(u, -1)
        um = np.roll(u, 1)
        u = u - 0.5 * lam * (up - um) + 0.5 * lam * lam * (up - 2.0 * u + um)
        t = n * dt
        det_x = cfg.v * t
        y = _periodic_interp(u, det_x, cfg.L)
        if t >= exclude_until:
            samples_t.append(t)
            samples_y.append(y)
            if len(samples_y) >= 3:
                a0, a1, a2 = samples_y[-3:]
                if a1 > a0 and a1 >= a2 and a1 > 0.15:
                    # Quadratic interpolation of the sampled detector peak.
                    y0, y1, y2 = a0, a1, a2
                    denom = (y0 - 2.0 * y1 + y2)
                    offset = 0.0 if abs(denom) < 1e-15 else 0.5 * (y0 - y2) / denom
                    offset = float(np.clip(offset, -1.0, 1.0))
                    t_peak = samples_t[-2] + offset * dt
                    return PDEResult(t_peak, nx, dt, cfl, float(y1))

    raise RuntimeError("PDE return peak not detected")
