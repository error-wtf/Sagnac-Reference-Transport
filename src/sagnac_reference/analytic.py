from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class SagnacConfig:
    L: float
    c: float
    v: float

    def __post_init__(self) -> None:
        if not math.isfinite(self.L) or self.L <= 0:
            raise ValueError("L must be finite and > 0")
        if not math.isfinite(self.c) or self.c <= 0:
            raise ValueError("c must be finite and > 0")
        if not math.isfinite(self.v) or abs(self.v) >= self.c:
            raise ValueError("v must be finite and satisfy |v| < c")

    @property
    def beta(self) -> float:
        return self.v / self.c

    @property
    def gamma(self) -> float:
        return 1.0 / math.sqrt(1.0 - self.beta * self.beta)

    @property
    def T0(self) -> float:
        return self.L / self.c


def directed_times(cfg: SagnacConfig) -> tuple[float, float]:
    return cfg.L / (cfg.c - cfg.v), cfg.L / (cfg.c + cfg.v)


def delta_t_axle(cfg: SagnacConfig) -> float:
    tp, tm = directed_times(cfg)
    return tp - tm


def delta_t_axle_closed(cfg: SagnacConfig) -> float:
    return 2.0 * cfg.L * cfg.v / (cfg.c * cfg.c - cfg.v * cfg.v)


def normalized_asymmetry(cfg: SagnacConfig) -> float:
    tp, tm = directed_times(cfg)
    return (tp - tm) / (tp + tm)


def delta_tau_detector(cfg: SagnacConfig) -> float:
    return delta_t_axle(cfg) / cfg.gamma


def series_partial(cfg: SagnacConfig, direction: int, N: int) -> float:
    if direction not in (+1, -1):
        raise ValueError("direction must be +1 (co-rotating) or -1 (counter-propagating)")
    if N < 0:
        raise ValueError("N must be >= 0")
    q = cfg.beta if direction == +1 else -cfg.beta
    total = 0.0
    term = 1.0
    for _ in range(N + 1):
        total += term
        term *= q
    return cfg.T0 * total


def series_remainder_bound(cfg: SagnacConfig, direction: int, N: int) -> float:
    q = cfg.beta if direction == +1 else -cfg.beta
    return abs(cfg.T0 * (q ** (N + 1)) / (1.0 - q))
