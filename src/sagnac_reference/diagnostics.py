from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from .analytic import (
    SagnacConfig,
    directed_times,
    delta_t_axle,
    delta_t_axle_closed,
    normalized_asymmetry,
    delta_tau_detector,
    series_partial,
)
from .segment_chain import return_time_segment_chain
from .true_chain import return_time_true_chain
from .transport_pde import return_time_pde
from .phase import detector_phase
from .inversion import velocity_from_times, velocity_from_delta_t


def evaluate(cfg: SagnacConfig, omega: float = 2.0, *, pde_nx: int = 2400) -> dict:
    """Evaluate all SAG-S gates. Tolerances FROZEN from the convergence
    campaign 2026-10-05 (artifacts/convergence/convergence_campaign.{csv,json}):
      * PDE at nx=2400: max observed error 4.225e-04 over the full beta matrix
        (13 values, |beta| up to 0.8) -> gate 5e-3 with ~12x margin.
      * true chain at 64k steps: max observed error 9.6e-12 -> gate 1e-10.
    Frozen after the campaign; do not loosen without a new campaign.
    """
    tp, tm = directed_times(cfg)
    chain_p = return_time_segment_chain(cfg, +1, 4000).return_time
    chain_m = return_time_segment_chain(cfg, -1, 4000).return_time
    true_p = return_time_true_chain(cfg, +1, steps=32000).return_time
    true_m = return_time_true_chain(cfg, -1, steps=32000).return_time
    pde_p = return_time_pde(cfg, +1, nx=pde_nx).return_time
    pde_m = return_time_pde(cfg, -1, nx=pde_nx).return_time

    tolerances = {
        "analytic": 1e-12,
        "chain_series": 2e-12,
        "chain_true": 1e-10,
        "pde": 5e-3,
        "inversion": 1e-10,
    }
    gates = {
        "SAG-S1-forward-exact": abs(delta_t_axle(cfg) - delta_t_axle_closed(cfg)) < tolerances["analytic"],
        "SAG-S2-delta-exact": abs((tp - tm) - delta_t_axle_closed(cfg)) < tolerances["analytic"],
        "SAG-S3-parity": abs(normalized_asymmetry(cfg) - cfg.beta) < tolerances["analytic"],
        "SAG-S4-series-convergence": abs(series_partial(cfg, +1, 50) - tp) < 1e-10 and abs(series_partial(cfg, -1, 50) - tm) < 1e-10,
        "SAG-S5-continuum-chain": abs(chain_p - tp) < tolerances["chain_series"] and abs(chain_m - tm) < tolerances["chain_series"],
        "SAG-S6-proper-time": abs(delta_tau_detector(cfg) - delta_t_axle(cfg) / cfg.gamma) < tolerances["analytic"],
        "SAG-S7-phase-readout": abs(detector_phase(cfg, omega) - omega * delta_tau_detector(cfg)) < tolerances["analytic"],
        "SAG-S8-inversion": abs(velocity_from_times(tp, tm, cfg.c) - cfg.v) < tolerances["inversion"] and abs(velocity_from_delta_t(delta_t_axle(cfg), cfg.L, cfg.c) - cfg.v) < tolerances["inversion"],
        "SAG-S9-pde-independent": abs(pde_p - tp) < tolerances["pde"] and abs(pde_m - tm) < tolerances["pde"],
        "SAG-S10-true-chain-independent": abs(true_p - tp) < tolerances["chain_true"] and abs(true_m - tm) < tolerances["chain_true"],
    }
    return {
        "status": "SAGNAC_REFERENCE_CLOSURE_PASS" if all(gates.values()) else "SAGNAC_REFERENCE_CLOSURE_FAIL",
        "config": asdict(cfg),
        "tolerances": tolerances,
        "gates": gates,
        "measurements": {
            "analytic": {"t_plus": tp, "t_minus": tm, "delta_t": delta_t_axle(cfg)},
            "chain_series": {"t_plus": chain_p, "t_minus": chain_m},
            "chain_true": {"t_plus": true_p, "t_minus": true_m},
            "pde": {"t_plus": pde_p, "t_minus": pde_m},
            "proper_time_delta": delta_tau_detector(cfg),
            "phase": detector_phase(cfg, omega),
        },
    }


def write_verdict(path: str | Path, cfg: SagnacConfig | None = None) -> dict:
    cfg = cfg or SagnacConfig(L=1.0, c=1.0, v=0.2)
    out = evaluate(cfg)
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return out
