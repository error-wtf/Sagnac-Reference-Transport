from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

from .analytic import (
    SagnacConfig,
    delta_t_axle,
    delta_t_axle_closed,
    delta_tau_detector,
    directed_times,
    normalized_asymmetry,
    series_partial,
)
from .inversion import velocity_from_delta_t, velocity_from_times
from .phase import detector_phase
from .segment_chain import return_time_segment_chain
from .transport_pde import return_time_pde
from .transport_pde_upwind import return_time_pde_upwind
from .true_chain import return_time_true_chain


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
    up_p = return_time_pde_upwind(cfg, +1, nx=pde_nx).return_time
    up_m = return_time_pde_upwind(cfg, -1, nx=pde_nx).return_time

    tolerances = {
        "analytic": 1e-12,
        "chain_series": 2e-12,
        "chain_true": 1e-10,
        "pde": 5e-3,
        "pde_upwind": 2e-4,
        "inversion": 1e-10,
    }
    gates = {
        "SAG-S1-forward-exact": (
            abs(delta_t_axle(cfg) - delta_t_axle_closed(cfg)) < tolerances["analytic"]),
        "SAG-S2-delta-exact": (
            abs((tp - tm) - delta_t_axle_closed(cfg)) < tolerances["analytic"]),
        "SAG-S3-parity": (
            abs(normalized_asymmetry(cfg) - cfg.beta) < tolerances["analytic"]),
        "SAG-S4-series-convergence": (
            abs(series_partial(cfg, +1, 50) - tp) < 1e-10
            and abs(series_partial(cfg, -1, 50) - tm) < 1e-10),
        "SAG-S5-continuum-chain": (
            abs(chain_p - tp) < tolerances["chain_series"]
            and abs(chain_m - tm) < tolerances["chain_series"]),
        "SAG-S6-proper-time": (
            abs(delta_tau_detector(cfg) - delta_t_axle(cfg) / cfg.gamma)
            < tolerances["analytic"]),
        "SAG-S7-phase-readout": (
            abs(detector_phase(cfg, omega) - omega * delta_tau_detector(cfg))
            < tolerances["analytic"]),
        "SAG-S8-inversion": (
            abs(velocity_from_times(tp, tm, cfg.c) - cfg.v) < tolerances["inversion"]
            and abs(velocity_from_delta_t(delta_t_axle(cfg), cfg.L, cfg.c) - cfg.v)
            < tolerances["inversion"]),
        "SAG-S9-pde-independent": (
            abs(pde_p - tp) < tolerances["pde"]
            and abs(pde_m - tm) < tolerances["pde"]
            and abs(up_p - tp) < tolerances["pde_upwind"]
            and abs(up_m - tm) < tolerances["pde_upwind"]),
        "SAG-S10-true-chain-independent": (
            abs(true_p - tp) < tolerances["chain_true"]
            and abs(true_m - tm) < tolerances["chain_true"]),
    }
    gate_meta = {
        "SAG-S1-forward-exact": {
            "statement": "delta_t(axle) computed per-direction equals the closed form",
            "formula": "t_+ - t_- == 2Lv/(c^2-v^2)",
            "implementation": "analytic.directed_times vs analytic.delta_t_axle_closed",
            "depends_on": [],
        },
        "SAG-S2-delta-exact": {
            "statement": "per-direction return times satisfy the closed form",
            "formula": "t_pm = L/(c -/+ v)",
            "implementation": "analytic.directed_times",
            "depends_on": [],
        },
        "SAG-S3-parity": {
            "statement": "normalized asymmetry identity",
            "formula": "(t_+ - t_-)/(t_+ + t_-) == v/c",
            "implementation": "analytic.normalized_asymmetry",
            "depends_on": [],
        },
        "SAG-S4-series-convergence": {
            "statement": "geometric correction series converges to the exact times",
            "formula": "t_pm^(N) = T0 sum_{n<=N} (+/-beta)^n",
            "implementation": "analytic.series_partial",
            "depends_on": [],
        },
        "SAG-S5-continuum-chain": {
            "statement": "baseline segment chain (relative-coordinate form) matches oracle",
            "formula": "sum_i w_i/(1 - s*beta) == t_s",
            "implementation": "segment_chain.return_time_segment_chain",
            "depends_on": ["SAG-S1"],
        },
        "SAG-S6-proper-time": {
            "statement": "detector proper time = axle time / gamma",
            "formula": "delta_tau_det = delta_t_axle * sqrt(1 - beta^2)",
            "implementation": "analytic.delta_tau_detector",
            "depends_on": [],
        },
        "SAG-S7-phase-readout": {
            "statement": "phase is readout AFTER transport: dphi = omega_det * tau_det",
            "formula": "delta_phi = omega_det * delta_tau_det",
            "implementation": "phase.detector_phase",
            "depends_on": ["SAG-S6"],
        },
        "SAG-S8-inversion": {
            "statement": "forward/inverse closure: v recovered from (t+,t-) and from dt",
            "formula": "v = c(t_+-t_-)/(t_++t_-); v = c^2 dt/(L + sqrt(L^2 + c^2 dt^2))",
            "implementation": "inversion.velocity_from_times, "
                              "inversion.velocity_from_delta_t",
            "depends_on": ["SAG-S1"],
        },
        "SAG-S9-pde-independent": {
            "statement": "two INDEPENDENT PDE schemes reproduce the oracle times",
            "formula": "u_t + s*c*u_x = 0 with moving detector x_D = v t mod L",
            "implementation": "transport_pde (Lax-Wendroff), "
                              "transport_pde_upwind (1st-order upwind)",
            "depends_on": ["SAG-S1"],
        },
        "SAG-S10-true-chain-independent": {
            "statement": "genuinely independent discrete ring simulation matches oracle",
            "formula": "x_light(t+h)=(x_light + s*c*h) mod L; "
                       "x_det(t+h)=(x_det + v*h) mod L; wrapped crossing",
            "implementation": "true_chain.return_time_true_chain",
            "depends_on": ["SAG-S1"],
        },
    }

    def _tolerance_for(gate_name: str) -> object:
        if "S9" in gate_name:
            return [tolerances["pde"], tolerances["pde_upwind"]]
        if "S10" in gate_name:
            return tolerances["chain_true"]
        if "S5" in gate_name:
            return tolerances["chain_series"]
        if "S8" in gate_name:
            return tolerances["inversion"]
        return tolerances["analytic"]

    gates_detailed = {
        name: {
            "pass": pass_value,
            **gate_meta.get(name, {}),
            "tolerance": _tolerance_for(name),
        }
        for name, pass_value in gates.items()
    }
    dependency_graph = {
        "SAG-S1": [], "SAG-S2": [], "SAG-S3": [], "SAG-S4": [],
        "SAG-S5": ["SAG-S1"], "SAG-S6": [], "SAG-S7": ["SAG-S6"],
        "SAG-S8": ["SAG-S1"], "SAG-S9": ["SAG-S1"], "SAG-S10": ["SAG-S1"],
    }
    return {
        "status": ("SAGNAC_REFERENCE_CLOSURE_PASS" if all(gates.values())
                    else "SAGNAC_REFERENCE_CLOSURE_FAIL"),
        "config": asdict(cfg),
        "tolerances": tolerances,
        "gates": gates,
        "gates_detailed": gates_detailed,
        "dependency_graph": dependency_graph,
        "measurements": {
            "analytic": {"t_plus": tp, "t_minus": tm, "delta_t": delta_t_axle(cfg)},
            "chain_series": {"t_plus": chain_p, "t_minus": chain_m},
            "chain_true": {"t_plus": true_p, "t_minus": true_m},
            "pde": {"t_plus": pde_p, "t_minus": pde_m},
            "pde_upwind": {"t_plus": up_p, "t_minus": up_m},
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
