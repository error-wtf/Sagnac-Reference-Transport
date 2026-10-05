#!/usr/bin/env python3
"""Falsifier battery: each intentional corruption MUST be detected.

Matrix output: artifacts/negative_controls.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sagnac_reference.analytic import (
    SagnacConfig, directed_times, delta_t_axle, delta_tau_detector,
    series_partial, normalized_asymmetry,
)
from sagnac_reference.inversion import velocity_from_delta_t, velocity_from_times
from sagnac_reference.phase import detector_phase

BETA = 0.5
cfg = SagnacConfig(L=1.0, c=1.0, v=BETA)
tp, tm = directed_times(cfg)
dt = delta_t_axle(cfg)
tau = delta_tau_detector(cfg)
omega = 2.0

results = []


def record(name, corrupted_value, clean_value, triggered, tolerance=None):
    if isinstance(corrupted_value, (int, float)) and isinstance(clean_value, (int, float)):
        rel = (abs(corrupted_value - clean_value) / max(abs(clean_value), 1e-300)
               if clean_value != 0 else abs(corrupted_value - clean_value))
    else:
        rel = None
    results.append({
        "corruption": name,
        "clean_value": clean_value,
        "corrupted_value": corrupted_value,
        "detected": bool(triggered),
        "relative_deviation": rel,
        "tolerance": tolerance,
    })


# 1. wrong branch c-v <-> c+v
wrong_tp = cfg.L / (cfg.c + cfg.v)
record("wrong_branch_swap_plus", wrong_tp, tp, abs(wrong_tp - tp) > 1e-3)

# 2. wrong sign in the alternating series (-beta -> +beta): the counter-
#    propagating branch then returns t_+ instead of t_- (branch swap at the
#    series level), which the harness must notice.
total_s, term_s = 0.0, 1.0
for _ in range(61):
    total_s += term_s
    term_s *= cfg.beta           # WRONG: must be -beta
t_minus_wrong_sign = cfg.T0 * total_s
record("wrong_series_sign", t_minus_wrong_sign, tm,
       abs(t_minus_wrong_sign - tp) < 1e-12 and abs(t_minus_wrong_sign - tm) > 0.1)

# 3. stationary detector (v=0 instead of v) — dt collapses to 0
dt_stationary = delta_t_axle(SagnacConfig(cfg.L, cfg.c, 0.0))
record("stationary_detector", dt_stationary, dt,
       abs(dt_stationary - dt) > 0.1)

# 4. both pulses same direction (t_- computed as t_+)
record("same_direction_both", tp, tm, abs(tp - tm) > 0.1)

# 5. wrong gamma factor: tau_wrong = gamma*dt (instead of dt/gamma)
tau_wrong = delta_t_axle(cfg) * cfg.gamma
record("wrong_gamma_multiplication", tau_wrong, tau,
       abs(tau_wrong - tau) > 0.05 * abs(tau))

# 6. phase from axle time (frame swap) at high beta
cfg_hi = SagnacConfig(L=1.0, c=1.0, v=0.8)
phi_wrong = omega * delta_t_axle(cfg_hi)
phi_right = detector_phase(cfg_hi, omega)
record("phase_from_axle_time", phi_wrong, phi_right,
       abs(phi_wrong - phi_right) / abs(phi_right) > 0.3)

# 7. wrong inversion root (other quadratic root: unstable / superluminal)
def invert_wrong_root(dt, L, c):
    root = (L * L + (c * dt) ** 2) ** 0.5
    return (c * c * dt) / (L - root)  # wrong sign of root -> blow-up
v_wrong = invert_wrong_root(dt, cfg.L, cfg.c)
record("inversion_wrong_root", v_wrong, cfg.v, abs(v_wrong) > 1.0 or abs(v_wrong - cfg.v) > 0.5)

# 8. alternating series read as physical velocity reversal:
# a naive "reversal" model would give t_- = L/(c-v) (identical to t_+) -> parity dies
naive_parity = (tp - tp) / (tp + tp)
record("alternation_as_reversal_parity_dies", naive_parity, cfg.beta,
       abs(naive_parity - cfg.beta) > 0.4)

# 9. domain guard fail-closed
try:
    SagnacConfig(L=1.0, c=1.0, v=1.0)
    record("domain_guard_v_ge_c", "no-raise", "ValueError", False)
except ValueError:
    record("domain_guard_v_ge_c", "ValueError", "ValueError", True)

# 10. normalized asymmetry identity (t+-t-)/(t++t-) = v/c — corruption check with wrong exponent
bad_asym = ((tp - tm) / (tp + tm)) ** 2
record("asymmetry_squared_detected", bad_asym, cfg.beta, abs(bad_asym - cfg.beta) > 0.05)

all_detected = all(r["detected"] for r in results)
out = {
    "verdict": "ALL_NEGATIVE_CONTROLS_DETECTED" if all_detected else "SOME_CONTROLS_MISSED",
    "n_controls": len(results),
    "n_detected": sum(r["detected"] for r in results),
    "controls": results,
}
p = Path("artifacts/negative_controls.json")
p.parent.mkdir(parents=True, exist_ok=True)
p.write_text(json.dumps(out, indent=1))
print(out["verdict"], f"({out['n_detected']}/{out['n_controls']} detected)")
for r in results:
    print(f"  {'OK ' if r['detected'] else 'MISS'} {r['corruption']}")
sys.exit(0 if all_detected else 1)
