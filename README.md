# Sagnac Reference Transport

Independent, analytically anchored validation harness for transport methods using the idealised circular Sagnac problem in flat spacetime.

This repository is intentionally **SSZ-independent**. Its purpose is to validate a transport architecture against a known reference problem before any later bridge to SSZ or another geometry model is attempted.

## Normative chain

\[
\text{paper reference}
\to \text{analytic oracle}
\leftrightarrow \text{discrete segment chain}
\leftrightarrow \text{transport PDE}
\to \text{proper time}
\to \text{phase readout}
\to \text{inverse reconstruction}.
\]

Core reference formulas:

\[
t_+=\frac{L}{c-v},\qquad t_-=\frac{L}{c+v},
\]

\[
\Delta t_{\rm axle}=\frac{2Lv}{c^2-v^2},
\]

\[
\Delta\tau_{\rm det}=\Delta t_{\rm axle}\sqrt{1-v^2/c^2},
\]

\[
\Delta\phi=\omega_{\rm det}\Delta\tau_{\rm det}.
\]

The geometric correction series is treated as a convergence construction, not as a physical back-and-forth oscillation of light.

## Repository layout

```text
src/sagnac_reference/
    analytic.py          exact oracle and frame conversion
    segment_chain.py     discrete moving-closure transport
    transport_pde.py     periodic finite-difference advection reference
    phase.py             detector phase readout
    inversion.py         reconstruct v from measured observables
    diagnostics.py       gate evaluation helpers

tests/
    test_analytic.py
    test_segment_chain.py
    test_pde.py
    test_phase_inversion.py
    test_negative_controls.py

docs/SAGNAC_REFERENCE_SPEC.md
reference/Sagnac.pdf
artifacts/
tools/evaluate_reference.py
```

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
pytest -q
python tools/evaluate_reference.py
```

The evaluator writes `artifacts/SAGNAC_REFERENCE_VERDICT.json`.

## Scientific scope

The baseline model assumes a prescribed circular optical path in flat spacetime, constant `c`, uniform detector/closure speed `v`, and `|v| < c`. It does not invoke an aether, a varying local light speed, frame dragging, or a physically rotating spacetime as the ordinary Sagnac mechanism. Those belong to separate extension layers.

## Status

`SAGNAC_REFERENCE_CLOSURE_PASS` — 10/10 gates (SAG-S1..S10) · 74/74 tests · 10/10 negative controls detected · CI green on Python 3.10/3.12/3.14 (ruff gate included).

**Hardening completed 2026-10-05** (convergence study, phase-4/5 rewrite):
- Convergence campaign (`artifacts/convergence/`) over 13 beta values: Lax-Wendroff route shows observed order ~2 after freezing the pulse width (the baseline `sigma = 6*dx` scaled with the grid and destroyed convergence — root cause documented in `transport_pde.py`).
- Frozen tolerances: PDE 5e-3 (nx=2400, ~12x campaign margin), true chain 1e-10 (64k steps, campaign max 9.6e-12).
- NEW `true_chain.py`: genuinely independent discrete ring simulation (light and detector positions evolved separately; return located by wrapped-crossing interpolation; no relative-speed formula inside).
- NEW negative-control battery: 10/10 intentional corruptions detected (`tools/run_negative_controls.py`, `artifacts/negative_controls.json`).
- Explorative phase-geometry block (`tools/run_phase_geometry.py`): the alternating correction series as damped phasor dynamics — labeled `MATHEMATICAL_CONVERGENCE_GEOMETRY / NOT_PHYSICAL_LIGHT_OSCILLATION`.
- Second independent PDE route: `transport_pde_upwind.py` (first-order upwind, observed order ~1) — gate SAG-S9 requires BOTH schemes to reproduce the oracle.
- Beta-scan cross-route regression (`tests/test_beta_scan_regression.py`): 8 betas x 4 routes against frozen tolerances.
- 50-digit `decimal` cross-check of the closed form, series remainder bounds, small-v limit, direction-swap identity (`tests/test_oracle_precision.py`).
- Dimensional scaling contracts: `L -> lambda*L` linear time scaling, unit-rescale invariance (`tests/test_scaling_invariance.py`).
- Gates carry machine-readable provenance in the verdict: statement, formula, implementation, dependencies, tolerance per gate (`gates_detailed` + `dependency_graph` in `artifacts/SAGNAC_REFERENCE_VERDICT.json`).

Commands:

```bash
pytest -q
python tools/evaluate_reference.py
python tools/run_convergence_campaign.py   # regenerate artifacts/convergence
python tools/run_negative_controls.py
python tools/run_phase_geometry.py
```

A future bridge to SSZ is documented as a contract only (`docs/TRANSPORT_CONTRACT.md`) and is deliberately not implemented.
