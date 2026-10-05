# Hermes handoff

The repository is a deliberately small, independent Sagnac calibration harness. Harden it without coupling it to SSZ.

Key non-negotiables:

1. Preserve analytic/PDE/segment-chain independence.
2. Do not import SSZ metric code, JIF internals, frozen members, or SSZ closure verdicts.
3. Preserve axle-frame coordinate time vs rotating-detector proper time.
4. Preserve phase as readout after transport, not as the transport mechanism.
5. Treat the alternating correction series as mathematical convergence structure, not physical reversal of light.
6. Keep negative controls and fail-closed domain guards.
7. Freeze tolerances only after a documented convergence study.
8. Any future SSZ bridge must be a separate package/layer and must not alter this repository's reference oracle.
