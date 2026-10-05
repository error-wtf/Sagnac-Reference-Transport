# SAGNAC_REFERENCE_SPEC.md

## Status

**Specification type:** independent external reference validation for transport architecture  
**Target repository:** standalone `Sagnac-Reference-Transport` reference repository  
**Layer:** reference transport / calibration layer  
**Relation to existing closure:** additive only; MUST NOT redefine or retroactively alter the meaning of the existing `TRUE_FULL_CLOSURE_PASS` verdict.

---

## 1. Purpose

This specification defines an analytically solvable Sagnac reference problem for validating a generic transport implementation before that implementation is compared with, abstracted toward, or reused by SSZ transport code.

The reference layer MUST remain independent of SSZ-specific assumptions, metric functions, frozen members, JIF internals, or any SSZ closure verdict.

The validation target is the **transport method**, not SSZ itself.

The required logic is:

\[
\boxed{
\text{paper reference}
\rightarrow
\text{analytic oracle}
\rightarrow
\text{PDE / discrete segment chain}
\rightarrow
\text{proper-time readout}
\rightarrow
\text{phase readout}
\rightarrow
\text{inversion}
}
\]

Only after this independent block passes may a separate downstream bridge compare or unify its abstract transport architecture with SSZ or another transport implementation.

---

## 2. Normative physical scope

The reference problem is the idealized circular Sagnac configuration in flat spacetime:

- prescribed circular optical path of radius \(R\),
- path length \(L = 2\pi R\),
- uniformly rotating detector/closure point with angular speed \(\Omega\),
- tangential detector speed \(v = \Omega R\),
- constant local light speed \(c\),
- no aether,
- no refractive-medium mechanism,
- no local variation of \(c\),
- no frame dragging,
- no physically rotating spacetime invoked as the cause of the ordinary Sagnac signal.

Any gravitational or background-rotational contribution belongs to a separate extension layer and MUST NOT be folded into this reference oracle.

Domain:

\[
L>0,\qquad c>0,\qquad |v|<c.
\]

Define

\[
\beta = \frac{v}{c},\qquad
\gamma = \frac{1}{\sqrt{1-\beta^2}}.
\]

---

## 3. Normative analytic oracle

### 3.1 Directed return times in the inertial axle frame

\[
t_+ = \frac{L}{c-v},
\qquad
 t_- = \frac{L}{c+v}.
\]

The co-rotating branch uses the receding closure speed \(c-v\); the counter-propagating branch uses the approaching closure speed \(c+v\).

### 3.2 Axle-frame return-time asymmetry

\[
\Delta t_{\rm axle}
= t_+ - t_-
= \frac{2Lv}{c^2-v^2}.
\]

For a circular path,

\[
L=2\pi R,\qquad v=\Omega R,\qquad A=\pi R^2,
\]

so

\[
\Delta t_{\rm axle}
= \frac{4A\Omega}{c^2-\Omega^2R^2}
= \frac{4A\Omega}{c^2}\gamma^2.
\]

Low-velocity limit:

\[
\Delta t_{\rm axle}
\approx \frac{4A\Omega}{c^2}
\qquad (|\beta|\ll 1).
\]

### 3.3 Even/odd decomposition under direction reversal

Let

\[
T_0=\frac{L}{c}.
\]

The two geometric correction series are

\[
t_+
= T_0\sum_{n=0}^{\infty}\beta^n,
\qquad
 t_-
= T_0\sum_{n=0}^{\infty}(-\beta)^n.
\]

Their sum isolates the even part:

\[
t_+ + t_-
= 2T_0(1+\beta^2+\beta^4+\cdots),
\]

while their difference isolates the odd, direction-dependent part:

\[
t_+ - t_-
= 2T_0(\beta+\beta^3+\beta^5+\cdots).
\]

The normalized asymmetry is exactly

\[
\boxed{
\frac{t_+-t_-}{t_++t_-}=\beta=\frac{v}{c}
}.
\]

### 3.4 Proper-time readout of the co-rotating detector

The detector proper time is

\[
d\tau_{\rm det}
= dt_{\rm axle}\sqrt{1-\beta^2}
= \frac{dt_{\rm axle}}{\gamma}.
\]

Therefore

\[
\Delta\tau_{\rm det}
= \frac{\Delta t_{\rm axle}}{\gamma}
= \frac{4A\Omega}{c^2}\gamma.
\]

The implementation MUST keep

\[
\Delta t_{\rm axle}
\neq
\Delta\tau_{\rm det}
\]

when relativistic higher-order corrections are retained.

### 3.5 Phase readout

Phase is a readout of the already-derived return-time asymmetry, not the mechanism represented by the geometric correction series.

For detector-referred angular frequency \(\omega_{\rm det}\),

\[
\boxed{
\Delta\phi
= \omega_{\rm det}\,\Delta\tau_{\rm det}
}.
\]

The reference implementation MUST derive transport first, then proper time, then phase.

---

## 4. Required implementation layout

Recommended isolated structure:

```text
src/sagnac_reference/
    __init__.py
    sagnac.py
    transport_pde.py
    segment_chain.py
    phase_readout.py
    inversion.py
    diagnostics.py

 tests/
    test_sagnac_forward_exact.py
    test_sagnac_delta_exact.py
    test_sagnac_parity.py
    test_sagnac_series_convergence.py
    test_sagnac_continuum.py
    test_sagnac_proper_time.py
    test_sagnac_phase_readout.py
    test_sagnac_inversion.py
    test_sagnac_falsifiers.py

artifacts/reference_transport/
    sagnac_validation.json
    sagnac_dependency_graph.json

 tools/
    evaluate_reference_transport.py
```

The following existing SSZ files MUST NOT be modified as part of the initial reference implementation:

```text
src/ssz_p5/postclosure/transport.py
artifacts/true_full_closure/physics_dependency_graph.json
tools/evaluate_true_closure.py
```

---

## 5. Independent numerical routes

Three routes MUST be available and compared.

### Route A — analytic oracle

Direct evaluation of the closed forms in Section 3.

This is the reference truth for this validation block.

### Route B — transport PDE

Use a periodic one-dimensional coordinate \(x\in[0,L)\) in the inertial axle frame.

Transport equations:

\[
\partial_t u_+ + c\,\partial_xu_+ = 0,
\]

\[
\partial_t u_- - c\,\partial_xu_- = 0.
\]

Periodic boundary condition:

\[
u_\pm(x+L,t)=u_\pm(x,t).
\]

Detector worldline:

\[
x_D(t)=vt\pmod L.
\]

The numerical return time is the first nontrivial recurrence of the transported pulse at the moving detector.

The production PDE test MUST NOT determine the return time by substituting the closed-form \(L/(c\mp v)\) into the event logic.

A conservative second-order periodic discretization or equivalent method-of-lines scheme is recommended. Convergence MUST be assessed under grid refinement.

### Route C — discrete segment chain

A spatially segmented ring uses

\[
\Delta x=\frac{L}{N}.
\]

The signal advances according to its transport direction while the detector advances simultaneously according to \(v\). A return is detected through event/crossing interpolation.

The method MUST satisfy

\[
t_\pm^{\rm chain}(N)\to t_\pm
\qquad (N\to\infty).
\]

A separate geometric-series implementation MAY also be included as an exact convergence control:

\[
t_+^{(N)}=T_0\sum_{n=0}^{N}\beta^n,
\qquad
 t_-^{(N)}=T_0\sum_{n=0}^{N}(-\beta)^n.
\]

The geometric-series path and the spatial segment-chain path MUST be labelled as distinct algorithms.

---

## 6. Gate matrix

Use an independent namespace. Recommended IDs: `SAG-S1` through `SAG-S9`.

### SAG-S1 — Forward exact

**Claim:** the directed analytic implementation reproduces

\[
t_+=\frac{L}{c-v},\qquad
 t_-=\frac{L}{c+v}.
\]

**Required checks:**

- multiple \(L,c,v\) tuples,
- positive and negative \(v\),
- \(v=0\),
- near-relativistic but domain-safe values such as \(|\beta|\le 0.9\).

**Pass criterion:** machine-precision agreement for the algebraic implementation.

---

### SAG-S2 — Delta exact

**Claim:**

\[
\Delta t_{\rm axle}
= t_+-t_-
= \frac{2Lv}{c^2-v^2}.
\]

**Pass criterion:** direct-route and difference-route agreement to floating-point precision for the oracle implementation; numerical PDE/chain agreement within their independently declared discretization tolerances.

---

### SAG-S3 — Parity and direction reversal

Under

\[
v\to -v,
\]

require

\[
t_+(v)=t_-(-v),
\]

\[
t_-(v)=t_+(-v),
\]

and

\[
\Delta t(-v)=-\Delta t(v).
\]

Also require

\[
\frac{t_+-t_-}{t_++t_-}=\frac{v}{c}.
\]

This gate is the primary sign-convention guard.

---

### SAG-S4 — Series convergence

For \(|\beta|<1\), require

\[
t_+^{(N)}\to \frac{L}{c-v},
\qquad
 t_-^{(N)}\to \frac{L}{c+v}.
\]

Use the exact remainders where convenient:

\[
R_+^{(N)}
= \frac{T_0\beta^{N+1}}{1-\beta},
\]

\[
R_-^{(N)}
= \frac{T_0(-\beta)^{N+1}}{1+\beta}.
\]

The alternating zig-zag of the counter-propagating partial sums MUST be documented as a convergence pattern in iteration index, not a physical oscillation of the light field.

---

### SAG-S5 — Continuum / loop-integral / segment-chain closure

Require agreement between:

\[
\int \frac{ds}{c-v}=\frac{L}{c-v},
\qquad
\int \frac{ds}{c+v}=\frac{L}{c+v},
\]

and the spatial segment-chain continuum limit.

At minimum:

\[
|t_\pm^{\rm chain}(2N)-t_\pm^{\rm exact}|
<
|t_\pm^{\rm chain}(N)-t_\pm^{\rm exact}|
\]

for a declared refinement ladder.

Recommended refinement ladder:

```text
N = 128, 256, 512, 1024, 2048
```

The gate SHOULD report the observed convergence order where meaningful.

---

### SAG-S6 — Proper-time readout

Require

\[
\Delta\tau_{\rm det}
= \Delta t_{\rm axle}\sqrt{1-\beta^2}
= \frac{\Delta t_{\rm axle}}{\gamma}
= \frac{4A\Omega}{c^2}\gamma.
\]

The implementation MUST preserve separate variable names and metadata for axle-frame coordinate time and detector proper time.

No automatic aliasing of these quantities is permitted.

---

### SAG-S7 — Phase readout

Require

\[
\Delta\phi
= \omega_{\rm det}\Delta\tau_{\rm det}.
\]

Dependency requirement:

```text
transport -> delta_t_axle -> delta_tau_det -> delta_phi
```

The phase module MUST NOT feed back into the return-time solver.

This gate verifies readout only.

---

### SAG-S8 — Inversion

Two inverse routes are required.

#### Route 1: from both return times

\[
\boxed{
 v
= c\frac{t_+-t_-}{t_++t_-}
}
\]

#### Route 2: from \(\Delta t\) only

Starting from

\[
\Delta t
=\frac{2Lv}{c^2-v^2},
\]

use the numerically stable positive/negative-sign-preserving form

\[
\boxed{
 v
=\frac{c^2\Delta t}
 {L+\sqrt{L^2+c^2\Delta t^2}}
}
\]

for the branch continuous through \(v=0\).

Require

\[
v_{\rm reconstructed}\approx v_{\rm input}.
\]

The inverse tests MUST include positive, negative, and near-zero velocities.

---

### SAG-S9 — Falsifiers and documentation guards

The harness MUST deliberately corrupt the reference implementation and confirm that the appropriate gates fail.

Minimum falsifier set:

1. **wrong co-rotating denominator**: use \(c+v\) instead of \(c-v\),
2. **wrong counter-propagating sign**,
3. **fixed detector**: replace \(x_D(t)=vt\) with \(x_D=0\),
4. **same propagation direction for both PDE branches**,
5. **wrong alternating-series sign**,
6. **wrong proper-time factor**: multiply by \(\gamma\) instead of divide by \(\gamma\), producing a spurious \(\gamma^3\)-type dependence when applied to the axle result,
7. **frame swap**: label \(\Delta t_{\rm axle}\) as \(\Delta\tau_{\rm det}\),
8. **wrong phase readout source**: compute detector phase from the wrong time variable,
9. **unstable/wrong inverse branch**,
10. **domain violation**: accept \(|v|\ge c\) instead of failing closed.

Documentation guard:

- wording that identifies the alternating correction-series zig-zag as a physical back-and-forth oscillation of light MUST fail lint/review checks.

The last item is a semantic/documentation falsifier, not a numerical physics test.

---

## 7. Tolerance policy

Tolerances MUST be method-specific and declared in machine-readable metadata.

### 7.1 Analytic oracle

Recommended default for float64 algebraic identities:

```text
rtol = 1e-12
atol = 1e-14 * characteristic_time
```

Higher precision MAY be used, but passing MUST NOT depend on arbitrary decimal overkill.

### 7.2 Geometric series

Prefer exact analytic remainder checks over a fixed tolerance whenever possible.

A pass requires the numerical error to be bounded by the known remainder plus floating-point slack.

### 7.3 PDE

Do not use a single unexplained absolute threshold.

Require both:

1. final refined-grid error below a declared tolerance, and
2. monotone or near-monotone convergence over the refinement ladder.

Suggested initial target after implementation tuning:

```text
relative return-time error <= 1e-5 on the finest registered grid
```

This is a starting engineering target, not a theorem. The final repository tolerance SHOULD be set from measured convergence and frozen only after reproducible benchmarking.

### 7.4 Segment chain

Suggested initial target:

```text
relative return-time error <= 1e-6 at N >= 2048
```

Again, the normative requirement is demonstrated convergence; the exact frozen threshold should follow measured implementation behavior.

### 7.5 Inversion

For oracle-derived inputs:

```text
relative velocity reconstruction error <= 1e-12
```

For PDE/chain-derived inputs, propagate the upstream transport tolerance rather than demanding oracle-level precision.

---

## 8. Machine-readable verdict

Recommended artifact:

```text
artifacts/reference_transport/sagnac_validation.json
```

Required top-level fields:

```json
{
  "status": "SAGNAC_REFERENCE_VALIDATION_PASS",
  "spec_version": "1.0",
  "scope": "flat_spacetime_circular_sagnac",
  "independent_of_ssz_member": true,
  "gates": {},
  "negative_controls": {},
  "tolerances": {},
  "provenance": {}
}
```

The evaluator MUST compute the verdict from gate results. The final status MUST NOT be manually hard-coded independent of evidence.

---

## 9. Dependency graph

Recommended graph:

```text
SAG-S1  Forward exact
  |
  +--> SAG-S2  Delta exact
  |      |
  |      +--> SAG-S3  parity / reversal
  |
  +--> SAG-S4  series convergence
  |      |
  |      +--> SAG-S5  continuum / chain / PDE closure
  |                |
  |                +--> SAG-S6  detector proper time
  |                          |
  |                          +--> SAG-S7  phase readout
  |                          |
  |                          +--> SAG-S8  inversion
  |
  +--------------------------------------> SAG-S9  falsifiers
```

Overall reference verdict:

\[
R_{\rm Sag}
=
\bigwedge_{i=1}^{9}\mathrm{SAG\text{-}S}i.
\]

Recommended string:

```text
SAGNAC_REFERENCE_VALIDATION_PASS
```

---

## 10. Relationship to existing SSZ closure

The existing SSZ closure verdict is historically and semantically frozen.

This new block MUST be additive:

\[
\boxed{
\text{existing TRUE FULL CLOSURE}
+
\text{independent Sagnac transport validation}
}
\]

Do not renumber, reinterpret, or retroactively make existing SSZ gates depend on the Sagnac block.

A later optional meta-gate may be introduced only after the reference block is complete and independently passing.

Example future meta-gate:

\[
G150
:
\text{REFERENCE-VALIDATED TRANSPORT ARCHITECTURE}
\]

with dependency

\[
G150
=
R_{\rm Sag}
\land
G130
\land
G140.
\]

Its meaning must be narrow:

> the transport architecture has passed an independent analytic Sagnac reference validation, while the registered SSZ forward chain and full closure retain their existing independent verdicts.

It MUST NOT be worded as empirical confirmation of SSZ by Sagnac.

---

## 11. Future bridge contract — specification only, not part of v1 implementation

After the reference block passes, define a generic transport contract of the form

\[
\mathcal S(x)
\rightarrow
A(x)
\rightarrow
\mathcal L[A]
\rightarrow
\Psi
\rightarrow
\mathcal O.
\]

### Sagnac instance

\[
\mathcal S_{\rm Sag}=(L,c,v),
\]

\[
\mathcal L_\pm
=\partial_t\pm c\partial_x,
\]

\[
\mathcal O_{\rm Sag}
=(t_+,t_-,\Delta t_{\rm axle},\Delta\tau_{\rm det},\Delta\phi).
\]

### SSZ instance

The existing SSZ transport side uses geometry-derived structure such as

\[
(f,h,f',h',\ldots)
\]

with covariant transport and its derived observables.

The bridge abstraction MUST be extracted only after both concrete sides are independently tested. It MUST NOT be designed first and then used to force both systems into a preconceived common interface.

---

## 12. Acceptance criteria for v1

Version 1 of the Sagnac reference layer is accepted only if all of the following hold:

- [ ] `SAG-S1` through `SAG-S9` all pass.
- [ ] PDE and segment-chain routes are algorithmically independent of the analytic oracle.
- [ ] Direction reversal is tested for positive and negative \(v\).
- [ ] Proper time and coordinate time are separate typed/named outputs.
- [ ] Phase is downstream of time transport and never used to generate the return-time result.
- [ ] At least one refinement study is stored as an artifact.
- [ ] At least one deliberate sign corruption fails.
- [ ] At least one deliberate frame/\(\gamma\) corruption fails.
- [ ] The \(|v|<c\) domain guard fails closed.
- [ ] No SSZ member file or SSZ metric function is imported by the reference layer.
- [ ] Existing `TRUE_FULL_CLOSURE_PASS` semantics remain unchanged.

---

## 13. Scientific interpretation rule

A passing reference block establishes only:

\[
\boxed{
\text{the implemented transport methodology reproduces the declared analytic Sagnac reference problem}
}
\]

It does **not** establish:

- that SSZ is empirically correct,
- that Sagnac uniquely favors SSZ,
- that the Sagnac correction-series zig-zag is a physical oscillator,
- that frame dragging or gravitational rotation has been modeled,
- that any later SSZ-to-Sagnac bridge is automatically valid.

Those are separate claims requiring separate evidence.

---

## 14. Recommended implementation order

1. `sagnac.py` analytic oracle + unit tests,
2. parity and inversion tests,
3. geometric correction-series convergence,
4. spatial segment chain,
5. independent PDE solver,
6. PDE ↔ chain ↔ oracle comparison,
7. proper-time conversion,
8. phase readout,
9. negative-control battery,
10. machine-readable evaluator and artifacts,
11. only then design the generic bridge contract.

This order keeps the reference problem independently falsifiable and minimizes accidental coupling to existing SSZ code.
