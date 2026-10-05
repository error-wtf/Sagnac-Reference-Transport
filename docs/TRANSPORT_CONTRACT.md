# Future abstract transport contract (documented, NOT implemented)

## The chain

    S (signal law)  ->  A (architecture)  ->  L[A] (transport operator)
        ->  Psi (propagated field)  ->  O (observables)

## Sagnac instantiation (this repository)

    S_Sag   = (L, c, v)
    A_Sag   = ring [0,L), emitter/receiver co-located at x=0,
              receiver moving at v (periodic)
    L_±     = (d/dt ± c d/dx)          two counter-propagating branches
    Psi     = advected pulse fields u_±(x,t)
    O_Sag   = (t_+, t_-, dt_axle, dt_det_proper, dphi)

## Later SSZ instantiation (out of scope here)

    S_SSZ   = (f(r), h(r), f', h', ...)
    A_SSZ   = spherically symmetric segment geometry
    L_SSZ   = nabla[g(f,h)]              Christoffel transport
    O_SSZ   = (geodesics, S_r phase, redshift, trapping, ...)

## Rule

Only AFTER both sides are independently validated (Sagnac via this repo,
SSZ via the registered TRUE FULL CLOSURE) may a shared `TransportProblem`
abstraction be extracted.  Designing the abstraction first risks shaping
it so that SSZ fits automatically — the circularity this program exists
to prevent.

The bridge itself will live in a SEPARATE package importing both
`sagnac_reference` (this repo) and the SSZ closure internals read-only.
