# yee-core

> Foundation crate: shared types, traits, errors, units, constants.

`yee-core` is the small, stable bottom of the Yee dependency graph. Every other Yee crate depends on it; it depends on no Yee crate. It has **no CUDA, no GUI, no I/O** — its compile must stay fast and its public surface must stay narrow.

## Scope

### Phase 0 (in progress)
- `units` — physical constants (c₀, ε₀, μ₀, η₀) and SI helpers
- `FreqRange` — linear frequency sweeps
- `Error` / `Result` — crate-wide error type with `thiserror`
- `Solver` trait — skeleton implemented by `yee-mom`, `yee-fdtd`

### Phase 1
- Material model traits (lossless, lossy dielectric, frequency-dispersive)
- Port abstractions (lumped, wave, modal) — concrete impls in solver crates
- S-parameter and field-data containers (`SParameters`, `FarFieldPattern`)
- Convergence and error-estimator traits

### Phase 2+
- Time-domain analog of `FreqRange` (`TimeSpan`, `Dt`)
- Dispersive-material trait family (Drude / Lorentz / Debye)

## Design rules

- **No `unsafe`** — `#![forbid(unsafe_code)]` is enforced.
- **No heavy deps.** Only `nalgebra`, `num-complex`, `thiserror`, `tracing`.
- **Documented public surface.** `missing_docs` warns; CI promotes to deny.
- **Stable ABI.** Breaking changes here cascade through the workspace; review carefully.

## Examples

```rust
use yee_core::{FreqRange, units::C0};

let band = FreqRange { start_hz: 2.0e9, stop_hz: 3.0e9, n_points: 201 };
let lambda_at_start = C0 / band.start_hz;
```

## Validation

See [`validation/README.md`](validation/README.md). Unit-level only — constants vs. CODATA, range arithmetic, error variant round-trips.

## Roadmap

See [`ROADMAP.md`](ROADMAP.md).
