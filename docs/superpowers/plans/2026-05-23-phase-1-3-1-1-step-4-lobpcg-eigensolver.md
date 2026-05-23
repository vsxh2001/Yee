# Phase 1.3.1.1 step 4 — implementation plan

**Spec:** `docs/superpowers/specs/2026-05-23-phase-1-3-1-1-step-4-lobpcg-eigensolver-design.md`
**Base SHA:** `8595df8` (post Phase 4.fem.eig.3.5.5 merge)
**Lane:** `crates/yee-fem/src/solve.rs`,
`crates/yee-fem/src/lib.rs` (module doc only),
`crates/yee-fem/tests/lobpcg_smoke.rs`, `ROADMAP.md`,
`docs/src/decisions/0050-*.md`.
**Out of lane** (surface as findings, do NOT fix): `crates/yee-mom/**`
(its `eigensolver_wr90` test is a *consumer* — run it read-only as a
gate, do not edit), `crates/yee-validation/**`, `crates/yee-py/**`,
the complex arm `ComplexInverseIterEigen` (step-4.1 follow-on).

## Step ladder

### L1 — `LobpcgEigen` skeleton + block helpers

`crates/yee-fem/src/solve.rs`:

1. Add `pub struct LobpcgEigen { max_iter, tol, guard }` + `new` +
   `Default` (`max_iter = 1000`, `tol = 1e-8`, `guard = 2`) mirroring
   the `InverseIterEigen` doc style.
2. Add block helpers as free fns next to the existing ones:
   `block_m_orthonormalize` (M-orthonormalise an `n × b` block via
   modified Gram-Schmidt using `csr_matvec` + existing `dot`),
   `block_seed` (deterministic `n × b` seed, fixed-seed — extend
   `seed_vector`). Reuse `build_shifted`, `lu_solve`, `csr_matvec`,
   `dot`; do not duplicate.

**Verification:** `cargo check -p yee-fem` exit 0.

### L2 — LOBPCG outer loop + Rayleigh-Ritz

1. Implement `impl SparseEigen for LobpcgEigen` per spec §3.1
   Algorithm 4.1: residual block `R = KX − MXΛ`, preconditioned
   `W = T·R` (LU solve per column), search space `S = [X|W|P]`,
   dense `3b×3b` Rayleigh-Ritz (`SᵀKS`, `SᵀMS` via faer; Cholesky-
   reduce per spec risk (b) if faer lacks a direct generalized
   symmetric path), update `X`,`P`, per-column residual convergence on
   the leading `num_eigs`.
2. Soft-locking / `P`-column drop on near-singular `SᵀMS` (spec §3.2).
3. Return leading `num_eigs` pairs sorted ascending, M-orthonormal —
   match the `EigenpairList` postcondition contract.

**Verification:** `cargo test -p yee-fem --lib solve` exit 0.

### L3 — Unit tests (mirror InverseIterEigen + add cluster test)

`crates/yee-fem/src/solve.rs` `#[cfg(test)]`:

1. Duplicate the three existing pencil tests for `LobpcgEigen`:
   known 4×4 eigenvalues {0.5, 1.2, 3.4, 7.8} to `1e-8`;
   scaled-identity; block M-orthogonality `eᵀMe ≈ I`.
2. **New DoD-V5 degenerate-cluster test:** synthetic pencil with a
   known double root; assert both members returned, each residual
   `< tol`, mutual M-orthonormality to `1e-6`. (`InverseIterEigen`
   may be added as a contrast assertion only if it does not flake.)

**Verification:** `cargo test -p yee-fem --lib solve` exit 0;
`cargo test -p yee-fem --test lobpcg_smoke` exit 0 (extend the smoke to
exercise `LobpcgEigen`).

### L4 — Consumer gates (read-only) + ROADMAP + ADR

1. Run, do not edit: `cargo test -p yee-mom --test eigensolver_wr90`
   and `cargo test -p yee-validation --test fem_eig_001_rectangular_cavity`.
   Record LOBPCG vs InverseIterEigen iteration counts in the L4 commit
   body. If either gate consumes a *fixed* `InverseIterEigen` and
   cannot be pointed at `LobpcgEigen` without an out-of-lane edit,
   surface as a finding — the in-lane DoD is the `yee-fem` self-tests
   (L3); the consumer swap is a follow-on.
2. `crates/yee-fem/src/lib.rs` + `solve.rs` header: document LOBPCG
   availability, the complex-arm boundary, the arpack deferral.
3. `ROADMAP.md`: mark Phase 1.3.1.1 step 4 shipped with merge SHA.
4. ADR-0050 `docs/src/decisions/0050-phase-1-3-1-1-step-4-lobpcg-eigensolver.md`:
   in-tree pure-Rust block LOBPCG over arpack-rs (spec §5 rationale).

**Verification:** `mdbook build docs/` exit 0;
`grep -n "step 4" ROADMAP.md` shows the shipped marker.

## Lint floor (every commit)

```bash
cargo fmt --check --all
cargo clippy --workspace --all-targets -- -D warnings
```

## Escape-hatch

Blocked > 20 min on the Rayleigh-Ritz dense generalized eigensolve
(faer API uncertainty, spec risk (b)) → fall back to the Cholesky-
reduction path and document it; if *that* is also blocked, commit
L1+L3-skeleton with the cluster test `#[ignore]`'d and surface the faer
generalized-eigensolve gap as a finding. Do NOT pull in an external
eigensolver crate to unblock — that reverses the spec §5 decision and
must go back through an ADR.

## Out-of-scope (surface as findings, do not fix)

* `ComplexLobpcgEigen` (step-4.1).
* Swapping the `fem-eig-001` / `eigensolver_wr90` *default* solver to
  LOBPCG (consumer-lane change; this phase only adds the impl + proves
  parity).
* `arpack` optional feature.
