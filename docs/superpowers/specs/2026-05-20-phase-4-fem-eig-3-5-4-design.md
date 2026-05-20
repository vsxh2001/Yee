# Phase 4.fem.eig.3.5.4 — multi-mode wave-port retirement of fem-eig-006

**Status:** Draft
**Owner:** TBD
**Phase:** 4.fem.eig.3.5.4 (retire the `fem_eig_006_magnitude_bounded`
strict gate `|S_{11}(30 GHz)| < 0.1` left in `#[ignore]` after the
VVVVVVVVV Phase 4.fem.eig.3.5.3 W1 single-mode wave-port attempt).
**Depends on:** Phase 4.fem.eig.3.5.3 (W1 single-mode `WavePort(1)`
on +x face landed; measured `|S_{11}|(30 GHz) = 0.925644 (-0.67 dB)`,
matching v3.5.2 CFS-PML 0.926 within numerical noise; merge SHA
`84a92e7`). ADR-0046 §Decision (5) blueprint.
**Blocks:** retirement of `#[ignore]` on
`fem_eig_006_magnitude_bounded`; closure of the fem-eig-006 line in
the Phase 4.fem.eig.3.5 chain.

## 1. Goal

Retire the **`fem_eig_006_magnitude_bounded`** strict-magnitude gate
`|S_{11}(30 GHz)| < 0.1` by extending the +x `PortDefinition` from a
TE_{10}-only modal projection to a **multi-mode termination
spanning {TE_{10}, TE_{20}, TE_{01}}**, per Jin §10.6 multi-mode
wave-port and per ADR-0046 §Decision (5).

Tolerance `< 0.1` is **not** weakened.

## 2. Background

### 2.1 v3.5.3 W1 measurement

Track VVVVVVVVV (commits `4b3316b` → `c89985d` → `223cddd`; merge
`84a92e7`) swapped the +x face of the high-aspect 100 : 10 : 1
cavity driver from a 14-cell CFS-PML shell to a single
`FaceKind::WavePort(1)` carrying a TE_{10} modal projection. The
native (16, 3, 2) Kuhn-6 mesh runs at **576 tets** (versus the
v3.5.2 ~580 extended PML mesh — net wash on assembly cost).

Measured at 30 GHz:

| Configuration                 | `|S_{11}|`  | dB     |
|-------------------------------|-------------|--------|
| v3.5.2 best H4 row CFS-PML    | 0.926       | -0.67  |
| v3.5.3 W1 TE_{10} wave-port   | 0.925644    | -0.67  |

The two terminations agree to four decimals — the v3.5.2 PML floor
and the v3.5.3 single-mode wave-port floor are **the same physics
expressed in different boundary operators**, both projecting the
field onto a single tangential basis vector and discarding the
orthogonal modal content.

### 2.2 TE_{20} cutoff at 30 GHz

For a rectangular waveguide of broad-wall `a` and narrow-wall `b`,
the TE_{mn} cutoff frequency is

```
f_c,mn = (c / 2) · sqrt((m/a)² + (n/b)²)
```

With `a = 100 mm` (+x face full broad-wall span) and `b = 10 mm`,
**TE_{20}** has `f_c = c·(2/0.100)/2 = c·10 m⁻¹ = 3.00 GHz`. The
**TE_{01}** mode has `f_c = c·(1/0.010)/2 = c·50 m⁻¹ = 15.00 GHz`.
Both are **propagating at 30 GHz** on the +x face cross-section
along with TE_{10}.

A single-mode wave-port `PortDefinition` whose `modal_e_t` returns
only the TE_{10} shape orthogonal-projects the field at the port
onto a one-dimensional modal subspace. Components in the TE_{20}
and TE_{01} directions of the modal Hilbert space are aliased into
reflection at the port — the floor measured in v3.5.3 (`|S_{11}|
= 0.926`) is the **modal content this projection cannot represent**.

### 2.3 ADR-0046 §Decision (5)

The ADR identifies multi-mode extension as the binding follow-on:

> Add `TE_{20}` and `TE_{01}` mode shapes to the +x `PortDefinition`,
> projecting the field onto a 3-D modal subspace at the port face.
> Expect `|S_{11}|` to drop by ~20 dB per added mode in the
> low-modal-content regime; sufficient at 30 GHz to clear the
> `< 0.1` (`-20 dB`) tolerance.

## 3. Approach

### 3.1 API extension: `PortDefinition` → `Vec<PortMode>`

The existing single-mode struct (`crates/yee-fem/src/open_boundary.rs:573`):

```rust
pub struct PortDefinition {
    pub beta_mode: Box<dyn Fn(f64) -> f64 + Send + Sync>,
    pub modal_e_t: Box<dyn Fn(Vector3<f64>) -> Vector3<f64> + Send + Sync>,
}
```

becomes a **modal-basis container**:

```rust
pub struct PortMode {
    pub beta_mode: Box<dyn Fn(f64) -> f64 + Send + Sync>,
    pub modal_e_t: Box<dyn Fn(Vector3<f64>) -> Vector3<f64> + Send + Sync>,
    /// Incident amplitude scaling for this mode. The TE_{10}
    /// driving mode carries the full incident amplitude; higher-
    /// order modes carry zero (they exist only as projection
    /// directions for outgoing scattering).
    pub a_inc: Complex64,
}

pub struct PortDefinition {
    pub modes: Vec<PortMode>,
}
```

`PortDefinition` becomes a thin wrapper over `Vec<PortMode>`. The
v3.5.3 single-mode call site collapses to
`PortDefinition { modes: vec![PortMode { beta_mode, modal_e_t, a_inc: Complex64::ONE }] }`.

### 3.2 Assembly path

`assemble_port_face_block` and `assemble_port_modal_rhs` (in
`crates/yee-fem/src/element/port.rs`) currently consume the single
closure pair. The multi-mode form sums over the modal basis:

```
K_port^p = Σ_m  K_port^{p, m}    (stiffness block per mode)
b_port^p = Σ_m  a_inc_m · b_port^{p, m}   (RHS per mode)
```

On the v3.5.3 side, the post-solve `S_{p,p}` extraction in
`OpenBoundarySolver::solve_at_frequency` reads a single inner
product `⟨E_h, e_t^{port}⟩`. The multi-mode S-parameter extraction
becomes a **modal-decomposition step**: for each mode `m`,
`S_{p,m}(ω) = ⟨E_h, e_t^{p,m}⟩ / ⟨e_t^{p,m}, e_t^{p,m}⟩`, and the
**driving-mode reflection** `S_{p,p}` is `S_{p, m₀}` where `m₀` is
the mode carrying `a_inc = 1` (the TE_{10} driving mode).

For v3.5.4 the **strict gate uses `|S_{p,m₀}|`** — same scalar the
v3.5.3 gate measured, now correctly projected against a
multi-dimensional modal subspace rather than a one-D collapse.

### 3.3 Caller-side construction (fem-eig-006 driver)

`run_fem_eig_006_high_aspect_pml` in
`crates/yee-validation/src/lib.rs` already constructs a single
`PortDefinition` for the +x face. v3.5.4 changes it to:

```rust
let port_modes = vec![
    PortMode {
        beta_mode: Box::new(fem_eig_006_beta_te10),
        modal_e_t: Box::new(fem_eig_006_modal_e_t_te10),
        a_inc: Complex64::ONE,        // driving mode
    },
    PortMode {
        beta_mode: Box::new(fem_eig_006_beta_te20),
        modal_e_t: Box::new(fem_eig_006_modal_e_t_te20),
        a_inc: Complex64::ZERO,       // outgoing-only
    },
    PortMode {
        beta_mode: Box::new(fem_eig_006_beta_te01),
        modal_e_t: Box::new(fem_eig_006_modal_e_t_te01),
        a_inc: Complex64::ZERO,       // outgoing-only
    },
];
```

The three modal `e_t(x, y, z)` closures follow the standard
TE_{mn} field-pattern recipe at the port face — sinusoidal in the
broad-wall coordinate, with peak orientation determined by `(m, n)`.

## 4. Risks

(a) **Cross-mode coupling at the port** — the assembly path
assumed an orthogonal modal basis. Verify the three TE shapes
remain numerically orthogonal on the Kuhn-6 face mesh; if not, add
a Gram-Schmidt normalisation pre-projection.

(b) **API churn** — `PortDefinition` is a public type re-exported
through `yee.fem.solve_open_cavity`. The Python kwarg shape
(`port_faces=[{"axis": "x", ..., "modal_e_t": (0.0, 1.0, 0.0)}]`)
needs a multi-mode list shape. Defer Python wiring to a follow-on
(v3.5.4.1) if shape design requires user input.

(c) **Insufficient improvement** — if `|S_{p,m₀}|` lands in
`[0.05, 0.10)`, the gate retires but with low margin. If it lands
`≥ 0.1` (still over threshold), escape-hatch: keep `#[ignore]`, log
the multi-mode measurement, queue (a) higher-order mode
augmentation (TE_{11}, TE_{30}) or (b) absorbing-mode termination
extension to v3.5.5.

## 5. Definition of done

DoD-1. `PortDefinition` carries `Vec<PortMode>`; existing
fem-eig-004 / fem-eig-005 / fem-eig-006 drivers compile against the
new shape (single-mode call sites use the `vec![PortMode { ... }]`
collapse).

DoD-2. `cargo test --workspace` green on the default path
(`fem_eig_006_magnitude_bounded` no longer `#[ignore]`'d if
DoD-3 succeeds).

DoD-3. `fem_eig_006_magnitude_bounded` with the multi-mode +x
port reports `|S_{11}|(30 GHz) < 0.1`. If not, DoD-3 invokes
escape-hatch: `#[ignore]` stays, measurement is logged in the
docstring with the same level of detail as the v3.5.3 record.

DoD-4. Tutorial `07-fem-open-cavity.md` and `ROADMAP.md` carry the
multi-mode wave-port subsection and the v3.5.4 ROADMAP line.

DoD-5. Lint floor (`cargo fmt --check --all` + `cargo clippy
--workspace --all-targets -- -D warnings`) clean on the merge
commit.

## 6. References

* ADR-0046, `docs/src/decisions/0046-phase-4-fem-eig-3-5-3-fem-eig-006-retire.md` §Decision (5).
* Jin, *FEM in EM*, 3rd ed., Chapter 10.6 "Wave-port termination" — modal-basis derivation.
* Pozar, *Microwave Engineering*, 4th ed., §3.3 TE_{mn} cutoff frequencies.
* `crates/yee-fem/src/open_boundary.rs:573` — current `PortDefinition`.
* `crates/yee-fem/src/element/port.rs` — current assembly path.
* `crates/yee-validation/src/lib.rs` — v3.5.3 fem-eig-006 driver.
* Phase 4.fem.eig.3.5.3 spec `docs/superpowers/specs/2026-05-20-phase-4-fem-eig-3-5-3-design.md`.
