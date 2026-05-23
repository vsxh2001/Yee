# Phase 4.fem.eig.3.5.5 — retire fem-eig-006 via frequency retune off the cutoff edge

**Status:** Draft
**Owner:** TBD
**Phase:** 4.fem.eig.3.5.5 (ADR-0048 Option (a): retune
`FEM_EIG_006_F_HZ` so the multi-mode wave-port basis landed in
v3.5.4 carries real propagating modal content).
**Depends on:** Phase 4.fem.eig.3.5.4 (multi-mode `PortDefinition`
API + 3-mode fem-eig-006 driver; merge `d6611d2`).
**Blocks:** retirement of `#[ignore]` on
`fem_eig_006_magnitude_bounded`; closure of the fem-eig-006 line.

## 1. Goal

Retire `fem_eig_006_magnitude_bounded` (`|S_{11}| < 0.1`) by moving
the test frequency from **30 GHz** (the exact TE_{20} cutoff, where
the v3.5.4 multi-mode basis collapses to single-mode) to **40 GHz**,
where TE_{20} propagates with non-trivial `β` and the multi-mode
wave-port has real modal content to terminate.

Tolerance `< 0.1` is **not** weakened.

## 2. Background

### 2.1 v3.5.4 cutoff-degeneracy finding

ADR-0048 closed the v3.5.4 measurement: at 30 GHz the multi-mode
basis `[TE_{10}, TE_{20}, TE_{01}]` collapses to single-mode
because:

| mode    | cutoff `f_c`         | at 30 GHz   |
|---------|----------------------|-------------|
| TE_{10} | `c/(2B) = 15.0 GHz`  | propagating |
| TE_{20} | `c/B   = 30.0 GHz`   | **at cutoff, β = 0** |
| TE_{01} | `c/(2D) = 150.0 GHz` | evanescent  |

(`B = 10 mm` port-face broad wall, `D = 1 mm` narrow wall.)

`β = 0` ⇒ the TE_{20} per-mode stiffness block vanishes identically;
the basis carries no second propagating direction.

### 2.2 Why 40 GHz

At 40 GHz:

| mode    | cutoff `f_c`         | at 40 GHz   | `β` (rad/m) |
|---------|----------------------|-------------|-------------|
| TE_{10} | 15.0 GHz             | propagating | `sqrt((ω/c)² − (π/B)²) ≈ 776` |
| TE_{20} | 30.0 GHz             | **propagating** | `sqrt((ω/c)² − (2π/B)²) ≈ 554` |
| TE_{01} | 150.0 GHz            | evanescent  | 0 |

TE_{20} now carries real propagating content (`β ≈ 554 rad/m`,
33% above cutoff). The multi-mode wave-port terminates **both**
TE_{10} and TE_{20}; the residual reflection measures how well the
two-mode modal termination matches the field at the +x face.

### 2.3 Expected result + benchmark provenance

fem-eig-006 is a **synthetic stress fixture**, not an external
published benchmark (cf. mom-001 NEC-4). The defensible physics
claim at 40 GHz is the **matched-modal-termination identity**: a
wave-port whose modal basis spans the propagating modal content of
the port face produces near-zero reflection of those modes. The
gate `|S_{11}| < 0.1` therefore tests "does the v3.5.4 multi-mode
wave-port correctly match a two-propagating-mode port" — a
known-physics regression check, not an external-reference accuracy
gate. This is the same class of self-consistency check as the
fem-eig-004 / fem-eig-005 single-mode matched-port gates already in
the suite.

If `|S_{11}|(40 GHz) ≥ 0.1`, the multi-mode termination has a real
defect (mode orthogonality, normalisation, or β sign) — escape-hatch
to v3.5.6 with the measurement, do **not** weaken the tolerance.

## 3. Approach

### 3.1 Frequency constant

`crates/yee-validation/src/lib.rs`: change
`FEM_EIG_006_F_HZ` from `30.0e9` to `40.0e9`. Update the
doc-comment with the cutoff table from §2.2.

### 3.2 Driver: TE_{20} now a propagating mode

The v3.5.4 driver already constructs the 3-mode basis. At 40 GHz
the `fem_eig_006_beta_te20` closure returns `β ≈ 554` (non-zero)
automatically — no closure change needed, the `arg > 0` branch now
fires. Verify the TE_{20} `a_inc` stays `Complex64::ZERO`
(outgoing-only); the driving mode remains TE_{10} with
`a_inc = ONE`.

Optionally drop TE_{01} from the basis (still evanescent at 40 GHz,
contributes nothing) — keep it for forward-compatibility documenting
the 150 GHz cutoff, or drop it for a cleaner two-mode basis. **Spec
recommendation: keep all three** so the basis is frequency-agnostic.

### 3.3 Gate disposition

If `|S_{11}|(40 GHz) < 0.1`:
- Remove `#[ignore]` from `fem_eig_006_magnitude_bounded`.
- Rewrite the docstring as a "Phase 4.fem.eig.3.5.5 retire"
  record citing the 40 GHz measurement and the matched-modal-
  termination rationale.

Else: escape-hatch (keep `#[ignore]`, log 40 GHz measurement,
queue v3.5.6 multi-mode-port defect investigation).

## 4. Risks

(a) **Mesh resolution at 40 GHz.** `λ₀(40 GHz) = 7.5 mm`; the
cavity `(16, 3, 2)` mesh has `B/3 ≈ 3.3 mm` transverse cells. At
40 GHz that is ~2.3 cells/λ transverse — coarse. If `|S_{11}|`
is dominated by discretisation error rather than modal mismatch,
the gate may not retire even with correct multi-mode physics.
Mitigation: bump `FEM_EIG_006_NY` / `FEM_EIG_006_NZ` if needed
(measure first, refine only if discretisation-limited).

(b) **TE_{20} field-pattern orientation.** Verify the
`fem_eig_006_modal_e_t_te20` closure's world-frame mapping matches
the actual +x-face cross-section orientation; a transposed pattern
projects onto the wrong subspace and inflates `|S_{11}|`.

## 5. Definition of done

DoD-1. `FEM_EIG_006_F_HZ = 40.0e9`; doc-comment carries the §2.2
cutoff table.

DoD-2. `cargo test --release -p yee-validation --test
fem_eig_006_high_aspect_pml -- --include-ignored --nocapture`
prints the 40 GHz `|S_{11}|`.

DoD-3. Gate disposition per §3.3 (retire if `< 0.1`, else
escape-hatch with logged measurement).

DoD-4. Tutorial `07-fem-open-cavity.md` + `ROADMAP.md` carry the
v3.5.5 subsection / line. If retired, mark the fem-eig-006 line
**closed**.

DoD-5. Lint floor clean (`cargo fmt --check --all` + `cargo clippy
--workspace --all-targets -- -D warnings`).

## 6. References

* ADR-0048 Option (a)
  `docs/src/decisions/0048-phase-4-fem-eig-3-5-5-disposition.md`.
* ADR-0047 — multi-mode wave-port API.
* Pozar, *Microwave Engineering*, 4th ed., §3.3 TE_{mn} cutoff +
  field patterns.
* `crates/yee-validation/src/lib.rs` — fem-eig-006 driver +
  `FEM_EIG_006_F_HZ`.
