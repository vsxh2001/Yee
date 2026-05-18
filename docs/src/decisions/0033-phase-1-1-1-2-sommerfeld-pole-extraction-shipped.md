# ADR-0033: Phase 1.1.1.2 Sommerfeld pole extraction — implementation shipped

## Status

Accepted — 2026-05-18.

## Context

ADR-0025 locked the spec and plan for Phase 1.1.1.2 Sommerfeld
surface-wave pole extraction; implementation was deferred to a
follow-up track. Track DDDDD landed the spec
(`2026-05-17-phase-1-1-1-2-sommerfeld-pole-extraction-design.md`)
and plan (same date, plans/); the implementation was scheduled
as the next non-trivial `yee-mom` sub-project.

Track JJJJJ (merge `a22d622`) ships that implementation. The
six steps from the ADR-0025 plan are all green: Newton–Raphson
root finder for `D_TM(k_ρ) = 0` and `D_TE(k_ρ) = 0` in the
complex `k_ρ` plane (commit f98bf14) with the closed-form
derivative; in-house Bessel `J_0` and Hankel `H_0^{(2)}`
evaluators; pole-subtracted GPOF on the smooth residual
`G̃_residual`; space-domain reconstruction
`Σ_n c_n · exp(−jk_0·r_n)/r_n + (j/4) · R_p · H_0^{(2)}(k_p · ρ)`;
a new `crates/yee-mom/src/sommerfeld.rs` module; and the
`GreensSpec::MicrostripSommerfeld` constructor (`n_surface_wave_poles = 0`
remains the bit-for-bit ADR-0020 tripwire).

Validation: 60/60 `yee-mom` tests pass, including the new
`sommerfeld_synthetic_pole_recovery` test (Newton converges to a
known synthetic pole within tolerance) and the
`sommerfeld_fr4_1ghz_finds_tm0_pole` smoke test (FR-4 at 1 GHz,
the TM_0 mode lands where Pozar §3.7 says it should).

## Decision

The Phase 1.1.1.2 implementation locks four load-bearing choices
the ADR-0025 spec was silent on:

1. **Sommerfeld machinery in a new `crates/yee-mom/src/sommerfeld.rs`
   module,** not extensions to `multilayer.rs`. The split makes
   pole-finding, residue-extraction, and Hankel-reconstruction
   testable in isolation; `multilayer.rs` consumes
   `sommerfeld::extract_poles` and `sommerfeld::reconstruct`
   through a narrow interface.
2. **Bessel `J_0` and Hankel `H_0^{(2)}` are implemented
   in-house,** not pulled from `special-functions` or a similar
   external crate. Small-argument power series + large-argument
   asymptotic; ~150 lines total, unit-tested against tabulated
   values, one fewer audit-surface dependency.
3. **Newton–Raphson with the analytic `D'(k_ρ)` derivative; no
   numerical-differentiation fallback.** The dispersion equation's
   derivative is closed form, and the ADR-0025 plan's analytic
   seed values land Newton in the right basin reliably for
   FR-4-1GHz-class problems. Müller's-method fallback from the
   spec is **not** wired here; a bad-basin case will surface as
   a `DegeneratePole` error rather than silently degrade.
4. **Pole-subtracted GPOF, not pole-aware GPOF.** Subtract the
   analytic pole contribution from the sampled spectral kernel
   *before* GPOF runs on the smooth residual; add
   `(j/4) · R_p · H_0^{(2)}(k_p · ρ)` back analytically in the
   space domain. This avoids GPOF trying to fit the slow
   `H_0^{(2)}` radial decay with complex exponentials — which it
   cannot do — and keeps the GPOF fit budget free for the truly
   smooth part of the kernel.

## Consequences

**What this unblocks.** `mom-002` (microstrip Z₀) is now eligible
for the tolerance retest against the Hammerstad–Jensen
`[35, 75] Ω` target; the loose tolerance noted in CLAUDE.md §4
and the ROADMAP "Outstanding validation gates" was gated on
Phase 1.1.1.2 landing — that gate has now cleared. `mom-003`
(2.4 GHz patch) inherits the same fix once the patch test is
re-run through `GreensSpec::MicrostripSommerfeld`.

**What is deferred.** Multi-pole extraction for higher-order
modes (TM_1, TE_2, …) and lossy-substrate poles (Sommerfeld
contour deformation around branch cuts) are queued for Phase
1.1.1.3. The current implementation handles the dominant TM_0
and TE_1 poles only; the ADR-0025 `n_surface_wave_poles`
parameter is plumbed end-to-end so 1.1.1.3 is an additive
change, not an API break.

## References

- ADR-0025 — Phase 1.1.1.2 Sommerfeld pole-extraction spec
  (the spec this ADR's implementation realises).
- `docs/superpowers/specs/2026-05-17-phase-1-1-1-2-sommerfeld-pole-extraction-design.md`
  and the matching plan under `docs/superpowers/plans/`.
- `docs/src/theory/multilayer-greens.md` — DCIM + surface-wave
  pole theory chapter.
- Commits f98bf14 (Newton–Raphson pole search), 6873302
  (`MicrostripSommerfeld` wiring); Track JJJJJ merge SHA
  `a22d622`.
- `crates/yee-mom/src/sommerfeld.rs` — pole-finding,
  residue extraction, Hankel reconstruction.
- ADR-0020 — multi-image DCIM via GPOF; the smooth-residual
  fit that this implementation's pole subtraction feeds into.
- CLAUDE.md §4 / §10 — `mom-002` / `mom-003` loose-tolerance
  caveat; now retestable post this ADR.
