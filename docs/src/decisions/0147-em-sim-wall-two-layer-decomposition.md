# ADR-0147: Full-board EM-sim wall — two-layer decomposition + the path through

**Status:** Accepted (research finding; productionization is a maintainer-funded decision)
**Date:** 2026-05-31
**Related:** ADR-0133 (the cavity-vs-CW wall), ADR-0134 (the F2.3 gate re-scope),
ADR-0124/0125 (the single-cell aperture-port reactance floor), ADR-0108 (CPML unstable
into substrate), ADR-0064 (planar-MoM microstrip-port ill-posedness), ADR-0070/fem-eig-006
(FEM wave-port modal saturation), [[lumped-lc-and-studio-redesign]],
[[project-filter-design-final-goal]]

---

## Context

The studio designs filters end-to-end (synthesis → board → Gerber/KiCad) but its
**Verify** stage is circuit-level only: full-board EM verification of a high-Q microstrip
filter was deferred as a wall (ADR-0133). The maintainer commissioned a **research team**
(team `em-wall`, 2026-05-31) — one agent per candidate direction + a skeptic **oracle**
validating every claim against the analytic ground truth (`yee_filter::ladder_s21` of a
3-pole Chebyshev 0.5 dB BP, f0 = 2 GHz, FBW = 10%: passband 0 dB, rejection 1.6 GHz
−41.8 dB / 2.4 GHz −36.3 dB, −3 dB edges 1.887/2.120 GHz; **geometrically asymmetric** —
low-side rejection deeper — a key discriminator). Directions: transient pole-extraction
(`prony`), substrate-stable CFS-PML (`cfspml`), in-FDTD TRL (`trl`), FEM driven-sweep
(`femmor`).

## Decision (the finding)

**The wall is not one problem — it decomposes into two independent layers, and no
single-spike direction broke both on the real 3-D microstrip.** This decomposition is the
team's central, oracle-validated output:

- **Layer 1 — cavity vs CW steady state (ADR-0133):** a high-Q microstrip filter must ring
  up (CW), but the only *stable* microstrip box is a PEC cavity whose own modes dominate
  (lossless thru reads ~1.3× over-unity at a box mode), and the cavity-killing absorber
  (CPML) was believed unstable into the substrate. **SOLVED two independent ways** — `prony`
  (extract-in-post: matrix-pencil pole separation) **and** `cfspml` (absorb-in-sim: a
  √εr-matched substrate CPML — the "CPML unstable into substrate" sub-claim is itself
  *superseded*, the binding knob was σ_max medium-calibration, not an intrinsic instability).
  `trl` proved de-embed *cannot* (the cavity is a non-cascadable parallel bypass), which is
  why removal — not de-embedding — is the route.
- **Layer 2 — the FDTD filter *realization* (ADR-0125):** even with the cavity removed,
  the embedded-resonator filter-DUT in FDTD does not reproduce the analytic filter — the
  **single-cell aperture port floors at a frequency-flat shunt capacitance** (ADR-0125;
  Q≈10 tanks give ~1.8 dB stopband vs 20 dB target). **OPEN.** This is the shared blocker
  that both FDTD de-embedding routes (prony, trl) hit; it is **NOT** the ADR-0064
  E_z-orthogonality wall (that is planar-MoM-specific — FDTD carries full 3-D E on the Yee
  grid and does not inherit it).

  **Layer 2 splits into two distinct sub-blockers** (oracle, evidence-backed):
  - **L2a — realization *shape*:** the FDTD filter's passband is too narrow / wrong-shape
    (~10 dB off at edges, FBW-independent) because the resonators carry uniform loaded-Q
    instead of **g-value-scaled J-inverters** (Pozar §8.8). **FIXABLE, ~hours, not
    fundamental** — prony's `ladder` control proves that embedding proper g-scaled coupled
    L,C through the same machinery reproduces the analytic curve *including* the asymmetry
    (analytic 5.51 dB low-side-deeper @1.6 vs 2.4 GHz; ladder 5.9 dB, within 0.4 dB). The
    3-D analog = proper distributed coupled-line synthesis.
  - **L2b — aperture-port co-location:** clustered single-cell aperture ports degenerate to
    one lumped short → flat 0.04, *no* band-pass at all (ADR-0125, "fundamental to the
    single-cell formulation"). The genuinely hard one — needs the **multi-cell aperture
    port** (the deferred F2.3 brick, ADR-0124/0125 "genuinely required", sketched but never
    built; unproven, multi-increment).

  The two are independent; a 3-D break must close **both**. Neither is the ADR-0064 wall.

## Ranked verdict (oracle-validated: source read, numbers reproduced, spikes re-run)

1. **FEM driven-sweep (`femmor`) — PASS, soundest path.** Sidesteps **both** layers
   (frequency-domain: one complex sparse LU per frequency with wave-port absorbing BCs —
   no cavity, no CW ring-up, no CPML instability). Independently reproduced on 3-D
   waveguide: air thru |S21| = −0.045 dB / |S11| = −53 dB / recip 2e-15; **dielectric**
   (εr 4.4) thru port phase-velocity 0.28% + dispersion β within 1.3% (εr back-out
   4.41–4.51) — the dielectric driven path + port fidelity are now *proven* (never tested
   before). Gap: no on-target microstrip-*filter* S21 — needs (a) interior-PEC trace
   support (~½ day; `with_extra_pec_edges`, verified missing from the public API),
   (b) quasi-TEM modal closure β(ω)+e_t(x) from a 2-D cross-section solve (~3–5 days, the
   real unknown), (c) full-board LU scaling (~3–6e5 tets — faer direct LU infeasible;
   needs **iterative/AMG** or Gmsh local refinement). **~2–3 weeks + a solver-scaling
   risk.** The only path with no unsolved-physics blocker. Gate before any filter geometry:
   a straight-microstrip-thru with ε_eff within 5% of 3.33 (HJ).
2. **Transient pole-extraction (`prony`) — PARTIAL PASS: Layer 1 SOLVED, Layer 2 open.**
   Matrix-pencil pole-fit of a broadband-pulse transient, separating filter poles from
   the cavity pole by Q + bare-thru match, then subtracting. On an **independently realized**
   coupled-resonator-on-a-line (tanks sized from target loaded-Q per Pozar §8.8 — *not*
   the synthesized L,C; graded vs an FDTD-clean run, **non-circular**, verified in source):
   **Q-insensitive** — 0.138 dB RMS flat from cavity-Q = 1000 to **∞ (lossless)**, the exact
   regime that makes CW de-embed read over-unity, 90× better than naive DFT-ratio
   (12.6 dB); **robust** to 5% FBW (12× Q gap holds). This genuinely retires Layer 1.
   **FAIL end-to-end:** graded vs the analytic reference the recon is 10–12 dB too low at
   the band edges and fails the geometric-asymmetry gate — because the embedded FDTD filter
   is the *wrong* filter (the Layer-2 realization gap). Still a 1-D telegrapher analog with
   an *injected* box mode (not the real 3-D PEC eigenspectrum). NB: prony's second "ladder"
   driver (synthesized L,C run as a circuit ODE) reproduces `ladder_s21` **by construction**
   — a tautology; it is **not** counted as evidence (only the independent-realization
   `spike` bin is).
3. **Substrate CFS-PML (`cfspml`) — PARTIAL PASS, UPGRADED: cavity-removal *demonstrated*
   for a line.** Built + verified a **stable √εr-matched substrate-CPML termination**:
   the library default σ_max is vacuum-η₀-calibrated (→ only −9.5 dB into εr=4.4); scaling
   σ_max by **√εr (×2.1)** + a 24-cell layer → **−35.2 dB**, stable (energy env_ratio = 1.00)
   across dx = 0.4→0.1 mm. **Decisively (oracle-validated):** a lossless straight-microstrip
   THRU |S21| in the matched box reads **flat 1.37 dB ripple over 1.6–2.4 GHz vs 7.96 dB in
   the PEC box** — ~6× flatter, **no over-unity, no box mode**; the oracle confirmed the
   +0.7 dB low-freq lift is *dispersion* (monotone, not a peak at the 2.0 GHz mode), not
   residual cavity. So a **stable + cavity-free + matched microstrip box EXISTS** — this
   **contradicts ADR-0133's premise** (high-Q-CW ⊥ stable-cavity-free-box) for the *box
   half*, and the fix is small (√εr-σ + ~24 PML cells). STANDING LIMIT: demonstrated for a
   **bare line / single tank**, not the filter — the filter-in-matched-box |S21| was the
   identified decisive next test but was **not run/validated this session** (lost to a
   teammate-shutdown race), so no end-to-end break is *claimed* here.
4. **In-FDTD TRL (`trl`) — FAIL as a break, but the BEST DIAGNOSTIC.** TRL algebra correct
   + non-circular (verified: `calibrate()` never sees `ladder_s21`); cavity over-unity
   reproduced (|S21|_thru = 1.14–1.29 at the box mode). FAIL end-to-end: a *known* lossless
   line-DUT de-embeds to |S21|≈0.44 (should be 1.0) with **healthy conditioning ~0.9** — so
   not an ill-conditioning failure. **Mechanism proof (decisive, pure-math):** the cavity
   couples port1↔port2 as an **additive parallel bypass through the box volume**, around the
   DUT reference planes; TRL removes only **cascadable (series)** fixture error boxes
   (R = X·D·Y has no term for a parallel leak), so it **structurally cannot** remove the
   bypass — a **topological** limit, not a de-embed-technique gap. Reproduced exactly by
   injecting a parallel leak b (|S21|_rec 0.998→0.57 as b 0→0.5, cond healthy throughout —
   cond is a *trap* here). This is what **exposed** that Layer 2 is the shared blocker and
   that every cascade de-embed (DUT/thru, 2-point standing-wave, TRL) cannot beat a cavity.

## Consequences / the path through

**The through-line (trl + cfspml together):** trl *proved* de-embed is topologically dead
(the cavity is a non-cascadable parallel bypass) → the only viable cavity route is to
**remove** the cavity, not de-embed it; and there are now **two demonstrated cavity-kills** —
`cfspml` **absorb-in-sim** (matched √εr-σ box, THRU flat, oracle-validated) and `prony`
**extract-in-post** (matrix-pencil, Q-insensitive). Both still need the Layer-2 filter
realization to read the analytic curve end-to-end.

**Cheapest decisive near-term test (newly enabled, ~1 FDTD run — recommended first):** drop
`synthesize_lumped`'s exact L,C tanks (correct components → **sidesteps L2a** by
construction) into `cfspml`'s matched √εr-σ box via `yee-voxel::simulate_lumped_board`, and
read |S21| vs the locked `ladder_s21` *including the geometric-asymmetry gate* (1.6 GHz
deeper than 2.4 GHz). This isolates the **last unknown**: does the production `LumpedRlcPort`
realize the band-pass, or does it floor like trl's aperture-coupled DUT (**L2b**)? GO →
end-to-end break on the lumped path; NO-GO → L2b confirmed as the sole remaining blocker
even with a perfect box + perfect components. **Requires a skeptic validator** (the original
`oracle`'s gate; the locked reference + grading tooling are preserved in
`crates/yee-filter/examples/oracle_reference.rs` / `oracle_grade.rs`).

**Larger fallback (if L2b floors):** the **multi-cell aperture port** (the deferred,
unproven, multi-increment F2.3 brick) + a g-scaled *distributed* coupled-line realization
(the 3-D analog of L2a) + a cavity-kill — graded vs the locked curve. Multi-week; L2b is
the long pole.

**Long-term (the blocker-free path):** the FEM driven-sweep (`femmor`'s scoping above) —
the only direction with no unsolved physics; ~2–3 weeks + the solver-scaling decision
(recommend an **iterative/AMG** solver over faer direct LU; gate on the straight-microstrip
ε_eff<5% milestone before any filter geometry).

**Bankable now:** `prony`'s Layer-1 matrix-pencil cavity-separation is a real, validated,
reusable advance (it defeats the cavity-Q→∞ regime that no CW de-embed survives) — worth
productionizing as a reusable extraction tool in `yee-fdtd` independent of the Layer-2 fix.
`cfspml`'s stable √εr substrate line-termination is a real component for that work.

**Sequencing:** run the cheap matched-box + correct-L,C test *first* (≈1 FDTD run + a fresh
validator) — it either yields an end-to-end break on the lumped path or pins L2b as the sole
remaining blocker. Only if L2b floors does the **multi-week, maintainer-funded fork** open
(multi-cell-port + distributed realization, or the FEM driven-sweep). The studio Verify
stage stays honestly circuit-level until an EM path is validated. **No EM result was merged;
nothing was faked; the oracle rejected every match-by-construction** (the ODE-ladder
tautology, the wrong-test FEM number). The cfspml THRU cavity-removal upgrade recorded above
*was* oracle-validated; the filter-in-box |S21| was **not** run/validated this session and
is **not** claimed as a break.

---

## References
- Team `em-wall` (2026-05-31): `prony` (transient pole-extraction), `cfspml` (substrate
  CFS-PML), `trl` (in-FDTD TRL), `femmor` (FEM driven-sweep), `oracle` (validator).
- `crates/yee-voxel/src/lumped_sim.rs` (`simulate_lumped_board`); `crates/yee-fdtd`
  (CPML); `crates/yee-fem` (driven sweep + wave ports); `crates/yee-filter`
  (`ladder_s21`). ADR-0133/0134/0124/0125/0108/0064/0070.
