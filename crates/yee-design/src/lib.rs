//! # yee-design
//!
//! Natural-language design surface for the Yee electromagnetic simulation
//! studio (Phase 3.nl.0 walking skeleton).
//!
//! At base SHA this crate carries only the type surface — the typed
//! [`DesignIntent`] produced by Stage 1 of the spec §5 pipeline and the
//! canonical substrate library it references. Stages 2–5 (initial-estimate
//! synthesis, deterministic TOML emit, offline-mode parser, CLI wiring) land
//! in subsequent tracks (R2 – R5) of the
//! `2026-05-18-phase-3-nl-0-design-surface.md` plan.
//!
//! The type surface is `serde`-round-trippable by contract: the integration
//! test in `tests/intent.rs` round-trips 100 randomised samples through
//! `serde_json::{to_string, from_str}` byte-identically. Downstream stages
//! (e.g. the spec §8 `<out>.intent.json` artefact, the Python sidecar's
//! schema-validated tool-use response in `yee-py`) depend on this property.
//!
//! See:
//! - `docs/superpowers/specs/2026-05-18-phase-3-nl-0-design-surface-design.md`
//!   §6 for the type surface and §7 for the JSON schema this crate is
//!   ultimately constrained to.
//! - `docs/superpowers/plans/2026-05-18-phase-3-nl-0-design-surface.md`
//!   Step R1 for this crate's scoping.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod emit;
pub mod estimate;
pub mod intent;

pub use emit::{ProjectFile, emit};
pub use estimate::{Error, InitialEstimate, ResolvedSubstrate};
pub use intent::{
    DesignIntent, GeometryFamily, NamedSubstrate, Provenance, Substrate, SubstrateLibrary,
    SubstrateOverride, SubstrateRecord, substrate_library,
};
