//! `yee-studio` binary entry — a no-GUI stub.
//!
//! The eframe/egui desktop+web view of the Yee Filter Studio was **retired** in
//! App.D.2 (ADR-0130): the pure-Rust **Dioxus** `yee-studio-web` is now the
//! studio (the goal's polished-UI component). This crate retains only the
//! egui-free, WASM-safe [`StudioState`](yee_studio::StudioState) logic layer
//! (spec → synthesis → ideal response → spec-mask verdict → dimensioning) so
//! the headless pipeline + its tests survive and remain reusable.
//!
//! The `[[bin]]` target therefore links to this stub `main`, which simply
//! points users at the Dioxus studio. Use the library's
//! [`StudioState`](yee_studio::StudioState) API to drive the headless flow.

fn main() {
    println!(
        "yee-studio is now a headless logic crate (the eframe view was retired in \
         App.D.2 / ADR-0130). The studio is `yee-studio-web` (Dioxus): \
         `dx serve` inside crates/yee-studio-web. Drive the headless flow via \
         the `yee_studio::StudioState` API."
    );
}
