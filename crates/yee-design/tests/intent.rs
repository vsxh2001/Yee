//! Property test: 100 randomly generated `DesignIntent` samples round-trip
//! through `serde_json::{to_string, from_str}` byte-identically (R1 DoD).
//!
//! Uses a deterministic linear-congruential generator so the test is
//! reproducible without pulling in `rand` (not in the workspace dev-dep set
//! at base SHA). Coverage spans:
//!   - both `Substrate::Named` (with and without `override_with`) and
//!     `Substrate::Explicit` arms
//!   - all `Option<f64>` fields populated / not populated
//!   - representative non-ASCII / multi-line `source_prompt` values

use yee_design::{
    DesignIntent, GeometryFamily, NamedSubstrate, Provenance, Substrate, SubstrateOverride,
    substrate_library,
};

/// Minimal deterministic LCG (Numerical Recipes constants). Deliberately not
/// cryptographic — we want reproducibility, not entropy.
struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        // glibc LCG constants.
        self.state = self.state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
        self.state
    }

    fn next_f64_in(&mut self, lo: f64, hi: f64) -> f64 {
        // Map the top 53 bits of state to [0, 1); scale to [lo, hi). Then
        // quantise to 6 significant digits via parse-after-format so the
        // value is in serde_json's ryu canonical form (its lexical parser
        // does not always reproduce a ryu emission bit-for-bit; quantising
        // here keeps the property test focused on the type surface rather
        // than on serde_json's `f64` parser quirks).
        let bits = self.next_u64() >> 11; // 53 random bits
        let u = (bits as f64) * (1.0 / ((1u64 << 53) as f64));
        let raw = lo + u * (hi - lo);
        // Snap to 6 significant decimal digits.
        let s = format!("{raw:.6e}");
        s.parse::<f64>().expect("ryu canonical f64")
    }

    fn next_choice<'a, T>(&mut self, items: &'a [T]) -> &'a T {
        let i = (self.next_u64() as usize) % items.len();
        &items[i]
    }

    fn next_bool(&mut self) -> bool {
        (self.next_u64() & 1) == 1
    }
}

fn gen_intent(rng: &mut Lcg, lib_version: &str) -> DesignIntent {
    let substrate_names = ["FR4", "RO4003C", "RO5880", "AluminaTC"];
    let prompts = [
        "2.4 GHz patch on FR4",
        "5.8 GHz patch on RO4003C with 200 MHz bandwidth",
        "915 MHz inset-fed patch on FR4, gain over 5 dBi",
        "design a 10 GHz patch on alumina",
        "patch antenna with non-ascii: αβγ ε_r ≈ 3.0",
        "multi\nline\nprompt with newlines and \"quotes\"",
    ];
    let sources = ["llm", "offline"];
    let models = [
        Some("claude-sonnet-4-5".to_string()),
        Some("claude-opus-4-7".to_string()),
        None,
    ];

    let substrate = if rng.next_bool() {
        let name = (*rng.next_choice(&substrate_names)).to_string();
        let override_with = if rng.next_bool() {
            Some(SubstrateOverride {
                eps_r: if rng.next_bool() {
                    Some(rng.next_f64_in(1.5, 12.0))
                } else {
                    None
                },
                h_mm: if rng.next_bool() {
                    Some(rng.next_f64_in(0.1, 5.0))
                } else {
                    None
                },
                loss_tangent: if rng.next_bool() {
                    Some(rng.next_f64_in(0.0, 0.05))
                } else {
                    None
                },
            })
        } else {
            None
        };
        Substrate::Named(NamedSubstrate {
            name,
            override_with,
        })
    } else {
        Substrate::Explicit {
            eps_r: rng.next_f64_in(1.5, 12.0),
            h_mm: rng.next_f64_in(0.1, 5.0),
            loss_tangent: rng.next_f64_in(0.0, 0.05),
        }
    };

    DesignIntent {
        family: GeometryFamily::RectangularPatch,
        target_frequency_hz: rng.next_f64_in(1.0e6, 1.0e12),
        substrate,
        gain_target_dbi: if rng.next_bool() {
            Some(rng.next_f64_in(-5.0, 30.0))
        } else {
            None
        },
        bandwidth_target_mhz: if rng.next_bool() {
            Some(rng.next_f64_in(0.0, 5_000.0))
        } else {
            None
        },
        source_prompt: (*rng.next_choice(&prompts)).to_string(),
        provenance: Provenance {
            source: (*rng.next_choice(&sources)).to_string(),
            model: rng.next_choice(&models).clone(),
            temperature: if rng.next_bool() {
                Some(rng.next_f64_in(0.0, 1.0))
            } else {
                None
            },
            schema_version: "1".to_string(),
            substrate_library_version: lib_version.to_string(),
        },
    }
}

#[test]
fn design_intent_serde_json_round_trips_byte_identically_over_100_samples() {
    // Spec §8 reproducibility invariant: the same `DesignIntent` value must
    // serialise to the same JSON bytes every time, and parsing those bytes
    // back must produce a value that serialises to those same bytes. This
    // is the contract the `<out>.intent.json` artefact rests on.
    //
    // Note on `f64` semantics: serde_json (via `ryu`) emits the shortest
    // string that round-trips to the same `f64` bit pattern, *but* the
    // second-pass parse is not required to be bit-identical to the original
    // arbitrary `f64` input — only to any other `f64` whose ryu shortest
    // representation is the same string. That is why we check the
    // "stabilises after one round-trip" property rather than bit-identity
    // against the original. After the first encode-decode the value is in
    // ryu's canonical form, and every subsequent encode is then truly
    // byte-identical — which is what spec §8 actually needs.
    let lib_version = substrate_library().version.clone();
    let mut rng = Lcg::new(0x59_45_45_44_53_47_4E_01); // "YEEDSGN\x01"
    for i in 0..100 {
        let intent = gen_intent(&mut rng, &lib_version);
        let json1 = serde_json::to_string(&intent).expect("serialize #1");
        let back: DesignIntent = serde_json::from_str(&json1).expect("deserialize");
        // Byte-identical re-serialisation: round-tripping the parsed value
        // through serde_json yields the same JSON string.
        let json2 = serde_json::to_string(&back).expect("serialize #2");
        assert_eq!(
            json1, json2,
            "sample {i}: re-serialisation diverged from first encode"
        );
        // Fixed-point: a second decode-encode pass is also stable.
        let back2: DesignIntent = serde_json::from_str(&json2).expect("deserialize #2");
        let json3 = serde_json::to_string(&back2).expect("serialize #3");
        assert_eq!(json2, json3, "sample {i}: second pass diverged");
        // And the value is also fixed-pointed (Eq on `back`/`back2`).
        assert_eq!(
            back, back2,
            "sample {i}: value not fixed-point after round-trip"
        );
    }
}
