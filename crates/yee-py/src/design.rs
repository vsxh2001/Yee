//! `yee.design` submodule — natural-language design surface (Phase 3.nl.0 R4).
//!
//! The submodule exposes two artefacts to Python:
//!
//! 1. `intent_schema_str() -> str` — the Draft-07 JSON schema baked into the
//!    `yee-design` crate via `include_str!`. The Python sidecar
//!    (`crates/yee-py/python/yee/design.py`) loads this string, deserialises
//!    it with `json.loads`, and uses it both as the Anthropic Messages-API
//!    tool's `input_schema` and as the `jsonschema.validate` argument for
//!    the returned `tool_use.input`. Sharing one schema string across both
//!    sides of the boundary is the spec §7 / plan R4 contract.
//! 2. `SchemaRejectedError` — a custom Python exception class raised by
//!    `yee.design.from_prompt_llm` when the LLM's tool-use response fails
//!    JSON-schema validation **twice** (once on the initial call, once after
//!    a re-prompt with the schema embedded in the user turn — see
//!    `yee.design.from_prompt_llm` in `python/yee/design.py`).
//!
//! The actual Anthropic Messages-API call lives entirely in Python; the Rust
//! side stays free of any LLM-client dependency (spec §11 / plan §"Tech-stack
//! additions"). This keeps the `cargo build` path green without network or
//! API keys, and keeps the `anthropic` SDK confined to the
//! `[project.optional-dependencies] llm` extra in `crates/yee-py/pyproject.toml`.
//!
//! ## Submodule registration
//!
//! Follows the same pattern as `yee.touchstone`, `yee.eigensolver`, and
//! `yee.fem`: `lib.rs` inserts this module into `sys.modules` so both
//! `from yee.design import intent_schema_str` and the attribute-access form
//! `yee.design.intent_schema_str` succeed.

use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;

use yee_design::intent::INTENT_SCHEMA;

// `SchemaRejectedError` lives at `yee.design.SchemaRejectedError`. Use
// `create_exception!` so the class is a real Python `Exception` subclass —
// `isinstance(err, Exception)` and `try: ... except yee.design.SchemaRejectedError`
// both work without further glue. The macro generates a `pub struct
// SchemaRejectedError;` in the local module along with the necessary
// `IntoPyErr` impl.
create_exception!(
    yee.design,
    SchemaRejectedError,
    PyException,
    "Raised when the Anthropic Messages-API tool-use response fails JSON-schema validation \
     twice in a row (once on the initial call, once on the re-prompt). The exception \
     carries the validation error message from the second failure; the offending JSON \
     itself is logged at debug level by the Python sidecar but is NOT included in the \
     exception payload, since it may contain partial prompts."
);

/// Register the `yee.design` submodule.
///
/// Adds the `SchemaRejectedError` exception class and the
/// `intent_schema_str` function. `from_prompt_llm` itself is implemented in
/// pure Python (see `crates/yee-py/python/yee/design.py`) — the Rust side
/// only exposes the schema string and the exception type.
pub fn register(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("SchemaRejectedError", py.get_type::<SchemaRejectedError>())?;
    m.add_function(wrap_pyfunction!(intent_schema_str, m)?)?;
    Ok(())
}

/// Return the Draft-07 JSON schema string for `DesignIntent`.
///
/// The string is the verbatim contents of
/// `crates/yee-design/src/intent_schema.json` baked at compile time. The
/// Python sidecar calls `json.loads(intent_schema_str())` once at import time
/// and re-uses the parsed `dict` as both the Anthropic tool's `input_schema`
/// and the `jsonschema.validate` argument.
#[pyfunction]
fn intent_schema_str() -> &'static str {
    INTENT_SCHEMA
}
