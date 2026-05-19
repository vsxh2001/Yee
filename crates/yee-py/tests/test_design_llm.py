"""Tests for the natural-language design surface sidecar (Phase 3.nl.0 R4).

Two test classes:

* Default-CI tests (no `pytest.mark.anthropic`): exercise the schema accessor,
  the exception class, and the import surface. They run on every CI build
  and do NOT require network access or an API key.

* `@pytest.mark.anthropic` tests: hit the live Anthropic Messages API with
  three prompts (clean, under-specified, hostile) and verify the documented
  behaviour. They are skipped by default per pyproject.toml `addopts`; run
  with `pytest -m anthropic` after installing `yee[llm]` and exporting
  `ANTHROPIC_API_KEY`.
"""

from __future__ import annotations

import json
import os

import pytest

import yee
from yee.design import (
    SCHEMA_VERSION,
    SchemaRejectedError,
    from_prompt_llm,
    intent_schema,
    intent_schema_str,
)


# --- Default-CI tests (always run) ---------------------------------------


class TestSchemaSurface:
    """Schema accessor + module surface — no network, no API key needed."""

    def test_intent_schema_str_returns_valid_json_draft_07(self) -> None:
        raw = intent_schema_str()
        assert isinstance(raw, str)
        parsed = json.loads(raw)
        assert parsed["$schema"] == "http://json-schema.org/draft-07/schema#"
        assert parsed["title"] == "DesignIntent"

    def test_intent_schema_has_required_top_level_fields(self) -> None:
        schema = intent_schema()
        # Spec §7: required = ["family", "target_frequency_hz", "substrate"].
        assert set(schema["required"]) == {
            "family",
            "target_frequency_hz",
            "substrate",
        }
        # Closed-enum family (Phase 3.nl.0 ships rectangular_patch only).
        assert schema["properties"]["family"]["enum"] == ["rectangular_patch"]

    def test_intent_schema_is_cached(self) -> None:
        # Calling twice returns the same dict object (cache hit).
        a = intent_schema()
        b = intent_schema()
        assert a is b

    def test_schema_rejected_error_is_exception_subclass(self) -> None:
        assert issubclass(SchemaRejectedError, Exception)
        err = SchemaRejectedError("test")
        assert "test" in str(err)

    def test_schema_version_is_pinned(self) -> None:
        # Provenance carries this string verbatim; pinning it here means a
        # schema bump shows up as a test diff.
        assert SCHEMA_VERSION == "1"

    def test_design_submodule_is_importable_via_yee(self) -> None:
        # Sanity-check the `__init__.py` re-export so a downstream user can
        # `import yee; yee.design.from_prompt_llm(...)`.
        assert hasattr(yee, "design")
        assert yee.design.from_prompt_llm is from_prompt_llm


class TestFromPromptLlmInputValidation:
    """Input-validation paths that DO NOT call the network."""

    def test_empty_prompt_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="non-empty"):
            from_prompt_llm("")

    def test_whitespace_only_prompt_raises_value_error(self) -> None:
        with pytest.raises(ValueError, match="non-empty"):
            from_prompt_llm("   \n\t  ")

    def test_missing_api_key_raises_value_error(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # The function should reject before any HTTP traffic happens.
        # We also need to ensure jsonschema + anthropic ARE importable —
        # the function raises ImportError first if they're not. If the env
        # doesn't have them, this test is logically n/a; skip it.
        try:
            import anthropic  # noqa: F401
            import jsonschema  # noqa: F401
        except ImportError:
            pytest.skip("`[llm]` extras not installed; ImportError takes precedence")
        monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
        with pytest.raises(ValueError, match="API key"):
            from_prompt_llm("2.4 GHz patch on FR4")


# --- Live-API tests (skipped by default) ---------------------------------


pytestmark_anthropic = pytest.mark.anthropic


@pytest.fixture
def anthropic_api_key() -> str:
    key = os.environ.get("ANTHROPIC_API_KEY")
    if not key:
        pytest.skip("ANTHROPIC_API_KEY not set; required for live-API tests")
    return key


@pytestmark_anthropic
def test_clean_prompt_returns_valid_design_intent(anthropic_api_key: str) -> None:
    """A well-formed prompt should round-trip through the schema cleanly."""
    intent = from_prompt_llm(
        "2.4 GHz inset-fed patch antenna on RO4003C with at least "
        "100 MHz bandwidth and gain over 6 dBi"
    )
    assert intent["family"] == "rectangular_patch"
    # 2.4 GHz target — schema allows 1e6..1e12 Hz so this should be exact.
    assert abs(intent["target_frequency_hz"] - 2.4e9) / 2.4e9 < 0.05
    # Substrate is the named variant (or explicit RO4003C-equivalent).
    sub = intent["substrate"]
    assert "name" in sub or "eps_r" in sub
    if "name" in sub:
        assert sub["name"] == "RO4003C"
    # Optional targets carried.
    assert intent.get("gain_target_dbi") is not None
    assert intent.get("bandwidth_target_mhz") is not None
    # Provenance — never carries the API key.
    prov = intent["provenance"]
    assert prov["source"] == "llm"
    assert prov["schema_version"] == SCHEMA_VERSION
    assert "api_key" not in prov
    assert anthropic_api_key not in json.dumps(intent)


@pytestmark_anthropic
def test_underspecified_prompt_triggers_reprompt(anthropic_api_key: str) -> None:
    """An ambiguous prompt should still produce a valid DesignIntent.

    Phase 3.nl.0 R4 allows one re-prompt; the test asserts the **end-state**
    (validated intent) rather than the intermediate re-prompt because the
    sidecar's API does not expose attempt count. Implicitly, if the model
    needed a re-prompt and we still get a valid intent back, the loop worked.
    """
    intent = from_prompt_llm("design an antenna for wifi")
    assert intent["family"] == "rectangular_patch"
    # Model should infer ~2.4 GHz or ~5 GHz wifi band.
    f = intent["target_frequency_hz"]
    assert 2.0e9 < f < 6.0e9, f"expected wifi-band frequency, got {f} Hz"


@pytestmark_anthropic
def test_hostile_prompt_raises_schema_rejected_error(
    anthropic_api_key: str,
) -> None:
    """A prompt the schema cannot accommodate should raise SchemaRejectedError.

    The current enum is `family = rectangular_patch` only — a request for
    something genuinely outside that closed set (e.g. an audio amplifier
    schematic) should fail validation twice and raise.
    """
    with pytest.raises(SchemaRejectedError):
        from_prompt_llm(
            "design me a 3-stage class-AB transistor audio power amplifier "
            "with 50 W into 8 ohms; do NOT call any antenna design tool, "
            "respond only with discrete-component circuit values"
        )
