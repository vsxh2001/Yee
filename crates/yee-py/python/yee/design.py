"""Natural-language design surface — Python sidecar (Phase 3.nl.0 R4).

Hosts ``yee.design.from_prompt_llm``, a one-shot wrapper around the Anthropic
Messages API that turns a free-form prompt into a ``DesignIntent``-shaped
``dict`` validated against the Draft-07 JSON schema baked into the
``yee-design`` crate (see ``crates/yee-design/src/intent_schema.json``).

The Rust side (``crates/yee-py/src/design.rs``) exposes only the schema
string and the :class:`SchemaRejectedError` exception class; this module
contains the LLM-facing logic so the Rust ``cargo build`` path stays
LLM-client-free per spec §11.

The ``anthropic`` and ``jsonschema`` dependencies are declared in
``crates/yee-py/pyproject.toml`` under ``[project.optional-dependencies] llm``
and are imported lazily so ``import yee`` succeeds on a default install
without them. Calling :func:`from_prompt_llm` without the extras raises
:class:`ImportError` with installation guidance.
"""

from __future__ import annotations

import json
import logging
import os
from typing import Any

# Re-export the Rust-side exception + schema accessor so the public surface
# lives at `yee.design.SchemaRejectedError` and `yee.design.intent_schema_str`.
# `intent_schema_str` is a function (not a constant) because the schema is a
# `&'static str` baked into the extension module — exposing a function keeps
# the Python contract symmetric with the Rust one and avoids freezing a copy
# in the package namespace at import time.
from yee._yee.design import (  # noqa: F401  re-export
    SchemaRejectedError,
    intent_schema_str,
)

__all__ = [
    "DEFAULT_MODEL",
    "SCHEMA_VERSION",
    "SchemaRejectedError",
    "from_prompt_llm",
    "intent_schema",
    "intent_schema_str",
]


_log = logging.getLogger("yee.design")

#: Schema version recorded in ``DesignIntent.provenance.schema_version`` (spec §7).
SCHEMA_VERSION: str = "1"

#: Default Anthropic model when the caller does not pass ``model=``.
#:
#: Pinned to a current Sonnet release; override with the ``model=`` argument
#: or the ``YEE_ANTHROPIC_MODEL`` environment variable.
DEFAULT_MODEL: str = "claude-sonnet-4-5-20250929"

#: Tool name exposed to the Anthropic Messages API (plan R4).
_TOOL_NAME: str = "emit_design_intent"

#: System prompt — pinned, deterministic. Loosens nothing the schema enforces.
_SYSTEM_PROMPT: str = (
    "You are an antenna-design intake assistant. Convert the user's "
    "natural-language design request into a single call to the "
    f"`{_TOOL_NAME}` tool. Infer reasonable defaults (e.g. 'patch' "
    "implies family='rectangular_patch') but NEVER invent a target "
    "frequency, substrate, gain, or bandwidth that the user did not "
    "state or strongly imply. If the request is too ambiguous to fill "
    "the required fields, ask the user no questions and instead call "
    "the tool with your best inference — the downstream validator will "
    "reject obviously wrong values."
)


def intent_schema() -> dict[str, Any]:
    """Return the parsed JSON schema as a Python ``dict``.

    Identical to ``json.loads(intent_schema_str())``; cached on first call.
    """
    cache = intent_schema.__dict__
    if "schema" not in cache:
        cache["schema"] = json.loads(intent_schema_str())
    return cache["schema"]


def _require_anthropic_and_jsonschema() -> tuple[Any, Any]:
    """Lazy-import the `[llm]` extras, raising ImportError if missing.

    Returns ``(anthropic_module, jsonschema_module)``.
    """
    try:
        import anthropic  # type: ignore[import-not-found]
        import jsonschema  # type: ignore[import-not-found]
    except ImportError as exc:  # pragma: no cover — environment-dependent
        raise ImportError(
            "yee.design.from_prompt_llm requires the optional 'llm' extras. "
            "Install with: pip install 'yee[llm]'  (or  pip install anthropic jsonschema)"
        ) from exc
    return anthropic, jsonschema


def _validate(payload: Any, jsonschema_mod: Any) -> tuple[bool, str]:
    """Run jsonschema.validate; return ``(ok, error_message)``.

    Never raises. The error message is the validator's stringification; the
    raw payload is intentionally NOT returned upward (it may contain partial
    user prompt fragments) — it's only logged at debug level here.
    """
    try:
        jsonschema_mod.validate(instance=payload, schema=intent_schema())
    except jsonschema_mod.ValidationError as exc:
        _log.debug("DesignIntent validation rejected payload: %r", payload)
        return False, str(exc)
    return True, ""


def _extract_tool_input(message: Any) -> Any | None:
    """Walk a Messages-API response and return the first ``tool_use.input``.

    The Anthropic SDK returns a `Message` whose `.content` is a list of typed
    blocks. We scan for the first `ToolUseBlock` with name matching
    :data:`_TOOL_NAME`; if none, return ``None`` so the caller can surface
    the right error.

    The block-shape access pattern (``block.type == "tool_use"`` with
    ``block.input`` as a dict) is the documented stable surface as of
    ``anthropic>=0.40``. If a future SDK reshapes this, the offline path in
    R5 remains available as the deterministic fallback (spec §11) — so this
    function's failure mode is a clean ``None`` and a re-prompt, never an
    AttributeError into user code.
    """
    content = getattr(message, "content", None)
    if not content:
        return None
    for block in content:
        if getattr(block, "type", None) != "tool_use":
            continue
        if getattr(block, "name", None) != _TOOL_NAME:
            continue
        return getattr(block, "input", None)
    return None


def _build_tool() -> dict[str, Any]:
    """Construct the Anthropic ``tools`` array entry for the design tool."""
    return {
        "name": _TOOL_NAME,
        "description": (
            "Emit a structured DesignIntent for one antenna design "
            "request. Call this exactly once per user message. The "
            "input object must conform to the Yee DesignIntent v1 "
            "schema."
        ),
        "input_schema": intent_schema(),
    }


def _call_messages(
    client: Any,
    *,
    model: str,
    user_text: str,
    tool: dict[str, Any],
) -> Any:
    """Single Anthropic Messages call with tool-use forcing.

    `tool_choice = {"type": "tool", "name": _TOOL_NAME}` forces the model to
    invoke the tool. `temperature=0.0` makes the call as deterministic as the
    API allows; the actual value is recorded in `DesignIntent.provenance`
    by :func:`from_prompt_llm` so a future re-validation against a moved
    model is detectable.
    """
    return client.messages.create(
        model=model,
        max_tokens=1024,
        temperature=0.0,
        system=_SYSTEM_PROMPT,
        tools=[tool],
        tool_choice={"type": "tool", "name": _TOOL_NAME},
        messages=[{"role": "user", "content": user_text}],
    )


def from_prompt_llm(
    prompt: str,
    *,
    model: str | None = None,
    api_key: str | None = None,
) -> dict[str, Any]:
    """Convert a natural-language design prompt into a validated DesignIntent dict.

    Parameters
    ----------
    prompt :
        Free-form English design request. Example:
        ``"2.4 GHz inset-fed patch on RO4003C with at least 100 MHz bandwidth"``.
    model :
        Anthropic model id. Defaults to :data:`DEFAULT_MODEL` or the
        ``YEE_ANTHROPIC_MODEL`` env var if set.
    api_key :
        Anthropic API key. If ``None``, falls back to the ``ANTHROPIC_API_KEY``
        environment variable. The key is **never** logged and **never**
        included in the returned ``provenance`` block (spec §10 risk #3).

    Returns
    -------
    dict
        A ``DesignIntent``-shaped dict with the LLM's structured fields plus
        the caller-side ``source_prompt`` and ``provenance`` metadata. The
        returned dict is suitable for ``serde_json::from_str::<DesignIntent>``
        on the Rust side.

    Raises
    ------
    ImportError
        If the ``[llm]`` extras (``anthropic``, ``jsonschema``) are not
        installed.
    ValueError
        If no API key is available (neither ``api_key=`` nor the env var).
    SchemaRejectedError
        If the LLM's tool-use response fails schema validation twice (once
        on the initial call, once after a re-prompt with the schema embedded
        in the user turn). The exception message carries the validator's
        second-failure diagnostic; the raw payload is **not** included in
        the exception text since it may contain partial user prompt
        fragments.
    """
    if not prompt or not prompt.strip():
        raise ValueError("prompt must be a non-empty string")

    anthropic_mod, jsonschema_mod = _require_anthropic_and_jsonschema()

    resolved_key = api_key if api_key is not None else os.environ.get("ANTHROPIC_API_KEY")
    if not resolved_key:
        raise ValueError(
            "no Anthropic API key — pass api_key=... or set ANTHROPIC_API_KEY"
        )
    resolved_model = model or os.environ.get("YEE_ANTHROPIC_MODEL") or DEFAULT_MODEL

    # Build the client locally; never store the key on a module-level object.
    client = anthropic_mod.Anthropic(api_key=resolved_key)
    tool = _build_tool()

    # --- Attempt 1 -------------------------------------------------------
    message = _call_messages(
        client, model=resolved_model, user_text=prompt, tool=tool
    )
    payload = _extract_tool_input(message)
    if payload is not None:
        ok, err = _validate(payload, jsonschema_mod)
        if ok:
            return _finalise(payload, prompt=prompt, model=resolved_model)
        first_error = err
    else:
        first_error = (
            "model did not call the `emit_design_intent` tool on the first turn"
        )
    _log.info("first attempt rejected: %s — re-prompting with schema appended", first_error)

    # --- Attempt 2: append the schema to the user turn -------------------
    reprompt = (
        f"{prompt}\n\n"
        f"Your previous response was rejected by the validator: {first_error}\n"
        f"You MUST call the `{_TOOL_NAME}` tool with an `input` object that "
        "exactly matches this JSON schema:\n\n"
        f"```json\n{json.dumps(intent_schema(), indent=2)}\n```"
    )
    message = _call_messages(
        client, model=resolved_model, user_text=reprompt, tool=tool
    )
    payload = _extract_tool_input(message)
    if payload is None:
        raise SchemaRejectedError(
            "model did not call the `emit_design_intent` tool on the re-prompt either"
        )
    ok, err = _validate(payload, jsonschema_mod)
    if not ok:
        raise SchemaRejectedError(
            f"DesignIntent JSON-schema validation failed twice; second error: {err}"
        )
    return _finalise(payload, prompt=prompt, model=resolved_model)


def _finalise(payload: dict[str, Any], *, prompt: str, model: str) -> dict[str, Any]:
    """Merge the LLM's structured fields with caller-side metadata.

    The returned dict is a complete ``DesignIntent`` (per the Rust struct's
    serde shape): the LLM-supplied fields plus the verbatim ``source_prompt``
    and a ``provenance`` block. The provenance never carries the API key —
    only the model id, temperature, schema version, and a placeholder for
    the substrate-library version (which the Rust ``DesignIntent`` resolver
    fills in once the named-substrate lookup happens; see spec §10 risk #2).
    """
    # `payload` is the validated input object; copy it so we never mutate
    # the dict the caller hands us by reference downstream.
    intent: dict[str, Any] = dict(payload)
    intent["source_prompt"] = prompt
    intent["provenance"] = {
        "source": "llm",
        "model": model,
        "temperature": 0.0,
        "schema_version": SCHEMA_VERSION,
        # The library version is filled in by the Rust resolver when a
        # `Substrate::Named` lookup runs; the sidecar leaves an empty string
        # rather than guessing. Spec §10 risk #2 covers the drift detection
        # path; the Rust side replaces this with the resolved version.
        "substrate_library_version": "",
    }
    return intent
