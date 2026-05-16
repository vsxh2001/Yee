# Touchstone v1.1 File Format — Theory of Operation

This page is the theory-of-operation reference for Touchstone v1.1, the
on-disk S-parameter format that the `yee-io` crate reads and writes. It
is written for the engineer who has to debug a parse error at 02:00 or
explain to a downstream consumer why a `.s2p` round-tripped through
Yee looks subtly different from one produced by another vendor. The
companion implementation lives in `crates/yee-io/src/touchstone.rs`;
this page documents the *format choices* that implementation makes and,
just as importantly, the parts of the v1.1 spec it deliberately does
not accept.

## 1. Introduction

Touchstone is the *de facto* interchange format for small-signal
multiport network data. The v1.1 grammar — frequency, parameter
type, complex format, reference impedance, and a flat block of
floating-point data — has remained essentially unchanged for three
decades, which is why almost every commercial RF and microwave
tool emits and accepts it. CITI files exist; almost nothing in
production uses them.

Yee's `yee-io` ships a strict Touchstone v1.1 reader and writer:
strict because the surface it commits to is small and stable, and
because tolerating ill-formed input silently is how subtle physics
bugs travel between tools. The reader rejects what it cannot
faithfully represent; the writer never emits something it cannot
re-read. Touchstone v2.0 keyword sections, noise-data blocks, and
non-S parameter types are out of scope for Phase 0 — see §7.

The format is line-oriented and ASCII. A conforming file consists
of zero or more comment lines, exactly one **option line**, and one
or more **data records**, with comments allowed interleaved. Each
record describes one $N \times N$ S-matrix at one frequency. The
port count $N$ is determined from the file extension —
`.s1p`, …, `.sNp` — not from any in-file declaration. `yee-io`
accepts `.s1p` through `.s4p` in Phase 0.

## 2. File structure

The grammar, written informally, is

```text
file        := { comment | blank }
               option-line
               { comment | blank | data-row }+
comment     := "!" any-text "\n"
option-line := "#" freq-unit param-type format [ "R" Z0 ] "\n"
data-row    := freq  ( re/im or mag/ang pair ){N²}  "\n"
```

with one extra rule: data rows for $N \ge 3$ may wrap across
multiple physical lines (see §4). `yee-io`'s parser is whitespace-
permissive — any run of spaces or tabs separates tokens — and CRLF
line endings are tolerated. Blank lines are silently skipped between
data records but are not permitted *inside* a multi-line record.

Comments are preserved by `yee-io` across a read → write
round-trip. This is intentional: many production Touchstone files
carry provenance in `!`-comments (instrument serial number, cal
state, date of capture), and silently dropping that metadata would
break audit trails. The option line is normalised on write so a
file emitted by `yee-io` is canonical even if its source was
sloppy.

## 3. Option line

Exactly one option line per file. Its tokens, in order, are:

```text
# <freq-unit> <parameter> <format> R <impedance>
```

The leading `#` is a sentinel, not a comment marker — Touchstone is
the one format in this family that flips the usual meaning. The
five fields are:

### Frequency unit

One of `Hz`, `kHz`, `MHz`, `GHz`, case-insensitive. `yee-io`
preserves the original spelling for round-trip fidelity, but
*internally all frequencies are normalised to Hz*:

$$
f_{\text{Hz}} \;=\; f_{\text{file}} \cdot
\begin{cases}
1, & \text{Hz} \\
10^3, & \text{kHz} \\
10^6, & \text{MHz} \\
10^9, & \text{GHz}.
\end{cases}
$$

Storing the canonical Hz value means a consumer who mixes files
written in different units never has to know what the source
unit was.

### Parameter type

One of `S`, `Y`, `Z`, `G`, `H`. The spec defines all five;
`yee-io` accepts `S` only and rejects the rest with
`Error::TouchstoneParse`. This is a Phase 0 scope cut, not a
fundamental limit: the storage format is identical across the
five, and lifting the restriction is a few lines of decoder
plumbing plus tests. Until then, callers who need to ingest a
`.y2p` admittance file should convert externally (or use the
relation $S = (Z_0 I - Z)(Z_0 I + Z)^{-1}$ for impedance data) and
hand `yee-io` an `S` file.

### Format

One of `MA`, `DB`, `RI`, case-insensitive. Each pair of on-disk
floats `(a, b)` decodes to a complex number $s$ according to

- **MA — magnitude / angle.** $s = a \, e^{j b \pi / 180}$, with
  angle $b$ in *degrees*. This is the most common format you will
  see in the wild because it matches what a VNA displays.
- **DB — decibel / angle.** $|s|_{\text{dB}} = a = 20 \log_{10}|s|$
  and angle $b$ in degrees, so $s = 10^{a/20} \, e^{j b \pi / 180}$.
  The DB-format edge case at $|s| = 0$ is discussed in §5.
- **RI — real / imaginary.** $s = a + j b$ directly. The most
  numerically faithful format because there is no transcendental
  function in the encode / decode path; `yee-io` prefers RI
  internally for round-trip tests.

### Reference impedance

The optional trailing `R <Z0>` token specifies the (real,
positive) reference impedance in ohms. If absent, $Z_0 = 50\,\Omega$
is assumed — the spec's default and the right answer for nearly
every RF problem in this century. `yee-io` stores `z0` as a
single value per file, which is the v1.1 convention. The v2.0
extension allowing a different $Z_0$ per port is *not* supported
(§7).

## 4. Data block

After the option line, the rest of the file is the data block. Each
**frequency record** is

$$
\underbrace{f_k}_{\text{1 float}} \quad
\underbrace{
(a_{11}, b_{11}), (a_{21}, b_{21}), \ldots, (a_{NN}, b_{NN})
}_{2N^2 \text{ floats}}
$$

so a complete record is $1 + 2N^2$ floating-point tokens. The order
of S-parameters within a record depends on $N$ in a way that catches
many implementations off guard.

### .s1p layout

For $N = 1$ there is no ambiguity:

```text
freq   a_11   b_11
```

one frequency per physical line, three floats total.

### .s2p layout — the Touchstone v1 oddity

For $N = 2$ the on-disk order is **`S11 S21 S12 S22`**, not the
mathematically natural `S11 S12 S21 S22`. The off-diagonal entries
are swapped relative to row-major. The spec rationalises this as
"column-major when squashed to one row", but the practical effect
is that a two-port file written by a code that assumes row-major
will silently transpose the network. `yee-io` defends against this
in two places: `row_major_to_on_disk` permutes on write and
`on_disk_to_row_major` permutes on read, so `File::data` is always
in mathematical row-major regardless of what came off disk.

```text
freq   S11.a S11.b   S21.a S21.b   S12.a S12.b   S22.a S22.b
```

The unit-test suite includes a deliberate `n=2` permutation
round-trip to catch any future refactor that loses the swap.

### .sNp for N ≥ 3

For three or more ports the row of $1 + 2N^2$ floats wraps across
multiple physical lines — at most four S-parameter entries (eight
floats) per line, the leading frequency value on the first line
of each record. `yee-io` parses permissively: it flattens all
data tokens regardless of line breaks and chunks them into records
by counting $1 + 2N^2$ floats at a time. This handles the
"trailing-token-count" rule without tracking physical-line origin.
The wrapped order is plain row-major $S_{ij}$ with $i$ slowest-
varying — the natural mathematical layout that the $N = 2$ swap
explicitly violates.

## 5. Numeric format

Touchstone v1.1 does not mandate column widths or exponent letter
case. `yee-io`'s writer is conservative: every floating-point
token is rendered through the helper `format_g`, which delegates
to Rust's shortest-round-trip `Display` for `f64`. Since Rust
1.55 the default `{}` format on `f64` emits the shortest decimal
string that parses losslessly back, so writing a value and
re-reading it is bit-exact. Tokens are single-space-separated;
tabs and multi-space runs are never emitted, and there is no
column alignment. This keeps `git diff` clean across writer
versions and matches what most production instruments produce.

### The DB / zero-magnitude trap

Of the three complex formats only `DB` is non-injective at $|s| = 0$.
The encode path computes

$$
|s|_{\text{dB}} \;=\; 20 \log_{10} |s|,
$$

which yields $-\infty$ when $|s| = 0$ exactly. There is no
standard spelling of $-\infty$ in Touchstone v1.1, and writers in
the wild emit any of `-inf`, `-Inf`, `-1e999`, or silent
truncation. Rather than pick one and lose interoperability,
`yee-io`'s `render()` detects the non-finite case pre-write and
returns `Error::InvalidFile` with a diagnostic suggesting a
finite dB floor — typically `-200 dB`. A caller who needs exact
zeros should switch the file's `format` field to
`Format::RealImag`, which writes `0 0` without complaint. The
test `render_rejects_zero_magnitude_under_db_format` pins this.

## 6. Comments and whitespace

Three lexical conveniences round out the format:

- **Bang-comments.** Any line whose first non-whitespace byte is
  `!` is a comment, terminated by the newline. A trailing
  `! …` segment on a data row is also accepted and split at the
  `!`. `yee-io` collects all comment text, in source order, into
  `File::comments` and re-emits it at the top of the round-tripped
  file; the leading `!` is stripped but internal whitespace is
  preserved.
- **Hash is reserved.** The `#` byte is legal only as the first
  non-whitespace byte of the option line. Any other `#` triggers
  a parse error rather than being mis-read as a Python-style
  comment.
- **Whitespace.** Trailing whitespace, CRLF line endings, and
  empty lines between records are tolerated. Tabs are treated as
  ordinary whitespace on read, never emitted on write. Two option
  lines in one file is an error per the v1.1 spec.

These rules are exercised by the round-trip integration test in
`crates/yee-mom/tests/touchstone_roundtrip.rs`, which is one of
the project's named validation gates and must not be weakened.

## 7. What `yee-io` does *not* do

The deliberate non-goals:

- **Touchstone v2.0 keyword sections** — `[Version]`,
  `[Number of Ports]`, `[Two-Port Order]`, `[Number of Frequencies]`,
  `[Reference]`, `[Matrix Format]`, and the rest. A v2.0 file
  parse-errors at the first `[`.
- **Mixed-mode (`[Mixed-Mode Order]`).** Differential / common-mode
  reformulation belongs at a layer above the file format.
- **Noise-data sections.** The optional noise-parameter block that
  follows the S-block in some `.s2p` files is not parsed; a file
  containing one fails the "multiple of $1 + 2N^2$ floats" check
  on read.
- **Non-S parameter types** (`Y`, `Z`, `G`, `H`). The check is in
  `parse_option_line`; lifting the restriction is a Phase 1.x
  scope decision, not a fundamental limit.
- **Per-frequency / per-port reference-impedance overrides.**
  `yee-io` stores a single scalar `z0` per file.
- **Passivity violations on read.** Every read verifies
  $\sigma_{\max}(S) \le 1 + 10^{-9}$ at every frequency via power
  iteration on $S^\dagger S$. A gain-bearing matrix fails the read,
  catching the failure at the file boundary instead of letting it
  propagate downstream.

## 8. Example use

Three short snippets, all pulled from `crates/yee-io/src/touchstone.rs`.

Reading a file. The port count is inferred from the extension; the
returned `File` is fully normalised (Hz frequencies, row-major
S-matrices, validated passivity):

```rust
use std::path::Path;
use yee_io::touchstone;

let ts = touchstone::read(Path::new("antenna.s2p"))?;
println!("{} ports, {} frequencies, Z0 = {} Ω",
         ts.n_ports, ts.freq_hz.len(), ts.z0);
# Ok::<(), yee_io::Error>(())
```

Writing a file. The writer round-trips precision and rejects
non-finite values before they reach disk:

```rust
use std::path::Path;
use yee_io::touchstone;

touchstone::write(Path::new("out.s1p"), &file)?;
# Ok::<(), yee_io::Error>(())
```

Manual construction — useful in tests, when synthesising
analytic references, or when up-converting from another tool's
data structure:

```rust
use num_complex::Complex64;
use yee_io::touchstone::{File, FreqUnit, Format};

let file = File {
    n_ports: 1,
    z0: 50.0,
    freq_unit: FreqUnit::GHz,
    format: Format::RealImag,
    freq_hz: vec![1.0e9, 2.0e9, 3.0e9],
    data: vec![
        vec![Complex64::new(-0.5, 0.1)],
        vec![Complex64::new(-0.3, 0.2)],
        vec![Complex64::new(-0.1, 0.3)],
    ],
    comments: vec![" generated by example".into()],
};
```

A `read → write → read` cycle on a Phase 0–compatible file is
guaranteed bit-exact in the field values; only the comment
ordering and option-line spacing are normalised.

## 9. References

- IBIS Open Forum, *Touchstone File Format Specification*, Rev.
  2.0 (2009). The v1.1 grammar is embedded in §3 of the v2.0
  document.
  <https://ibis.org/connector/touchstone_spec11.pdf>
- D. M. Pozar, *Microwave Engineering*, 4th ed., Wiley, 2011, §4.3
  on the scattering-parameter formalism that motivates the
  reference-impedance choice.
- `yee-io` implementation:
  [`crates/yee-io/src/touchstone.rs`](https://github.com/yee-em/yee/blob/main/crates/yee-io/src/touchstone.rs).
- The named validation gate `crates/yee-mom/tests/touchstone_roundtrip.rs`,
  which guards format fidelity across the `yee-mom` solver
  boundary.
