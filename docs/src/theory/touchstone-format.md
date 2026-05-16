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
multiport network data. It originated at EEsof in the late 1980s,
passed through Hewlett-Packard / Agilent / Keysight, and is now
maintained by the IBIS Open Forum. The v1.1 grammar — frequency,
parameter type, complex format, reference impedance, and a flat block
of floating-point data — has remained essentially unchanged across
those three decades, which is why almost every commercial RF and
microwave tool emits and accepts it. CITI files exist; almost nothing
in production uses them.

Yee's `yee-io` ships a strict Touchstone v1.1 reader and writer:
strict because the surface it commits to is small and stable, and
because tolerating ill-formed input silently is how subtle physics
bugs travel between tools. The reader rejects what it cannot
faithfully represent; the writer never emits something it cannot
re-read. Touchstone v2.0 keyword sections, noise-data blocks,
per-frequency reference-impedance overrides, and Y/Z/G/H parameter
types are out of scope for Phase 0 — see §7 for the explicit
non-goals list.

The format is line-oriented and ASCII. A conforming file consists of,
in order: zero or more comment lines, exactly one **option line**,
and one or more **data records**, with comments allowed interleaved.
A record describes the full $N \times N$ S-matrix at one frequency.
The port count $N$ is determined from the file extension —
`.s1p`, `.s2p`, …, `.sNp` — not from any in-file declaration.
`yee-io` accepts `.s1p` through `.s4p` in Phase 0.

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
preserves the original spelling so a `read → write` round-trip
emits the same unit string, but *internally all frequencies are
normalised to Hz*. The conversion is

$$
f_{\text{Hz}} \;=\; f_{\text{file}} \cdot
\begin{cases}
1, & \text{unit} = \text{Hz} \\
10^3, & \text{unit} = \text{kHz} \\
10^6, & \text{unit} = \text{MHz} \\
10^9, & \text{unit} = \text{GHz}
\end{cases}
$$

Storing the canonical Hz value means that a downstream consumer who
mixes files written in different units — a not-uncommon case — never
has to know what the original file's unit was.

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

For three or more ports the row of $1 + 2N^2$ floats is wrapped
across multiple physical lines, with **at most four S-parameter
entries (eight floats) per line** per the spec, the leading
frequency value on the first line of each record. `yee-io` parses
permissively — it flattens all data tokens regardless of line
breaks, then chunks them into records by counting $1 + 2N^2$
floats at a time. This handles the "trailing-token-count" rule
without requiring the parser to track which physical line a token
came from. The order within the wrapped record is plain row-major
`S_{ij}` with $i$ slowest-varying, the natural mathematical layout
that the $N = 2$ swap explicitly violates.
