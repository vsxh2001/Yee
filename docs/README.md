# Yee — Documentation

User and developer documentation for Yee. Source for the eventual mdBook site lives here.

## Layout

```
docs/
├── README.md          # this file
├── source/            # archived input artifacts (initial content package, design notes)
├── theory/            # theory-of-operation pages, one per solver, with citations
│   ├── planar-mom.md
│   ├── fdtd.md
│   └── surrogates.md
├── tutorials/         # hands-on guides
│   ├── 01-microstrip-line.md
│   ├── 02-patch-antenna.md
│   └── 03-hairpin-filter.md
├── reference/         # API reference (auto-generated rustdoc + pydoc rolls up here)
└── decisions/         # ADRs: dependency choices, licensing, design tradeoffs
```

## Conventions

- Every theory page **cites its sources** — paper title, author, year, DOI where possible.
- Every tutorial is **runnable end-to-end** with a fixed Yee version pin at the top.
- ADRs (Architecture Decision Records) follow the Michael Nygard one-page template: Context, Decision, Status, Consequences.

## Building the site (Phase 1+)

```bash
cargo install mdbook
mdbook serve docs/
```

CI publishes to GitHub Pages on every push to `main`.

## Phase 0 status

Stub. Source artifact archived under `docs/source/`. Tutorials and theory pages arrive incrementally as the corresponding solver features ship.
