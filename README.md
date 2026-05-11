# REM — Rust Electromagnetic Solver

[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20%2F%20Apache--2.0-blue.svg)](#licensing)
[![Rust](https://img.shields.io/badge/rust-stable-orange.svg)](https://www.rust-lang.org)

**REM** is an open-source, Palace-compatible electromagnetic finite-element simulation toolkit written in pure Rust.  
It targets computational EM workflows where safety, portability, and WebAssembly deployment matter.

---

## Goals

- **Palace-compatible** — accepts the same JSON/YAML problem configs used by the Palace FEM solver
- **Solver breadth** — electrostatic, magnetostatic, eigenmode, frequency-domain driven, transient, and boundary-element problems in one workspace
- **Browser-ready** — the `wasm` crate compiles to WebAssembly and powers a Yew-based UI
- **Extensible split** — Community edition keeps a clean public API; advanced Pro capabilities (MoM, SBR, FEBI, ROM) build on top without forking

---

## Quick Start

### Install the CLI

```bash
git clone git@github.com:rem-rs/rem.git
cd rem
cargo build -p rem-cli --release
# binary lands at: target/release/rem
```

### Run a simulation

```bash
# Electrostatic simulation with a Palace-format config
./target/release/rem run --config examples/cpw/cpw_electrostatic.json

# Or with cargo during development
cargo run -p rem-cli -- run --config examples/cpw/cpw_electrostatic.json
```

### Build the WASM target

```bash
wasm-pack build crates/wasm --target web
```

---

## Capability Matrix

| Feature | Community (this repo) | Pro (private) |
|---|:---:|:---:|
| Electrostatic solver | ✓ | ✓ |
| Magnetostatic solver | ✓ | ✓ |
| Eigenmode solver | ✓ | ✓ |
| Frequency-domain driven | ✓ | ✓ |
| Transient solver | ✓ | ✓ |
| Boundary element (BEM) | ✓ | ✓ |
| WebAssembly / browser UI | ✓ | ✓ |
| Surface geometry primitives | ✓ | ✓ |
| Sonnet 19 / Palace convert | — | ✓ |
| 3-D MoM (full-wave) | — | ✓ |
| Planar full-wave solver | — | ✓ |
| FE-BI hybrid solver | — | ✓ |
| DDM domain decomposition | — | ✓ |
| SBR shooting-bouncing rays | — | ✓ |
| Driven ROM acceleration | — | ✓ |
| Advanced S-parameter analysis | — | ✓ |
| Touchstone matrix conversion | — | ✓ |
| Design / parametric optimization | — | ✓ |

---

## Crate Architecture

```
rem/
├── crates/
│   ├── core          — shared types, error, traits
│   ├── config        — Palace-format JSON/YAML deserialization
│   ├── mesh          — FEM mesh data structures
│   ├── surface       — shared BEM/MoM surface primitives & quadrature
│   ├── materials     — material property definitions
│   ├── bc            — boundary condition types
│   ├── electrostatic — electrostatic FEM solver
│   ├── magnetostatic — magnetostatic FEM solver
│   ├── eigenmode     — eigenmode FEM solver
│   ├── driven        — frequency-domain driven FEM solver
│   ├── transient     — time-domain transient FEM solver
│   ├── bem           — boundary element solver
│   ├── touchstone    — Touchstone SNP file I/O
│   ├── parallel      — MPI/Rayon parallel execution helpers
│   ├── result        — post-processing and field export
│   ├── cli           — Community command-line binary (`rem`)
│   ├── wasm          — WebAssembly entry point
│   └── yew-app       — Yew browser UI (excluded from default workspace)
└── vendor/
    ├── rmetis        — Rust bindings for METIS partitioner
    └── rmsh-*/       — Rust mesh I/O utilities
```

**Note:** Project format converters (Sonnet19, Ansys, ADS) have been moved to the private **rem-pro** workspace for exclusive Pro feature development.
```

---

## Examples

Example problem configs live in the parent workspace `examples/` directory.  
Topics covered: coplanar waveguide (CPW), near-field scanning, Palace benchmark problems, and Sonnet 19 round-trip validation.

---

## Licensing

REM Community Edition is dual-licensed under either of the following, at your option:

- [Apache License, Version 2.0](LICENSE-APACHE)
- [MIT License](LICENSE-MIT)

SPDX-License-Identifier: `MIT OR Apache-2.0`
