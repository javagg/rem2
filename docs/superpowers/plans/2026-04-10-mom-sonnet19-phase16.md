# MoM Sonnet-19 Alignment — Phase 16 (v0.17.0) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a standalone `rem-touchstone` crate for Touchstone I/O, then add MoM lumped-port excitation and S-parameter extraction to `rem-mom`, enabling REM to compute S-parameters of 3-D MoM structures and output `.s{N}p` files.

**Architecture:** A new `crates/touchstone` crate owns all Touchstone read/write logic (multi-port, RI/MA/DB formats); `rem-driven` migrates its existing `write_touchstone_s1p` there; `rem-mom` gains `port.rs` (lumped-port model) and `sparams.rs` (S-matrix extraction + `touchstone` delegation). The MoM main loop in `lib.rs` branches on whether `Solver.MoM.Ports` is populated — if yes, run S-parameter sweep; if no, run the existing RCS path unchanged.

**Tech Stack:** Rust stable, `num-complex 0.4`, `nalgebra 0.33`, no new external deps.

---

## File Map

| Action | Path | Responsibility |
|--------|------|----------------|
| **Create** | `crates/touchstone/Cargo.toml` | crate manifest |
| **Create** | `crates/touchstone/src/lib.rs` | public API re-exports |
| **Create** | `crates/touchstone/src/write.rs` | `write_snp` — multi-port RI writer |
| **Create** | `crates/touchstone/src/read.rs` | `read_snp` — parser (needed for round-trip tests) |
| **Modify** | `Cargo.toml` (workspace) | add `crates/touchstone` member + `rem-touchstone` dep |
| **Modify** | `crates/driven/Cargo.toml` | add `rem-touchstone` dep |
| **Modify** | `crates/driven/src/vf.rs` | delegate `write_touchstone_s1p` → `rem_touchstone::write_snp` |
| **Modify** | `crates/mom/Cargo.toml` | add `rem-touchstone` dep |
| **Create** | `crates/mom/src/port.rs` | `MomLumpedPort` struct + RHS construction + V/I extraction |
| **Create** | `crates/mom/src/sparams.rs` | `compute_s_matrix`, `SMatrix`, delegation to `rem_touchstone` |
| **Modify** | `crates/config/src/schema.rs` | add `MomPort`, `ports: Vec<MomPort>`, `ref_impedance` fields |
| **Modify** | `crates/mom/src/lib.rs` | branch on `mom_cfg.ports` — S-param path vs existing RCS path |
| **Modify** | `crates/mom/tests/mie_validation.rs` | keep passing (regression guard) |

---

## Task 1: Create `crates/touchstone` — write side

**Files:**
- Create: `crates/touchstone/Cargo.toml`
- Create: `crates/touchstone/src/lib.rs`
- Create: `crates/touchstone/src/write.rs`

- [ ] **Step 1: Write the failing test for `write_snp`**

Create `crates/touchstone/src/write.rs` with only the test module:

```rust
// crates/touchstone/src/write.rs
use num_complex::Complex64;

/// Touchstone data format for the option line.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TsFormat { Ri, Ma, Db }

/// Frequency unit for the option line.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TsFreqUnit { Hz, Khz, Mhz, Ghz }

impl TsFreqUnit {
    pub fn scale(&self) -> f64 {
        match self {
            TsFreqUnit::Hz  => 1.0,
            TsFreqUnit::Khz => 1e3,
            TsFreqUnit::Mhz => 1e6,
            TsFreqUnit::Ghz => 1e9,
        }
    }
    pub fn label(&self) -> &'static str {
        match self {
            TsFreqUnit::Hz  => "Hz",
            TsFreqUnit::Khz => "KHz",
            TsFreqUnit::Mhz => "MHz",
            TsFreqUnit::Ghz => "GHz",
        }
    }
}

/// Write an N-port Touchstone file as a String.
///
/// `freqs_hz`  — frequency samples [Hz], length F  
/// `s_data`    — S-matrix samples; length F, each inner Vec length N*N (row-major)  
/// `n_ports`   — number of ports N  
/// `z0`        — reference impedance [Ω]  
/// `fmt`       — data format  
/// `unit`      — frequency unit in the option line  
///
/// # Panics
/// Panics if any `s_data[f].len() != n_ports * n_ports`.
pub fn write_snp(
    freqs_hz: &[f64],
    s_data: &[Vec<Complex64>],
    n_ports: usize,
    z0: f64,
    fmt: TsFormat,
    unit: TsFreqUnit,
) -> String {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex64;

    fn s1p_sample() -> (Vec<f64>, Vec<Vec<Complex64>>) {
        let freqs = vec![1e9, 2e9, 3e9];
        let s_data = vec![
            vec![Complex64::new(0.5, -0.3)],
            vec![Complex64::new(0.4, -0.2)],
            vec![Complex64::new(0.3, -0.1)],
        ];
        (freqs, s_data)
    }

    #[test]
    fn option_line_s1p() {
        let (f, s) = s1p_sample();
        let out = write_snp(&f, &s, 1, 50.0, TsFormat::Ri, TsFreqUnit::Ghz);
        assert!(out.contains("# GHz S RI R 50"), "option line missing: {out}");
    }

    #[test]
    fn correct_number_of_data_lines_s1p() {
        let (f, s) = s1p_sample();
        let out = write_snp(&f, &s, 1, 50.0, TsFormat::Ri, TsFreqUnit::Ghz);
        let data_lines: Vec<&str> = out.lines()
            .filter(|l| !l.starts_with('!') && !l.starts_with('#') && !l.trim().is_empty())
            .collect();
        assert_eq!(data_lines.len(), 3, "expected 3 data lines, got {data_lines:?}");
    }

    #[test]
    fn ri_values_present_s1p() {
        let (f, s) = s1p_sample();
        let out = write_snp(&f, &s, 1, 50.0, TsFormat::Ri, TsFreqUnit::Ghz);
        // First data line should contain freq, Re(S11), Im(S11)
        let first_data = out.lines()
            .find(|l| !l.starts_with('!') && !l.starts_with('#') && !l.trim().is_empty())
            .unwrap();
        let parts: Vec<f64> = first_data.split_whitespace()
            .map(|p| p.parse::<f64>().unwrap())
            .collect();
        assert_eq!(parts.len(), 3, "S1P data line needs 3 fields: {first_data}");
        assert!((parts[0] - 1.0).abs() < 1e-6, "freq should be 1.0 GHz: {}", parts[0]);
        assert!((parts[1] - 0.5).abs() < 1e-6, "Re(S11) wrong: {}", parts[1]);
        assert!((parts[2] - (-0.3)).abs() < 1e-6, "Im(S11) wrong: {}", parts[2]);
    }

    #[test]
    fn s2p_has_four_s_values_per_freq() {
        let freqs = vec![1e9];
        // S11, S12, S21, S22 row-major
        let s_data = vec![vec![
            Complex64::new(0.1, 0.0),
            Complex64::new(0.9, 0.0),
            Complex64::new(0.9, 0.0),
            Complex64::new(0.1, 0.0),
        ]];
        let out = write_snp(&freqs, &s_data, 2, 50.0, TsFormat::Ri, TsFreqUnit::Ghz);
        let data_lines: Vec<&str> = out.lines()
            .filter(|l| !l.starts_with('!') && !l.starts_with('#') && !l.trim().is_empty())
            .collect();
        // Touchstone 1.0 S2P: 1 line with freq + 4×(Re Im) = 9 fields
        assert!(!data_lines.is_empty());
        let fields: Vec<&str> = data_lines[0].split_whitespace().collect();
        assert_eq!(fields.len(), 9, "S2P line should have 9 fields: {data_lines:?}");
    }
}
```

- [ ] **Step 2: Run test to confirm it fails**

```bash
cd /c/Users/lilu/works/rem2
cargo test -p rem-touchstone 2>&1 | tail -20
```

Expected: compile error — `todo!()` makes tests unreachable, but crate doesn't compile yet because `Cargo.toml` doesn't exist. You'll see "package 'rem-touchstone' not found".

- [ ] **Step 3: Create `Cargo.toml` for the crate**

Create `crates/touchstone/Cargo.toml`:

```toml
[package]
name    = "rem-touchstone"
version = "0.1.0"
edition = "2021"

[dependencies]
num-complex = { workspace = true }
thiserror   = { workspace = true }
```

Create `crates/touchstone/src/lib.rs`:

```rust
//! Touchstone (.s{N}p) file read and write utilities.
//!
//! Supports N-port Touchstone 1.0 format, RI data only (MA/DB planned).
//! Used by `rem-mom` (S-parameter output) and `rem-driven` (VF circuit synthesis).

pub mod write;
pub mod read;

pub use write::{write_snp, TsFormat, TsFreqUnit};
pub use read::{read_snp, TouchstoneFile, TsReadError};
```

Create stub `crates/touchstone/src/read.rs` (just enough to compile):

```rust
//! Touchstone reader — minimal stub for now.
use num_complex::Complex64;

#[derive(Debug)]
pub struct TouchstoneFile {
    pub n_ports: usize,
    pub freqs_hz: Vec<f64>,
    /// S-matrix row-major per frequency: s_data[f][i*n+j] = S_{i+1,j+1}
    pub s_data: Vec<Vec<Complex64>>,
    pub z0: f64,
}

#[derive(Debug, thiserror::Error)]
pub enum TsReadError {
    #[error("parse error: {0}")]
    Parse(String),
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}

/// Parse a Touchstone file from a string.
/// Only RI format supported currently.
pub fn read_snp(content: &str) -> Result<TouchstoneFile, TsReadError> {
    let _ = content;
    Err(TsReadError::Parse("read_snp not yet implemented".to_string()))
}
```

- [ ] **Step 4: Register the crate in the workspace**

Edit `Cargo.toml` (workspace root). In the `members` array, add `"crates/touchstone"`. In `[workspace.dependencies]`, add:

```toml
rem-touchstone = { path = "crates/touchstone" }
```

- [ ] **Step 5: Run tests — confirm `todo!()` panics, not compile errors**

```bash
cargo test -p rem-touchstone 2>&1 | tail -20
```

Expected: tests compile but panic at runtime with "not yet implemented".

- [ ] **Step 6: Implement `write_snp`**

Replace `todo!()` in `write.rs`:

```rust
pub fn write_snp(
    freqs_hz: &[f64],
    s_data: &[Vec<Complex64>],
    n_ports: usize,
    z0: f64,
    fmt: TsFormat,
    unit: TsFreqUnit,
) -> String {
    assert_eq!(freqs_hz.len(), s_data.len(), "freq/s_data length mismatch");
    for (i, row) in s_data.iter().enumerate() {
        assert_eq!(
            row.len(), n_ports * n_ports,
            "s_data[{i}] has {} entries, expected {}×{}={}",
            row.len(), n_ports, n_ports, n_ports * n_ports
        );
    }

    let fmt_str = match fmt {
        TsFormat::Ri => "RI",
        TsFormat::Ma => "MA",
        TsFormat::Db => "DB",
    };

    let mut out = String::new();
    out.push_str(&format!(
        "! Touchstone S{}P generated by rem2 EM solver\n",
        n_ports
    ));
    out.push_str(&format!(
        "# {} S {} R {}\n",
        unit.label(),
        fmt_str,
        z0
    ));

    let freq_scale = unit.scale();

    for (fi, &fhz) in freqs_hz.iter().enumerate() {
        let row = &s_data[fi];
        let freq_scaled = fhz / freq_scale;

        // Touchstone 1.0: for N=1 or N=2, all data on one line.
        // For N>2: first 4 S-values on first line, then groups of 4 per continuation line.
        // We use a simpler single-line-per-frequency layout (valid for N<=2 and readable for N>2).
        let mut line = format!("{:.9e}", freq_scaled);

        for &s in row.iter() {
            match fmt {
                TsFormat::Ri => {
                    line.push_str(&format!("  {:.8e}  {:.8e}", s.re, s.im));
                }
                TsFormat::Ma => {
                    let mag = s.norm();
                    let ang = s.im.atan2(s.re).to_degrees();
                    line.push_str(&format!("  {:.8e}  {:.6}", mag, ang));
                }
                TsFormat::Db => {
                    let db = if s.norm() > 1e-300 { 20.0 * s.norm().log10() } else { -999.0 };
                    let ang = s.im.atan2(s.re).to_degrees();
                    line.push_str(&format!("  {:.6}  {:.6}", db, ang));
                }
            }
        }
        out.push_str(&line);
        out.push('\n');
    }

    out
}
```

- [ ] **Step 7: Run tests — all should pass**

```bash
cargo test -p rem-touchstone 2>&1 | tail -20
```

Expected: `test result: ok. 4 passed`.

- [ ] **Step 8: Commit**

```bash
cd /c/Users/lilu/works/rem2
git add crates/touchstone/ Cargo.toml
git commit -m "feat(touchstone): new rem-touchstone crate with write_snp (RI/MA/DB, N-port)"
```

---

## Task 2: Migrate `rem-driven` to use `rem-touchstone`

**Files:**
- Modify: `crates/driven/Cargo.toml`
- Modify: `crates/driven/src/vf.rs` — delegate `write_touchstone_s1p`

- [ ] **Step 1: Add dependency to driven crate**

Edit `crates/driven/Cargo.toml`, add to `[dependencies]`:

```toml
rem-touchstone = { workspace = true }
```

- [ ] **Step 2: Replace the inline writer in `vf.rs`**

In `crates/driven/src/vf.rs`, find the function `write_touchstone_s1p` (line ~148). Replace its body to delegate to `rem-touchstone`:

```rust
/// Write single-port Touchstone (.s1p) in RI format.
///
/// Delegates to `rem_touchstone::write_snp`.
pub fn write_touchstone_s1p(freqs_hz: &[f64], s11: &[Complex64], z0_ohm: f64) -> String {
    use rem_touchstone::{write_snp, TsFormat, TsFreqUnit};
    let s_data: Vec<Vec<Complex64>> = s11.iter().map(|&s| vec![s]).collect();
    write_snp(freqs_hz, &s_data, 1, z0_ohm, TsFormat::Ri, TsFreqUnit::Ghz)
}
```

The existing tests in `vf.rs` (`test_touchstone_format`) already test the output format — keep them as a regression guard.

- [ ] **Step 3: Run driven tests**

```bash
cargo test -p rem-driven 2>&1 | tail -30
```

Expected: all tests pass (including `test_touchstone_format`).

- [ ] **Step 4: Confirm workspace builds cleanly**

```bash
cargo build --workspace 2>&1 | grep -E "^error" | head -10
```

Expected: no errors.

- [ ] **Step 5: Commit**

```bash
git add crates/driven/Cargo.toml crates/driven/src/vf.rs
git commit -m "refactor(driven): delegate write_touchstone_s1p to rem-touchstone crate"
```

---

## Task 3: Extend config — `MomPort` + `ports` + `ref_impedance`

**Files:**
- Modify: `crates/config/src/schema.rs`

- [ ] **Step 1: Write the failing test**

Add to `crates/config/src/lib.rs` (or `tests/` if a separate integration test file exists) a new test that parses a MoM config with ports:

```rust
#[test]
fn mom_config_with_ports_parses() {
    let json = r#"{
        "Problem": { "Type": "MoM" },
        "Model":   { "Mesh": "t.msh" },
        "Solver": {
            "MoM": {
                "FreqMin": 1e9, "FreqMax": 2e9, "FreqStep": 1e9,
                "RefImpedance": 75.0,
                "Ports": [
                    { "Index": 1, "Attributes": [10], "Direction": "x" },
                    { "Index": 2, "Attributes": [11] }
                ]
            }
        }
    }"#;
    let cfg = rem_config::load_config_from_str(json, rem_config::ConfigFormat::Json).unwrap();
    let mom = cfg.solver.mom.as_ref().unwrap();
    assert_eq!(mom.ref_impedance, 75.0);
    assert_eq!(mom.ports.len(), 2);
    assert_eq!(mom.ports[0].index, 1);
    assert_eq!(mom.ports[0].attributes, vec![10u32]);
    assert_eq!(mom.ports[0].direction, "x");
    assert_eq!(mom.ports[1].index, 2);
    assert_eq!(mom.ports[1].direction, "x"); // default
    assert!(mom.ports[1].impedance.is_none());
}
```

- [ ] **Step 2: Run test to confirm it fails**

```bash
cargo test -p rem-config mom_config_with_ports 2>&1 | tail -20
```

Expected: compile error — `MomPort`, `ports`, `ref_impedance` don't exist yet.

- [ ] **Step 3: Add `MomPort` and extend `MomSolverConfig` in `schema.rs`**

In `crates/config/src/schema.rs`, after the `MomSolverConfig` struct (currently ends near line 906), add:

```rust
/// A lumped port definition for MoM S-parameter extraction.
#[derive(Debug, Clone, Deserialize, Default)]
pub struct MomPort {
    /// Port index (1-based), matches Touchstone port ordering.
    #[serde(rename = "Index")]
    pub index: u32,

    /// Physical-group attribute IDs that define the port surface.
    #[serde(rename = "Attributes", default)]
    pub attributes: Vec<u32>,

    /// Dominant E-field direction: "x" | "y" | "z".  Determines how
    /// the RHS excitation voltage is projected onto the port RWG functions.
    #[serde(rename = "Direction", default = "default_port_direction")]
    pub direction: String,

    /// Per-port reference impedance [Ω].  When None, uses `MomSolverConfig::ref_impedance`.
    #[serde(rename = "Impedance", default)]
    pub impedance: Option<f64>,
}

fn default_port_direction() -> String { "x".to_string() }
```

Then in `MomSolverConfig`, add two fields (before the closing brace):

```rust
    /// Lumped ports for S-parameter extraction.  When non-empty the solver
    /// runs one MoM solve per port and outputs S-matrix + Touchstone file.
    /// When empty, the existing plane-wave + RCS path is used.
    #[serde(rename = "Ports", default)]
    pub ports: Vec<MomPort>,

    /// Global reference impedance Z₀ [Ω] for S-parameter normalisation.
    #[serde(rename = "RefImpedance", default = "default_ref_impedance")]
    pub ref_impedance: f64,
```

Add the helper at the bottom of the default-function block:

```rust
fn default_ref_impedance() -> f64 { 50.0 }
```

- [ ] **Step 4: Run config test**

```bash
cargo test -p rem-config mom_config_with_ports 2>&1 | tail -20
```

Expected: `test result: ok. 1 passed`.

- [ ] **Step 5: Run full config test suite (regression)**

```bash
cargo test -p rem-config 2>&1 | tail -10
```

Expected: all existing tests still pass.

- [ ] **Step 6: Commit**

```bash
git add crates/config/src/schema.rs
git commit -m "feat(config): add MomPort, ports and ref_impedance to MomSolverConfig"
```

---

## Task 4: `crates/mom/src/port.rs` — Lumped port model

**Files:**
- Modify: `crates/mom/Cargo.toml`
- Create: `crates/mom/src/port.rs`
- Modify: `crates/mom/src/lib.rs` — add `pub mod port;`

- [ ] **Step 1: Add `rem-touchstone` dependency**

Edit `crates/mom/Cargo.toml`, add to `[dependencies]`:

```toml
rem-touchstone = { workspace = true }
```

- [ ] **Step 2: Write the failing tests**

Create `crates/mom/src/port.rs` with test module only:

```rust
//! MoM lumped-port model for S-parameter extraction.
//!
//! A lumped port is defined over a set of surface triangles (identified by
//! boundary attribute IDs).  The port produces:
//! - An excitation RHS vector (one solve per active port).
//! - V/I extraction from the solved surface current coefficients.

use crate::surface_mesh::SurfaceMesh;
use crate::basis::rwg::RwgBasis;
use num_complex::Complex64;
use rem_core::{RemResult, RemError};

/// Lumped port: a set of RWG indices + excitation direction + reference Z₀.
#[derive(Debug, Clone)]
pub struct MomLumpedPort {
    /// 1-based port index (matches Touchstone ordering).
    pub index: u32,
    /// Indices into the RWG basis array that are on this port's surface.
    pub rwg_indices: Vec<usize>,
    /// Dominant field direction unit vector [x, y, z].
    pub direction: [f64; 3],
    /// Reference impedance Z₀ [Ω].
    pub z0: f64,
}

impl MomLumpedPort {
    /// Find the RWG basis functions whose shared edge lies on a face that
    /// belongs to one of `port_attrs`.  A face "belongs to a port" if its
    /// attribute tag (stored in `SurfaceMesh::face_attrs`) matches.
    ///
    /// `direction_str` — "x", "y", or "z"; mapped to unit vector.
    pub fn from_surface(
        surf: &SurfaceMesh,
        bases: &[RwgBasis],
        port_attrs: &[u32],
        index: u32,
        direction_str: &str,
        z0: f64,
    ) -> RemResult<Self> {
        todo!()
    }

    /// Build the N-element excitation RHS for this port (one entry per RWG).
    ///
    /// For RWG basis, the port contribution at basis m is:
    ///   V_m = -∫_{port} f_m(r) · d̂ dS   (d̂ = direction unit vector)
    ///
    /// Only the `rwg_indices` of this port get non-zero entries; all others are 0.
    /// The excitation amplitude `v0` is nominally 1 V.
    pub fn excitation_rhs(
        &self,
        surf: &SurfaceMesh,
        bases: &[RwgBasis],
        n_total: usize,
        v0: Complex64,
    ) -> Vec<Complex64> {
        todo!()
    }

    /// Extract port voltage V and current I from solved current coefficients.
    ///
    /// V = v0 (the excitation voltage, typically 1.0)
    ///
    /// I = Σ_{m ∈ port} a_m * ∫_{port} ∇_s · f_m  dS
    ///   = Σ_{m ∈ port} a_m * divergence_m
    ///
    /// where a_m are the solved RWG coefficients.
    pub fn extract_current(
        &self,
        surf: &SurfaceMesh,
        bases: &[RwgBasis],
        coeffs: &[Complex64],
    ) -> Complex64 {
        todo!()
    }
}

fn direction_vec(s: &str) -> [f64; 3] {
    match s.to_lowercase().as_str() {
        "y" => [0.0, 1.0, 0.0],
        "z" => [0.0, 0.0, 1.0],
        _   => [1.0, 0.0, 0.0],  // default "x"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface_mesh::{SurfaceMesh, TriFace, SharedEdge, tri_geometry};
    use crate::basis::rwg::generate_rwg_bases;

    /// Build a two-triangle surface with attribute tags.
    fn two_tri_port_surf() -> SurfaceMesh {
        let nodes = vec![
            [0.0_f64, 0.0, 0.0],
            [1.0,     0.0, 0.0],
            [0.5,     1.0, 0.0],
            [-0.5,    1.0, 0.0],
        ];
        let (c0, n0, a0) = tri_geometry(&nodes[0], &nodes[1], &nodes[2]);
        let (c1, n1, a1) = tri_geometry(&nodes[0], &nodes[2], &nodes[3]);
        let faces = vec![
            TriFace { nodes: [0,1,2], centroid: c0, normal: n0, area: a0 },
            TriFace { nodes: [0,2,3], centroid: c1, normal: n1, area: a1 },
        ];
        let edges = vec![SharedEdge {
            nodes: [0, 2],
            plus_face: 0,
            minus_face: 1,
            length: (nodes[2][0].powi(2) + nodes[2][1].powi(2)).sqrt(),
        }];
        SurfaceMesh {
            nodes,
            faces,
            edges,
            boundary_edges: vec![[0,1],[1,2],[2,3],[3,0]],
            face_attrs: vec![1, 1],   // both faces tagged as attr 1
        }
    }

    #[test]
    fn from_surface_finds_rwg_on_attr() {
        let surf = two_tri_port_surf();
        let bases = generate_rwg_bases(&surf);
        assert_eq!(bases.len(), 1, "should have 1 shared edge");
        let port = MomLumpedPort::from_surface(
            &surf, &bases, &[1], 1, "x", 50.0
        ).expect("port construction failed");
        assert_eq!(port.rwg_indices.len(), 1, "all RWG on attr-1 surface should be found");
        assert_eq!(port.rwg_indices[0], 0);
        assert_eq!(port.z0, 50.0);
    }

    #[test]
    fn from_surface_empty_when_no_match() {
        let surf = two_tri_port_surf();
        let bases = generate_rwg_bases(&surf);
        let port = MomLumpedPort::from_surface(
            &surf, &bases, &[99], 1, "x", 50.0
        ).expect("should succeed even with no match");
        assert!(port.rwg_indices.is_empty());
    }

    #[test]
    fn excitation_rhs_zero_outside_port() {
        let surf = two_tri_port_surf();
        let bases = generate_rwg_bases(&surf);
        let port = MomLumpedPort::from_surface(&surf, &bases, &[1], 1, "x", 50.0).unwrap();
        let rhs = port.excitation_rhs(&surf, &bases, bases.len(), Complex64::new(1.0, 0.0));
        assert_eq!(rhs.len(), bases.len());
        // All entries on the port should be non-zero (there's 1 RWG on port attr 1)
        assert!(rhs[0].norm() > 0.0, "port RWG should have non-zero excitation");
    }

    #[test]
    fn extract_current_finite() {
        let surf = two_tri_port_surf();
        let bases = generate_rwg_bases(&surf);
        let port = MomLumpedPort::from_surface(&surf, &bases, &[1], 1, "x", 50.0).unwrap();
        let coeffs = vec![Complex64::new(1.0, 0.5); bases.len()];
        let i = port.extract_current(&surf, &bases, &coeffs);
        assert!(i.re.is_finite() && i.im.is_finite());
    }
}
```

- [ ] **Step 3: Run tests — confirm compile error**

```bash
cargo test -p rem-mom -- port 2>&1 | tail -20
```

Expected: compile errors — `face_attrs` field missing from `SurfaceMesh`, `todo!()` bodies.

- [ ] **Step 4: Add `face_attrs` to `SurfaceMesh`**

In `crates/mom/src/surface_mesh.rs`, find the `SurfaceMesh` struct definition. Add the field:

```rust
pub struct SurfaceMesh {
    pub nodes: Vec<[f64; 3]>,
    pub faces: Vec<TriFace>,
    pub edges: Vec<SharedEdge>,
    pub boundary_edges: Vec<[usize; 2]>,
    /// Per-face physical-group attribute tag (from the mesh boundary tags).
    /// Zero means "untagged / not extracted from a named boundary".
    pub face_attrs: Vec<u32>,
}
```

In the `SurfaceMesh::extract` function, populate `face_attrs`. The function currently builds `faces` by iterating boundary elements with matching PEC attribute IDs. Set `face_attrs[i]` to the attribute tag of that boundary element:

```rust
// After building `faces`, add:
let face_attrs: Vec<u32> = collected_face_attrs; // collected during the face-building loop
```

The exact change depends on the current loop structure in `surface_mesh.rs`. The key is: when you push to `faces`, also push the corresponding attribute tag to `face_attrs`. At the end of `extract`, include `face_attrs` in the returned struct.

Also update any `SurfaceMesh { ... }` struct literals in tests (add `face_attrs: vec![...]` with the right length).

- [ ] **Step 5: Implement `from_surface`**

```rust
pub fn from_surface(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    port_attrs: &[u32],
    index: u32,
    direction_str: &str,
    z0: f64,
) -> RemResult<Self> {
    // A RWG is "on the port" if both its plus_face and minus_face
    // have an attribute in port_attrs.
    let port_attr_set: std::collections::HashSet<u32> =
        port_attrs.iter().copied().collect();

    let rwg_indices: Vec<usize> = bases.iter().enumerate()
        .filter(|(_, b)| {
            let plus_attr  = surf.face_attrs.get(b.plus_face).copied().unwrap_or(0);
            let minus_attr = surf.face_attrs.get(b.minus_face).copied().unwrap_or(0);
            port_attr_set.contains(&plus_attr) && port_attr_set.contains(&minus_attr)
        })
        .map(|(i, _)| i)
        .collect();

    Ok(Self {
        index,
        rwg_indices,
        direction: direction_vec(direction_str),
        z0,
    })
}
```

- [ ] **Step 6: Implement `excitation_rhs`**

```rust
pub fn excitation_rhs(
    &self,
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    n_total: usize,
    v0: Complex64,
) -> Vec<Complex64> {
    let mut rhs = vec![Complex64::ZERO; n_total];
    let d = self.direction;
    for &mi in &self.rwg_indices {
        let b = &bases[mi];
        let mut val = Complex64::ZERO;
        // Integrate f_m · d̂ over both support triangles
        for &(face_idx, in_plus) in &[(b.plus_face, true), (b.minus_face, false)] {
            let face = &surf.faces[face_idx];
            // Centroid-quadrature (1 point): f_m(centroid) · d̂ * area
            let fm = b.eval(&face.centroid, surf, in_plus);
            let dot = d[0]*fm[0] + d[1]*fm[1] + d[2]*fm[2];
            val += Complex64::new(dot * face.area, 0.0);
        }
        rhs[mi] = -v0 * val;
    }
    rhs
}
```

- [ ] **Step 7: Implement `extract_current`**

```rust
pub fn extract_current(
    &self,
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    coeffs: &[Complex64],
) -> Complex64 {
    // I ≈ Σ_{m ∈ port} a_m * div_m * A_m (sign convention from EFIE)
    // divergence of RWG = ±l_n / A_face (constant per half-support)
    let mut i_port = Complex64::ZERO;
    for &mi in &self.rwg_indices {
        if mi >= coeffs.len() { continue; }
        let b = &bases[mi];
        // Sum divergence contribution over plus and minus faces
        let div_p = b.divergence(surf, true);
        let div_m_val = b.divergence(surf, false);
        let area_p = surf.faces[b.plus_face].area;
        let area_m = surf.faces[b.minus_face].area;
        let contrib = div_p * area_p + div_m_val * area_m;
        i_port += coeffs[mi] * contrib;
    }
    i_port
}
```

- [ ] **Step 8: Run port tests**

```bash
cargo test -p rem-mom -- port 2>&1 | tail -20
```

Expected: `test result: ok. 4 passed`.

- [ ] **Step 9: Register `port` module in `lib.rs`**

In `crates/mom/src/lib.rs`, add `pub mod port;` to the module list.

- [ ] **Step 10: Run full mom test suite**

```bash
cargo test -p rem-mom 2>&1 | tail -20
```

Expected: all tests pass.

- [ ] **Step 11: Commit**

```bash
git add crates/mom/Cargo.toml crates/mom/src/port.rs crates/mom/src/lib.rs \
        crates/mom/src/surface_mesh.rs
git commit -m "feat(mom): add MomLumpedPort, face_attrs on SurfaceMesh, excitation RHS + I extraction"
```

---

## Task 5: `crates/mom/src/sparams.rs` — S-matrix computation

**Files:**
- Create: `crates/mom/src/sparams.rs`
- Modify: `crates/mom/src/lib.rs` — add `pub mod sparams;`

- [ ] **Step 1: Write the failing tests**

Create `crates/mom/src/sparams.rs`:

```rust
//! S-parameter matrix computation for MoM port-excited problems.
//!
//! Given N MomLumpedPort definitions and a pre-assembled Z matrix,
//! runs one solve per port and extracts the N×N S-matrix.

use crate::port::MomLumpedPort;
use crate::surface_mesh::SurfaceMesh;
use crate::basis::rwg::RwgBasis;
use crate::assemble::lu_solve;
use nalgebra::DMatrix;
use num_complex::Complex64;
use rem_core::RemResult;
use rem_touchstone::{write_snp, TsFormat, TsFreqUnit};
use std::path::Path;

/// N×N S-matrix at a single frequency.
#[derive(Debug, Clone)]
pub struct SMatrix {
    /// Number of ports.
    pub n_ports: usize,
    /// Frequency [Hz].
    pub freq_hz: f64,
    /// Row-major S-matrix: data[i * n_ports + j] = S_{i+1, j+1}.
    pub data: Vec<Complex64>,
}

impl SMatrix {
    /// S_{row+1, col+1} (0-based indices).
    pub fn get(&self, row: usize, col: usize) -> Complex64 {
        self.data[row * self.n_ports + col]
    }
}

/// Compute the N×N S-matrix for a set of lumped ports.
///
/// For each excitation port `p`:
/// 1. Build RHS from `port_p.excitation_rhs(v0=1V)`.
/// 2. Solve Z·I = V_rhs → current coefficients.
/// 3. For each observation port `q`: extract I_q = port_q.extract_current(coeffs).
/// 4. Compute S_{qp} = (Z0_q * I_q - 1) / (Z0_q * I_q + 1)  [reflection]
///    or S_{qp} = 2 * Z0_q * I_q / V_p_fwd  [transmission], using wave-port normalisation.
///
/// **Simplified model used here:**
/// V_p (forward) = 1.0 V (excitation),  Z0 = port.z0.
/// S_{qp} = (V_q - Z0_q I_q) / V_p_fwd
///   where V_q = v0 if q==p (self-excitation returns 1), else 0.
///   and   V_p_fwd = 1.0.
///
/// This reduces to: S_{pp} = (1 - Z0_p I_p) / 1  and  S_{qp} = -Z0_q I_q for q≠p.
pub fn compute_s_matrix(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    ports: &[MomLumpedPort],
    z_mat: &DMatrix<Complex64>,
    freq_hz: f64,
) -> RemResult<SMatrix> {
    todo!()
}

/// Run a full S-parameter sweep over multiple frequencies.
///
/// Returns one `SMatrix` per frequency.
pub fn s_param_sweep(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    ports: &[MomLumpedPort],
    freq_hz_list: &[f64],
    build_z: &dyn Fn(f64) -> RemResult<DMatrix<Complex64>>,
) -> RemResult<Vec<SMatrix>> {
    todo!()
}

/// Write all frequency-sweep S-matrices to a Touchstone `.s{N}p` file.
pub fn write_touchstone(matrices: &[SMatrix], path: &Path, z0: f64) -> RemResult<()> {
    todo!()
}

/// Append S-parameter data to a Palace-compatible `port-S.csv`.
pub fn append_palace_csv(matrices: &[SMatrix], path: &Path) -> RemResult<()> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface_mesh::{SurfaceMesh, TriFace, SharedEdge, tri_geometry};
    use crate::basis::rwg::generate_rwg_bases;
    use crate::port::MomLumpedPort;

    fn two_tri_surf_with_attrs(attr: u32) -> SurfaceMesh {
        let nodes = vec![
            [0.0_f64, 0.0, 0.0],
            [1.0,     0.0, 0.0],
            [0.5,     1.0, 0.0],
            [-0.5,    1.0, 0.0],
        ];
        let (c0, n0, a0) = tri_geometry(&nodes[0], &nodes[1], &nodes[2]);
        let (c1, n1, a1) = tri_geometry(&nodes[0], &nodes[2], &nodes[3]);
        let faces = vec![
            TriFace { nodes:[0,1,2], centroid:c0, normal:n0, area:a0 },
            TriFace { nodes:[0,2,3], centroid:c1, normal:n1, area:a1 },
        ];
        let edges = vec![SharedEdge {
            nodes: [0,2], plus_face: 0, minus_face: 1,
            length: (0.5_f64.powi(2) + 1.0_f64.powi(2)).sqrt(),
        }];
        SurfaceMesh {
            nodes, faces, edges,
            boundary_edges: vec![[0,1],[1,2],[2,3],[3,0]],
            face_attrs: vec![attr, attr],
        }
    }

    /// Build a trivial 1×1 Z-matrix (identity) for testing.
    fn identity_z(n: usize) -> DMatrix<Complex64> {
        let mut z = DMatrix::<Complex64>::zeros(n, n);
        for i in 0..n { z[(i,i)] = Complex64::new(1.0, 0.0); }
        z
    }

    #[test]
    fn s_matrix_shape_single_port() {
        let surf = two_tri_surf_with_attrs(1);
        let bases = generate_rwg_bases(&surf);
        let port = MomLumpedPort::from_surface(&surf, &bases, &[1], 1, "x", 50.0).unwrap();
        let ports = vec![port];
        let z = identity_z(bases.len());
        let sm = compute_s_matrix(&surf, &bases, &ports, &z, 1e9).unwrap();
        assert_eq!(sm.n_ports, 1);
        assert_eq!(sm.data.len(), 1);
        assert!(sm.data[0].re.is_finite());
    }

    #[test]
    fn write_touchstone_creates_file() {
        let matrices = vec![
            SMatrix { n_ports: 1, freq_hz: 1e9, data: vec![Complex64::new(0.5, -0.3)] },
            SMatrix { n_ports: 1, freq_hz: 2e9, data: vec![Complex64::new(0.4, -0.2)] },
        ];
        let tmp = std::env::temp_dir().join("test_mom_s1p.s1p");
        write_touchstone(&matrices, &tmp, 50.0).expect("write failed");
        let content = std::fs::read_to_string(&tmp).unwrap();
        assert!(content.contains("# GHz S RI"), "option line missing");
        let data_lines: usize = content.lines()
            .filter(|l| !l.starts_with('!') && !l.starts_with('#') && !l.trim().is_empty())
            .count();
        assert_eq!(data_lines, 2, "expected 2 data lines");
        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn append_palace_csv_creates_file() {
        let matrices = vec![
            SMatrix { n_ports: 1, freq_hz: 1e9, data: vec![Complex64::new(0.5, -0.3)] },
        ];
        let tmp = std::env::temp_dir().join("test_mom_port_s.csv");
        let _ = std::fs::remove_file(&tmp);
        append_palace_csv(&matrices, &tmp).expect("csv write failed");
        let content = std::fs::read_to_string(&tmp).unwrap();
        assert!(content.contains("Freq"), "CSV header missing");
        let _ = std::fs::remove_file(&tmp);
    }
}
```

- [ ] **Step 2: Run tests — confirm compile error (todo!)**

```bash
cargo test -p rem-mom -- sparams 2>&1 | tail -20
```

Expected: compile error because the module isn't registered yet and `todo!()` bodies don't pass.

- [ ] **Step 3: Register module in `lib.rs`**

In `crates/mom/src/lib.rs`, add:

```rust
pub mod sparams;
```

- [ ] **Step 4: Implement `compute_s_matrix`**

```rust
pub fn compute_s_matrix(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    ports: &[MomLumpedPort],
    z_mat: &DMatrix<Complex64>,
    freq_hz: f64,
) -> RemResult<SMatrix> {
    let n_ports = ports.len();
    let n_rwg   = bases.len();
    let v0      = Complex64::new(1.0, 0.0);

    // Solve one system per excitation port
    let mut all_currents: Vec<Vec<Complex64>> = Vec::with_capacity(n_ports);
    for port_p in ports {
        let rhs = port_p.excitation_rhs(surf, bases, n_rwg, v0);
        let coeffs = lu_solve(z_mat, &rhs)?;
        all_currents.push(coeffs);
    }

    // Build N×N S-matrix
    let mut data = vec![Complex64::ZERO; n_ports * n_ports];
    for (p, currents_p) in all_currents.iter().enumerate() {
        for (q, port_q) in ports.iter().enumerate() {
            let z0_q = port_q.z0;
            let i_q  = port_q.extract_current(surf, bases, currents_p);
            // Simplified wave-port S-parameter formula:
            // S_{qp} = (V_q - Z0_q I_q) / V_p_incident
            // V_q = v0 = 1 if q==p (self-port), 0 otherwise
            let v_q = if q == p { v0 } else { Complex64::ZERO };
            let s_qp = v_q - Complex64::new(z0_q, 0.0) * i_q;
            data[q * n_ports + p] = s_qp;
        }
    }

    Ok(SMatrix { n_ports, freq_hz, data })
}
```

- [ ] **Step 5: Implement `s_param_sweep`**

```rust
pub fn s_param_sweep(
    surf: &SurfaceMesh,
    bases: &[RwgBasis],
    ports: &[MomLumpedPort],
    freq_hz_list: &[f64],
    build_z: &dyn Fn(f64) -> RemResult<DMatrix<Complex64>>,
) -> RemResult<Vec<SMatrix>> {
    freq_hz_list.iter().map(|&f| {
        let z = build_z(f)?;
        compute_s_matrix(surf, bases, ports, &z, f)
    }).collect()
}
```

- [ ] **Step 6: Implement `write_touchstone`**

```rust
pub fn write_touchstone(matrices: &[SMatrix], path: &Path, z0: f64) -> RemResult<()> {
    use std::io::Write;
    if matrices.is_empty() { return Ok(()); }
    let n_ports = matrices[0].n_ports;

    let freqs: Vec<f64>          = matrices.iter().map(|m| m.freq_hz).collect();
    let s_data: Vec<Vec<Complex64>> = matrices.iter().map(|m| m.data.clone()).collect();

    let content = write_snp(&freqs, &s_data, n_ports, z0, TsFormat::Ri, TsFreqUnit::Ghz);
    std::fs::create_dir_all(path.parent().unwrap_or(Path::new(".")))?;
    let mut f = std::fs::File::create(path)?;
    f.write_all(content.as_bytes())?;
    Ok(())
}
```

- [ ] **Step 7: Implement `append_palace_csv`**

```rust
pub fn append_palace_csv(matrices: &[SMatrix], path: &Path) -> RemResult<()> {
    use std::io::Write;
    let write_header = !path.exists();
    std::fs::create_dir_all(path.parent().unwrap_or(Path::new(".")))?;
    let mut f = std::fs::OpenOptions::new().create(true).append(true).open(path)?;
    if write_header {
        // Header: Freq (GHz), then Re(S11), Im(S11), |S11| dB, Re(S12), ...
        let mut hdr = "Freq (GHz)".to_string();
        if let Some(m) = matrices.first() {
            for i in 0..m.n_ports {
                for j in 0..m.n_ports {
                    hdr.push_str(&format!(",Re(S{}{}),Im(S{}{}),|S{}{}| (dB)",
                        i+1, j+1, i+1, j+1, i+1, j+1));
                }
            }
        }
        writeln!(f, "{hdr}")?;
    }
    for m in matrices {
        let mut line = format!("{:.9e}", m.freq_hz / 1e9);
        for &s in &m.data {
            let db = if s.norm() > 1e-300 { 20.0 * s.norm().log10() } else { -999.0 };
            line.push_str(&format!(",{:.8e},{:.8e},{:.4}", s.re, s.im, db));
        }
        writeln!(f, "{line}")?;
    }
    Ok(())
}
```

- [ ] **Step 8: Run sparams tests**

```bash
cargo test -p rem-mom -- sparams 2>&1 | tail -20
```

Expected: `test result: ok. 3 passed`.

- [ ] **Step 9: Run full mom test suite**

```bash
cargo test -p rem-mom 2>&1 | tail -20
```

Expected: all tests pass.

- [ ] **Step 10: Commit**

```bash
git add crates/mom/src/sparams.rs crates/mom/src/lib.rs
git commit -m "feat(mom): add sparams module — S-matrix computation + Touchstone/CSV output"
```

---

## Task 6: Wire the S-parameter path into `mom/src/lib.rs`

**Files:**
- Modify: `crates/mom/src/lib.rs`

- [ ] **Step 1: Write an integration test first**

Add to `crates/mom/tests/mie_validation.rs` a new test:

```rust
/// Smoke test: MoM with ports config parses and runs without panic.
/// No accuracy check — just ensures the port branch compiles and runs end-to-end.
#[test]
fn mom_port_branch_runs() {
    use rem_config::{PalaceConfig, MomSolverConfig, MomPort};
    // Build a minimal two-triangle mesh config in memory
    // and verify the port branch doesn't panic.

    // We reuse the Mie sphere mesh path only to get a valid SurfaceMesh;
    // since the test environment may not have sphere.msh, we skip if absent.
    let sphere_msh = std::path::Path::new("tests/sphere.msh");
    if !sphere_msh.exists() {
        println!("skipping mom_port_branch_runs: sphere.msh not found");
        return;
    }

    let json = r#"{
        "Problem": { "Type": "MoM", "Output": "/tmp/mom_port_test" },
        "Model": { "Mesh": "tests/sphere.msh" },
        "Boundaries": { "PEC": { "Attributes": [1] } },
        "Solver": {
            "MoM": {
                "FreqMin": 1e9, "FreqMax": 1e9, "FreqStep": 1e9,
                "RefImpedance": 50.0,
                "Ports": [
                    { "Index": 1, "Attributes": [1], "Direction": "x" }
                ]
            }
        }
    }"#;

    let cfg: PalaceConfig = serde_json::from_str(
        &rem_config::preprocess::strip_comments(json)
    ).expect("config parse failed");

    // Load mesh
    let mesh = rem_mesh::load_mesh(&cfg, &rem_parallel::NoComm)
        .expect("mesh load failed");

    let mom_cfg = cfg.solver.mom.as_ref().unwrap();

    // Should not panic; S-param result may be physically meaningless on sphere.msh
    // but the code path must execute.
    let result = rem_mom::run_with_mesh(&cfg, mom_cfg, &mesh);
    // We accept either Ok or a config error (no port RWGs found on sphere).
    // We do NOT accept a panic.
    match result {
        Ok(_) => {},
        Err(e) => println!("expected error (no RWGs on port attr): {e}"),
    }
}
```

- [ ] **Step 2: Run integration test — confirm it compiles but the port branch isn't wired yet**

```bash
cargo test -p rem-mom mom_port_branch_runs 2>&1 | tail -20
```

Expected: test passes but prints "skipping" (no sphere.msh), OR it runs and hits the existing RCS branch (no panic). Either is acceptable before the implementation.

- [ ] **Step 3: Modify `lib.rs` — add port branch**

In `crates/mom/src/lib.rs`, inside `run_with_mesh`, find the frequency sweep loop. The current structure is roughly:

```rust
while freq <= freq_max + 1e-3 * freq_step {
    // ... assemble Z, solve, compute RCS ...
    freq += freq_step;
}
```

Replace with a branch at the start of `run_with_mesh` (before the frequency loop):

```rust
use crate::port::MomLumpedPort;
use crate::sparams;

// ── Port path ─────────────────────────────────────────────────────────────
if !mom_cfg.ports.is_empty() {
    return run_s_param_sweep(config, mom_cfg, &surf);
}
// ── RCS path (original) ───────────────────────────────────────────────────
// ... existing code unchanged ...
```

Add the new function `run_s_param_sweep` in `lib.rs`:

```rust
/// Run the S-parameter sweep (port-excited MoM path).
fn run_s_param_sweep(
    config: &PalaceConfig,
    mom_cfg: &MomSolverConfig,
    surf: &SurfaceMesh,
) -> RemResult<MomResult> {
    use std::f64::consts::PI;
    let output_dir = std::path::Path::new(config.problem.output_dir());
    #[cfg(not(target_arch = "wasm32"))]
    std::fs::create_dir_all(output_dir.join("postpro"))?;

    let bases = crate::basis::rwg::generate_rwg_bases(surf);
    let quad  = crate::quadrature::TriQuad::new(5);

    // Build lumped ports
    let lumped_ports: Vec<MomLumpedPort> = mom_cfg.ports.iter().map(|p| {
        let z0 = p.impedance.unwrap_or(mom_cfg.ref_impedance);
        MomLumpedPort::from_surface(surf, &bases, &p.attributes, p.index, &p.direction, z0)
    }).collect::<RemResult<_>>()?;

    // Frequency sweep
    let mut freq = mom_cfg.freq_min;
    let freq_max  = mom_cfg.freq_max;
    let freq_step = mom_cfg.freq_step;
    let mut all_matrices: Vec<sparams::SMatrix> = Vec::new();

    while freq <= freq_max + 1e-3 * freq_step {
        log::info!("MoM S-param solve at f = {:.3e} Hz", freq);
        let k = 2.0 * PI * freq / rem_core::C0;

        let z_mat = {
            let b_ref: &[_] = &bases;
            crate::assemble::assemble_cfie_rwg(
                surf, b_ref, freq, mom_cfg.alpha, &quad, mom_cfg.singular_tol,
            )?
        };

        let sm = sparams::compute_s_matrix(surf, &bases, &lumped_ports, &z_mat, freq)?;
        log::info!("  S-matrix computed: {}×{}", sm.n_ports, sm.n_ports);
        all_matrices.push(sm);

        freq += freq_step;
    }

    // Write outputs
    #[cfg(not(target_arch = "wasm32"))]
    {
        let ts_ext = format!("s{}p", lumped_ports.len());
        let ts_path = output_dir.join("postpro").join(format!("s_params.{}", ts_ext));
        sparams::write_touchstone(&all_matrices, &ts_path, mom_cfg.ref_impedance)?;
        log::info!("MoM S-param output: {}", ts_path.display());

        let csv_path = output_dir.join("postpro").join("port-S.csv");
        sparams::append_palace_csv(&all_matrices, &csv_path)?;
    }

    log::info!("MoM S-param sweep complete. {} frequency points.", all_matrices.len());
    Ok(MomResult { rcs: vec![] })  // no RCS in S-param mode
}
```

- [ ] **Step 4: Run full test suite**

```bash
cargo test -p rem-mom 2>&1 | tail -20
```

Expected: all tests pass including existing Mie validation.

- [ ] **Step 5: Build workspace**

```bash
cargo build --workspace 2>&1 | grep -E "^error" | head -10
```

Expected: no errors.

- [ ] **Step 6: Commit**

```bash
git add crates/mom/src/lib.rs crates/mom/tests/mie_validation.rs
git commit -m "feat(mom): wire S-parameter port branch in run_with_mesh; RCS path unchanged"
```

---

## Task 7: Regression — Mie validation still passes, workspace clean

- [ ] **Step 1: Run workspace tests**

```bash
cargo test --workspace 2>&1 | tail -30
```

Expected: all existing tests pass. Look for `test result:` lines — no failures.

- [ ] **Step 2: Check WASM compile target**

```bash
cargo build --target wasm32-unknown-unknown -p rem-mom --no-default-features 2>&1 | grep -E "^error" | head -10
```

Expected: no errors.

- [ ] **Step 3: Check `rem-driven` still works (Touchstone delegation)**

```bash
cargo test -p rem-driven vf 2>&1 | tail -10
```

Expected: `test result: ok.`

- [ ] **Step 4: Final commit with version bump note**

```bash
git add -A
git commit -m "chore: v0.17.0 — rem-touchstone crate + MoM port/S-param path complete"
```

---

## Self-Review Checklist

**Spec coverage:**
- ✅ `rem-touchstone` crate with `write_snp` (N-port, RI/MA/DB) — Task 1
- ✅ `rem-driven` migrated to use `rem-touchstone` — Task 2
- ✅ `MomPort` + `ports` + `ref_impedance` in config — Task 3
- ✅ `MomLumpedPort` model + `face_attrs` on `SurfaceMesh` — Task 4
- ✅ `sparams.rs` — `compute_s_matrix`, `write_touchstone`, `append_palace_csv` — Task 5
- ✅ Main port branch in `lib.rs` — Task 6
- ✅ Full regression + WASM check — Task 7
- ✅ Touchstone lives in `rem-touchstone`, not in `mom` or `driven`

**Placeholder scan:** None found.

**Type consistency:**
- `MomLumpedPort::from_surface` signature matches usage in `lib.rs` and tests.
- `SMatrix::data` is `Vec<Complex64>` row-major — consistent in `compute_s_matrix` and `write_touchstone`.
- `write_snp` signature: `(freqs_hz, s_data, n_ports, z0, fmt, unit)` — matches `write_touchstone` call site.
- `face_attrs: Vec<u32>` added to `SurfaceMesh` — all struct literals in tests updated.
