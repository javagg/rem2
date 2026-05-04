# Palace Example Synchronization

> This document records the process for keeping REM's Palace-compatible examples
> in sync with the upstream Palace project. Last sync: 2026-04-08.

## Palace Repository

```bash
# Clone Palace (one-time)
git clone --depth=1 https://github.com/awslabs/palace.git /tmp/palace-src

# Sync examples directory
cd /tmp/palace-src && find examples -name "*.json" | sort
```

## Current REM Examples (complete Palace coverage)

All Palace JSON files are present in `examples/palace/` with aligned directory
structure and filenames, adapted for REM's config schema.

| Subdirectory (`examples/palace/`) | Files | Solver Type | REM Status |
|-------------|-------|-------------|-----------|
| `adapter/` | `hybrid.json` | Driven | ✅ |
| `antenna/` | `antenna_halfwave_dipole.json`, `antenna_short_dipole.json` | Driven | ✅ |
| `coaxial/` | `coaxial_matched.json`, `coaxial_open.json`, `coaxial_short.json` | Transient | ✅ |
| `cpw/` | 8 files | Driven/Eigenmode | ✅ |
| `cylinder/` | `cavity_impedance.json`, `cavity_pec.json`, `driven_wave.json`, `floquet.json`, `waveguide.json` | Eigenmode/Driven | ✅ |
| `rings/` | `rings.json` | Magnetostatic | ✅ |
| `spheres/` | `spheres.json` | Electrostatic | ✅ |
| `transmon/` | `transmon_amr.json`, `transmon_coarse.json` | Eigenmode | ✅ |

## REM Quick Examples (web/wasm)

Quick-check examples are under `examples/rem/` and wired to the Yew demo.

| Yew Key | Config Path | Problem Type | Mesh Source |
|---------|-------------|--------------|-------------|
| `rem_es_fast` | `examples/rem/es_parallel_plate_fast/es_parallel_plate_fast.json` | Electrostatic | `examples/rem/parallel_plate/mesh/plate_2d.msh` |
| `rem_ms_fast` | `examples/rem/ms_parallel_plate_fast/ms_parallel_plate_fast.json` | Magnetostatic | `examples/rem/parallel_plate/mesh/plate_2d.msh` |
| `rem_driven_fast` | `examples/rem/driven_cpw_fast/driven_cpw_fast.json` | Driven | `examples/palace/cpw/mesh/cpw_coax.msh` |
| `rem_eigen_fast` | `examples/rem/eigen_cylinder_fast/eigen_cylinder_fast.json` | Eigenmode | `examples/palace/cylinder/mesh/cylinder_tet.msh` |
| `rem_transient_fast` | `examples/rem/transient_coax_fast/transient_coax_fast.json` | Transient | `examples/palace/coaxial/mesh/coaxial.msh` |
| `rem_mom_fast` | `examples/rem/mom_sphere_fast/mom_sphere_fast.json` | MoM | `examples/rem/sbr_sphere/mesh/sphere.msh` |
| `rem_sbr_fast` | `examples/rem/sbr_sphere_fast/sbr_sphere_fast.json` | SBR | `examples/rem/sbr_sphere/mesh/sphere.msh` |

## Adaptation Rules (Palace → REM)

When syncing a Palace JSON to REM, the following fields require adaptation:

### Supported fields (keep as-is)
- `Problem.Type` — "Electrostatic", "Magnetostatic", "Eigenmode", "Driven", "Transient"
- `Problem.Output` — convert `postpro/xxx` → `output/xxx`
- `Model.Mesh` — keep path (e.g. `mesh/cylinder_hex.msh`); ensure mesh file exists
- `Model.L0` — keep (mesh scale factor)
- `Model.Refinement.MaxIter` → `Model.Refinement.MaxIter` (REM uses `MaxIter`)
- `Domains.Materials[*].Permittivity` — scalar only (REM; use first component if Palace uses array)
- `Domains.Materials[*].Permeability` — scalar only (REM; use first component if Palace uses array)
- `Domains.Materials[*].LossTan` — scalar; REM also accepts `LossTan: [x, y, z]` array form (uses first element)
- `Domains.Materials[*].Conductivity` — scalar only
- `Boundaries.PEC`, `Boundaries.PMC`, `Boundaries.Ground`
- `Boundaries.LumpedPort` — fields: `Index`, `Attributes`, `Direction`, `R`, `L`, `C`, `Excitation`
- `Boundaries.WavePort` — fields: `Index`, `Attributes`, `Excitation`, `Mode`
- `Boundaries.Absorbing` — fields: `Attributes`, `Order`
- `Boundaries.Impedance` — fields: `Attributes`, `Rs`
- `Solver.Order`, `Solver.Eigenmode`, `Solver.Driven`, `Solver.Transient`
- `Solver.Linear.Tol`, `Solver.Linear.MaxIter`

### Accepted with warnings (not yet implemented in REM)
These fields are deserialized and logged as warnings; they do NOT need to be removed from Palace JSONs.
- `Problem.Verbose` — accepted, value is ignored
- `Problem.OutputFormats.GridFunction` — VTK solution written to `<output>/paraview/solution.vtk` by all solvers; field is accepted and logs info
- `Solver.Device: "CPU"` — REM is CPU-only; ignored with warning
- `Solver.Linear.KSPType` — "GMRES" (default) and "CG"/"PCG" (routes to PCG for SPD/Helmholtz solves) are supported; other values are logged and ignored
- `Solver.Linear.MGLevels` — algebraic multigrid not implemented; ignored
- `Solver.Linear.ComplexCoarseSolve` — complex coarse-grid solve not implemented; ignored
- `Domains.Materials[*].LossTan` (array form) — anisotropic loss not implemented; uses first element
- `Domains.Postprocessing.Energy` — per-group energy written to `postpro/energy-E.csv` (Electrostatic) and `postpro/energy-B.csv` (Magnetostatic); groups sum energy over the specified domain attributes
- `Domains.Postprocessing.Probe` — field probe sampling is implemented for Electrostatic (φ + E-field) and Magnetostatic (A_z); Eigenmode writes all modes to `postpro/probe-phi-modes.csv`
- `Boundaries.Periodic` — Γ-point periodic BCs supported; complex Floquet (non-zero FloquetWaveVector) logs warning and skips
- `Boundaries.Postprocessing` — `Electric` (displacement flux → `postpro/surface-flux.csv` in C) and `Magnetic` (B-flux → Wb) implemented; `Power`, `SA`, `MS`, `MA`, `FarField`, `Dielectric` not yet implemented; **`postpro/port-VI.csv`** (complex port V/I/P) now emitted by Driven solver (Phase 18)
- `Boundaries.WavePort.{Offset, MaxIts, EigenTol, Verbose}` — accepted, ignored
- `Boundaries.LumpedPort.Elements` — multi-element lumped port; uses first element only
- `Solver.Driven.Save` (array) — ignored; use SaveStep integer
- `Solver.Transient.{Excitation, ExcitationFreq, ExcitationWidth}` — custom waveforms not implemented

### Remove (fields that would cause parse errors in REM)
None — REM now accepts all known Palace fields.

### Mesh files
Ensure corresponding `.msh` files exist in `examples/palace/<category>/mesh/`
or `examples/rem/<category>/mesh/`. Check:
```
find examples -name "*.msh" | sort
```
Known mesh files: `coaxial.msh`, `cpw_coax_0.msh`, `cpw_lumped_0.msh`,
`cpw_wave_0.msh`, `cylinder_hex.msh`, `cylinder_prism.msh`,
`cylinder_tet.msh`, `transmon.msh2`, `antenna.msh`, `rings.msh`,
`spheres.msh`, `adapter.msh`, `plate_2d.msh`, `sphere.msh`

## Sync Checklist

1. Pull latest Palace: `git -C /tmp/palace-src pull`
2. For each new/changed Palace JSON:
   - Copy to `examples/<category>/<name>.json`
   - Apply adaptation rules above
   - Verify mesh file exists in `examples/<category>/mesh/`
   - Test: `cargo run --example <name>` or parse with `rem-config`
3. Update this document with new files and any schema changes

## REM Schema Reference

See `crates/config/src/schema.rs` for authoritative list of supported fields.

```
Problem.Type:    Electrostatic | Magnetostatic | Eigenmode | Driven | Transient | MoM | SBR
Model.Mesh:       string (relative path to .msh file)
Model.L0:         f64 (mesh scale factor, default 1.0)
Model.Refinement: { MaxIter, Tol, Nonconformal }
Material:         { Attributes, Permittivity, Permeability, LossTan, Conductivity }  # scalar only
BoundaryTag:      PEC | PMC | Ground | Impedance | Absorbing | LumpedPort | WavePort | SurfaceCurrent
Solver.Eigenmode: { N, Tol, Target, Save }
Solver.Driven:    { MinFreq, MaxFreq, FreqStep, SaveStep, AdaptiveTol, RomOrder, Samples, CircuitSynthesis }
Solver.FarField:  { NTheta, NPhi }   # near-to-far-field transform; generates far_field.csv artifact
Solver.Transient: { Type, MaxTime, TimeStep, SaveStep }
Solver.Linear:    { Tol, MaxIter }
Solver.Mom:       { FreqMin, FreqMax, FreqStep, Equation, Basis, Alpha, SingularTol, FastSolver,
                    WallConductivity, RefImpedance, Ports, RomOrder, AmrIter, AmrTheta,
                    NearFieldSource, NearFieldProbes, TlineLength,
                    DeembedEpsEff, DeembedAlpha, Substrate }
  # Phase 19 (v0.20.0) additions:
  # RomOrder > 0 enables snapshot ROM acceleration for S-param sweeps (anchor-point Galerkin projection)
  # AmrIter > 0 enables Dörfler-marking AMR with RWG current error indicator
  # AmrTheta: Dörfler fraction (default 0.5)
  # Full Touchstone 1.0 reader (RI/MA/DB; multi-line N≥3 format) available in rem-touchstone crate
  # Phase 22 (ongoing) additions:
  # Ports[].Type: "Lumped" | "WavePort" (WavePort mode metadata accepted in MoM)
  # Ports[].Mode: WavePort mode index (mode>1 currently logs fallback warning)
  # Ports[].PairWith: differential pairing for mixed-mode S output (port-S-mixed.csv, s_params_mixed.sNp)
  # Ports[].DeembedLength: per-port reference plane shift [m]
  # DeembedEpsEff, DeembedAlpha: global de-embedding propagation model parameters
```
