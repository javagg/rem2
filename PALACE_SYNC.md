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

All Palace JSON files are present in `examples/` with identical directory
structure and filenames, adapted for REM's config schema.

| Subdirectory | Files | Solver Type | REM Status |
|-------------|-------|-------------|-----------|
| `adapter/` | `hybrid.json` | Driven | ✅ |
| `antenna/` | `antenna_halfwave_dipole.json`, `antenna_short_dipole.json` | Driven | ✅ |
| `coaxial/` | `coaxial_matched.json`, `coaxial_open.json`, `coaxial_short.json` | Transient | ✅ |
| `cpw/` | 8 files | Driven/Eigenmode | ✅ |
| `cylinder/` | `cavity_impedance.json`, `cavity_pec.json`, `driven_wave.json`, `floquet.json`, `waveguide.json` | Eigenmode/Driven | ✅ |
| `rings/` | `rings.json` | Magnetostatic | ✅ |
| `spheres/` | `spheres.json` | MoM | ✅ |
| `transmon/` | `transmon_amr.json`, `transmon_coarse.json` | Eigenmode | ✅ |

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
- `Problem.OutputFormats.GridFunction` — accepted, not implemented
- `Solver.Device: "CPU"` — REM is CPU-only; ignored with warning
- `Solver.Linear.KSPType` — only GMRES is supported; ignored with warning
- `Solver.Linear.MGLevels` — algebraic multigrid not implemented; ignored
- `Solver.Linear.ComplexCoarseSolve` — complex coarse-grid solve not implemented; ignored
- `Domains.Materials[*].MaterialAxes` — anisotropic εᵣ/μᵣ not implemented; uses first component
- `Domains.Materials[*].LossTan` (array form) — anisotropic loss not implemented; uses first element
- `Domains.Postprocessing.Energy`, `Domains.Postprocessing.Probe` — not implemented
- `Domains.CurrentDipole` — Hertzian dipole source not implemented
- `Boundaries.Periodic` — periodic/Floquet BCs not implemented
- `Boundaries.Postprocessing.{SurfaceFlux, FarField, Dielectric}` — not implemented
- `Boundaries.WavePort.{Offset, MaxIts, EigenTol, Verbose}` — accepted, ignored
- `Boundaries.LumpedPort.Elements` — multi-element lumped port; uses first element only
- `Solver.Driven.Samples` — complex sampling schedules; ignored; use MinFreq/MaxFreq/FreqStep
- `Solver.Driven.Save` (array) — ignored; use SaveStep integer
- `Solver.Transient.{Excitation, ExcitationFreq, ExcitationWidth}` — custom waveforms not implemented

### Remove (fields that would cause parse errors in REM)
None — REM now accepts all known Palace fields.

### Mesh files
Ensure corresponding `.msh` files exist in `examples/<category>/mesh/`. Check:
```
find examples -name "*.msh" | sort
```
Known mesh files: `coaxial.msh`, `cpw_coax_0.msh`, `cpw_lumped_0.msh`,
`cpw_wave_0.msh`, `cylinder_hex.msh`, `cylinder_prism.msh`,
`cylinder_tet.msh`, `transmon.msh2`, `antenna.msh`, `rings.msh`,
`spheres.msh`, `adapter.msh`, `plate_2d.msh`, `slab_2d.msh`,
`coaxial_2d.msh`, `sphere.msh`

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
Solver.Driven:    { MinFreq, MaxFreq, FreqStep, SaveStep, AdaptiveTol }
Solver.Transient: { Type, MaxTime, TimeStep, SaveStep }
Solver.Linear:    { Tol, MaxIter }
```
