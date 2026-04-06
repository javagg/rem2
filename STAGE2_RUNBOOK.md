# Phase 2 Runbook

This stage turns the first-round fixes into a repeatable regression pipeline.

## What it covers

- Native run of all example configs in `examples/*/*.json`
- Numerical baseline tests for electrostatic, magnetostatic, MoM and SBR
- MPI validation in `vendor/fem-rs` with 4 ranks
- WASM build validation for `rem-wasm` and `crates/yew-app`

## One-command execution

```bash
bash scripts/phase2_run_all.sh
```

The script writes logs and report files into:

- `output/phase2/logs/*.log`
- `output/phase2/report.md`

## Optional toggles

Disable a section by environment variable:

```bash
RUN_NATIVE_EXAMPLES=0 bash scripts/phase2_run_all.sh
RUN_NUMERIC_TESTS=0 bash scripts/phase2_run_all.sh
RUN_FEMRS_MPI=0 bash scripts/phase2_run_all.sh
RUN_WASM_BUILD=0 bash scripts/phase2_run_all.sh
```

## Expected prerequisites

- Rust toolchain and cargo available
- `trunk` installed for yew build
- OpenMPI/MPICH available for MPI section (`mpirun`, `mpicc`)
