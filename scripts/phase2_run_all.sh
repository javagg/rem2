#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$ROOT_DIR/output/phase2"
LOG_DIR="$OUT_DIR/logs"
REPORT="$OUT_DIR/report.md"

RUN_NATIVE_EXAMPLES="${RUN_NATIVE_EXAMPLES:-1}"
RUN_NUMERIC_TESTS="${RUN_NUMERIC_TESTS:-1}"
RUN_FEMRS_MPI="${RUN_FEMRS_MPI:-1}"
RUN_WASM_BUILD="${RUN_WASM_BUILD:-1}"

mkdir -p "$LOG_DIR"

run_step() {
  local name="$1"
  local cmd="$2"
  local log_file="$LOG_DIR/${name}.log"

  echo "[phase2] >>> $name"
  if bash -lc "$cmd" >"$log_file" 2>&1; then
    echo "[phase2] PASS $name"
    echo "- [PASS] $name" >>"$REPORT"
    return 0
  fi

  echo "[phase2] FAIL $name (see $log_file)"
  echo "- [FAIL] $name" >>"$REPORT"
  return 1
}

run_native_examples() {
  local examples=(
    "examples/adapter/adapter.json"
    "examples/antenna/antenna.json"
    "examples/coaxial/coaxial.json"
    "examples/cpw/cpw.json"
    "examples/cylinder/cylinder.json"
    "examples/parallel_plate/parallel_plate.json"
    "examples/rings/rings.json"
    "examples/sbr_sphere/sbr_sphere.json"
    "examples/spheres/spheres.json"
    "examples/transmon/transmon.json"
  )

  for cfg in "${examples[@]}"; do
    local name
    name="native_$(basename "$cfg" .json)"
    run_step "$name" "cd '$ROOT_DIR' && cargo run --release -p rem-cli -- '$cfg' -v"
  done
}

run_numeric_tests() {
  run_step "test_electrostatic_palace_spheres" \
    "cd '$ROOT_DIR' && cargo test --release -p rem-electrostatic --test palace_spheres -- --nocapture"
  run_step "test_electrostatic_parallel_plate" \
    "cd '$ROOT_DIR' && cargo test --release -p rem-electrostatic palace_parallel_plate -- --nocapture"
  run_step "test_magnetostatic_palace_rings" \
    "cd '$ROOT_DIR' && cargo test --release -p rem-magnetostatic --test palace_rings -- --nocapture"
  run_step "test_mom_mie_validation" \
    "cd '$ROOT_DIR' && cargo test --release -p rem-mom --test mie_validation -- --nocapture"
  run_step "test_sbr_mie_validation" \
    "cd '$ROOT_DIR' && cargo test --release -p rem-sbr --test mie_validation -- --nocapture"
}

run_femrs_mpi_examples() {
  local femrs_dir="$ROOT_DIR/vendor/fem-rs"
  local examples=("pex1_poisson" "pex2_mixed_darcy" "pex3_maxwell" "pex5_darcy")

  run_step "mpi_env_check" "command -v mpirun >/dev/null && command -v mpicc >/dev/null && mpirun --version | head -n 1"

  for ex in "${examples[@]}"; do
    run_step "mpi_${ex}" \
      "cd '$femrs_dir' && mpirun -n 4 cargo run --release -p fem-examples --example '$ex' --features fem-parallel/mpi"
  done
}

run_wasm_builds() {
  run_step "wasm_target" "rustup target add wasm32-unknown-unknown"
  run_step "build_rem_wasm" "cd '$ROOT_DIR' && cargo build --release --target wasm32-unknown-unknown -p rem-wasm"
  run_step "build_yew_app" "cd '$ROOT_DIR/crates/yew-app' && trunk build --release"
}

{
  echo "# Phase 2 Verification Report"
  echo ""
  echo "- Date: $(date '+%Y-%m-%d %H:%M:%S %z')"
  echo "- Root: $ROOT_DIR"
  echo ""
  echo "## Summary"
} >"$REPORT"

status=0

if [[ "$RUN_NATIVE_EXAMPLES" == "1" ]]; then
  run_native_examples || status=1
fi

if [[ "$RUN_NUMERIC_TESTS" == "1" ]]; then
  run_numeric_tests || status=1
fi

if [[ "$RUN_FEMRS_MPI" == "1" ]]; then
  run_femrs_mpi_examples || status=1
fi

if [[ "$RUN_WASM_BUILD" == "1" ]]; then
  run_wasm_builds || status=1
fi

echo "" >>"$REPORT"
if [[ $status -eq 0 ]]; then
  echo "Overall: PASS" >>"$REPORT"
else
  echo "Overall: FAIL" >>"$REPORT"
fi

echo "[phase2] report: $REPORT"
exit $status
