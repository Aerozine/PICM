#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PICM_ROOT="${PICM_ROOT:-$(cd "${SCRIPT_DIR}/../../.." && pwd)}"

CONFIG_DIR="${CONFIG_DIR:-test/tension/section-8-8-3}"
PIC_BIN="${PIC_BIN:-${PICM_ROOT}/build-release/bin/PIC}"
LOG_DIR="${LOG_DIR:-${PICM_ROOT}/results/tension/section-8-8-3/slurm-logs}"

CPUS_PER_TASK="${CPUS_PER_TASK:-8}"
MAX_PARALLEL="${MAX_PARALLEL:-4}"
TIME_LIMIT="${TIME_LIMIT:-02:00:00}"
MEMORY="${MEMORY:-8G}"
JOB_NAME="${JOB_NAME:-picm-contact-angle}"

load_configs() {
  local -n out=$1
  mapfile -t out < <(
    find "${PICM_ROOT}/${CONFIG_DIR}" -maxdepth 1 -type f -name 'angle-*.json' \
      -printf '%P\n' | sort -V
  )
}

check_inputs() {
  if [[ ! -d "${PICM_ROOT}/${CONFIG_DIR}" ]]; then
    echo "Config directory not found: ${PICM_ROOT}/${CONFIG_DIR}" >&2
    exit 1
  fi

  if [[ ! -x "${PIC_BIN}" ]]; then
    echo "PIC executable not found or not executable: ${PIC_BIN}" >&2
    echo "Build PIC first, or run with PIC_BIN=/path/to/PIC $0" >&2
    exit 1
  fi
}

maybe_load_modules() {
  # NIC5 course environment. If the module is unavailable, keep going: the
  # executable may already have been built with all runtime libraries visible.
  if type module >/dev/null 2>&1; then
    module load Info0939Tools >/dev/null 2>&1 || true
  fi
}

run_one() {
  local -a configs=()
  load_configs configs

  local task_id="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set}"
  if (( task_id < 0 || task_id >= ${#configs[@]} )); then
    echo "Invalid SLURM_ARRAY_TASK_ID=${task_id}; ${#configs[@]} configs found." >&2
    exit 1
  fi

  local config_path="${CONFIG_DIR}/${configs[$task_id]}"

  cd "${PICM_ROOT}"
  maybe_load_modules

  export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-${CPUS_PER_TASK}}"
  export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
  export OMP_PLACES="${OMP_PLACES:-cores}"

  echo "Job ${SLURM_JOB_ID:-local}.${SLURM_ARRAY_TASK_ID:-0}"
  echo "Node(s): ${SLURM_JOB_NODELIST:-local}"
  echo "Config: ${config_path}"
  echo "PIC_BIN: ${PIC_BIN}"
  echo "OMP_NUM_THREADS: ${OMP_NUM_THREADS}"
  echo

  "${PIC_BIN}" "${config_path}"
}

submit_array() {
  check_inputs
  local -a configs=()
  load_configs configs

  if (( ${#configs[@]} == 0 )); then
    echo "No angle-*.json files found in ${PICM_ROOT}/${CONFIG_DIR}" >&2
    exit 1
  fi

  mkdir -p "${LOG_DIR}"

  if ! command -v sbatch >/dev/null 2>&1; then
    echo "sbatch not found. Run this script on the cluster login node." >&2
    echo "For a local sequential run, use: $0 --local" >&2
    exit 1
  fi

  local last_index=$(( ${#configs[@]} - 1 ))

  echo "Submitting ${#configs[@]} contact-angle simulations from ${CONFIG_DIR}"
  echo "Resources per simulation: ${CPUS_PER_TASK} CPU(s), ${MEMORY}, ${TIME_LIMIT}"
  echo "Array throttle: ${MAX_PARALLEL} concurrent simulation(s)"
  echo "Logs: ${LOG_DIR}"
  echo

  export PICM_ROOT CONFIG_DIR PIC_BIN LOG_DIR
  export CPUS_PER_TASK MAX_PARALLEL TIME_LIMIT MEMORY JOB_NAME

  sbatch \
    --job-name="${JOB_NAME}" \
    --ntasks=1 \
    --cpus-per-task="${CPUS_PER_TASK}" \
    --mem="${MEMORY}" \
    --time="${TIME_LIMIT}" \
    --export=ALL \
    --array="0-${last_index}%${MAX_PARALLEL}" \
    --output="${LOG_DIR}/%x_%A_%a.out" \
    --error="${LOG_DIR}/%x_%A_%a.err" \
    "$0" --run-array
}

run_local() {
  check_inputs
  local -a configs=()
  load_configs configs

  if (( ${#configs[@]} == 0 )); then
    echo "No angle-*.json files found in ${PICM_ROOT}/${CONFIG_DIR}" >&2
    exit 1
  fi

  cd "${PICM_ROOT}"
  export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${CPUS_PER_TASK}}"

  for config in "${configs[@]}"; do
    echo
    echo "==> Running ${CONFIG_DIR}/${config}"
    "${PIC_BIN}" "${CONFIG_DIR}/${config}"
  done

  echo
  echo "All contact-angle simulations completed."
}

case "${1:-}" in
  --run-array)
    check_inputs
    run_one
    ;;
  --local)
    run_local
    ;;
  ""|--submit)
    submit_array
    ;;
  *)
    echo "Usage: $0 [--submit|--local|--run-array]" >&2
    echo "Override resources with CPUS_PER_TASK=8 MAX_PARALLEL=4 TIME_LIMIT=02:00:00 MEMORY=8G." >&2
    exit 2
    ;;
esac
