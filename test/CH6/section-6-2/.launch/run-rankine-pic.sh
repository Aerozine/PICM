#!/bin/bash -l
#SBATCH --job-name=rankine-pic
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --output=rankine-pic_%j.out
#SBATCH --error=rankine-pic_%j.err

set -euo pipefail

if [ -z "${SLURM_JOB_ID:-}" ]; then
  echo "Error: submit this script with sbatch, for example:"
  echo "  sbatch PICM/test/CH6/section-6-2/.launch/run-rankine-pic.sh"
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PICM_ROOT="$(cd "${SCRIPT_DIR}/../../../../.." && pwd)"

CONFIG="${CONFIG:-test/CH6/section-6-2/rankine-pic.json}"

if type module >/dev/null 2>&1; then
  module purge
  module load releases/2021b
  module load Info0939Tools
fi

if [ -z "${PIC_BIN:-}" ]; then
  for candidate in \
    "${PICM_ROOT}/build-release/bin/PIC" \
    "${PICM_ROOT}/bin/PIC" \
    "${HOME}/project_info0939/bin/PIC"; do
    if [ -x "$candidate" ]; then
      PIC_BIN="$candidate"
      break
    fi
  done
fi

PIC_BIN="${PIC_BIN:-${PICM_ROOT}/build-release/bin/PIC}"

if [ ! -x "$PIC_BIN" ]; then
  echo "PIC executable not found or not executable: $PIC_BIN" >&2
  echo "Build first with: cmake --build build-release" >&2
  echo "Or submit with: sbatch --export=ALL,PIC_BIN=/path/to/PIC $0" >&2
  exit 1
fi

cd "$PICM_ROOT"

if [ ! -f "$CONFIG" ]; then
  echo "Config not found: ${PICM_ROOT}/$CONFIG" >&2
  echo "Override it with:" >&2
  echo "  sbatch --export=ALL,CONFIG=test/CH6/section-6-2/rankine-pic.json $0" >&2
  exit 1
fi

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
export OMP_PLACES="${OMP_PLACES:-cores}"

echo "Job info"
echo "--------"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node(s): ${SLURM_JOB_NODELIST:-unknown}"
echo "Submit dir: ${SLURM_SUBMIT_DIR:-unknown}"
echo "PICM root: ${PICM_ROOT}"
echo "Executable: ${PIC_BIN}"
echo "Config: ${CONFIG}"
echo "OMP_NUM_THREADS: ${OMP_NUM_THREADS}"
echo "Start time: $(date +"%Y-%m-%d %H:%M:%S")"
echo

srun --ntasks=1 \
  --cpus-per-task="$OMP_NUM_THREADS" \
  --cpu-bind=cores \
  "$PIC_BIN" -c "$CONFIG"

echo
echo "End time: $(date +"%Y-%m-%d %H:%M:%S")"
