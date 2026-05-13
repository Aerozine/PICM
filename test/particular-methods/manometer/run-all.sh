#!/bin/bash -l
#SBATCH --job-name=FLAPIC-manometer
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=02:00:00
#SBATCH --output=FLAPIC-manometer-%j.out
#SBATCH --error=FLAPIC-manometer-%j.err

set -euo pipefail

# -----------------------------------------------------------------------------
# Safety guard: if this script is executed directly on the login node, it only
# submits itself to Slurm and exits. The simulations never run on the login node.
# -----------------------------------------------------------------------------
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  if ! command -v sbatch >/dev/null 2>&1; then
    echo "Error: this script must be submitted with Slurm, but sbatch was not found." >&2
    exit 1
  fi

  echo "Not inside a Slurm allocation. Submitting this script with sbatch..."
  sbatch "$0" "$@"
  exit 0
fi

scriptDir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
picmRoot="${PICM_ROOT:-$(cd "${scriptDir}/../../.." && pwd)}"
startDir="$(pwd)"

numCpuCores="${SLURM_CPUS_PER_TASK:-64}"
executable="${PIC_BIN:-${picmRoot}/build-release/bin/PIC}"
configDir="test/FLAPIC/manometer"
logDir="${picmRoot}/results/FLAPIC/manometer/slurm-logs/${SLURM_JOB_ID}"

configs=(
  "manometer-pic.json"
  "manometer-flip.json"
  "manometer-apic.json"
)

if type module >/dev/null 2>&1; then
  module purge
  module load releases/2021b
  module load Info0939Tools
fi

if [[ ! -x "${executable}" ]]; then
  echo "PIC executable not found or not executable: ${executable}" >&2
  echo "Build it first with: cmake --build build-release" >&2
  echo "Or submit with: PIC_BIN=/path/to/PIC sbatch $0" >&2
  exit 1
fi

mkdir -p "${logDir}"

export OMP_NUM_THREADS="${numCpuCores}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
export OMP_PLACES="${OMP_PLACES:-cores}"

echo "Job info"
echo "--------"
echo "    Job ID: ${SLURM_JOB_ID}"
echo " Node list: ${SLURM_JOB_NODELIST:-unknown}"
echo " cpus-per-task: ${numCpuCores}"
echo " executable: ${executable}"
echo " log dir: ${logDir}"
echo "Start time: $(date +"%d-%m-%Y %H:%M:%S")"
echo

cd "${picmRoot}"

for config in "${configs[@]}"; do
  configPath="${configDir}/${config}"
  name="${config%.json}"
  output="${logDir}/${name}.out"
  error="${logDir}/${name}.err"

  if [[ ! -f "${configPath}" ]]; then
    echo "Config not found: ${configPath}" >&2
    exit 1
  fi

  echo "[$(date +"%d-%m-%Y %H:%M:%S")] Running ${configPath}"
  echo "  stdout: ${output}"
  echo "  stderr: ${error}"

  startTime="$(date +%s.%N)"

  srun --ntasks=1 \
       --cpus-per-task="${numCpuCores}" \
       --cpu-bind=cores \
       "${executable}" -c "${configPath}" >"${output}" 2>"${error}"

  endTime="$(date +%s.%N)"
  elapsed="$(awk -v start="${startTime}" -v end="${endTime}" 'BEGIN { printf "%.6f", end - start }')"

  echo "[$(date +"%d-%m-%Y %H:%M:%S")] Done ${config} in ${elapsed} seconds"
  echo
done

cd "${startDir}"

echo "End time: $(date +"%d-%m-%Y %H:%M:%S")"
echo "All FLAPIC manometer simulations completed."
