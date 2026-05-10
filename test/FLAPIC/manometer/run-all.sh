#!/bin/bash -l
#SBATCH --job-name=FLAPIC-manometer
#SBATCH --partition=hmem
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --time=02:00:00

set -euo pipefail

scriptDir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
picmRoot="${PICM_ROOT:-$(cd "${scriptDir}/../../.." && pwd)}"
startDir="$(pwd)"

numCpuCores="${SLURM_CPUS_PER_TASK:-64}"
executable="${PIC_BIN:-${picmRoot}/build-release/bin/PIC}"
configDir="test/FLAPIC/manometer"
logDir="${picmRoot}/results/FLAPIC/manometer/slurm-logs/${SLURM_JOB_ID:-local}"

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
  echo "Or run with: PIC_BIN=/path/to/PIC $0" >&2
  exit 1
fi

mkdir -p "${logDir}"

export OMP_NUM_THREADS="${numCpuCores}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
export OMP_PLACES="${OMP_PLACES:-cores}"

echo "Job info"
echo "--------"
echo
echo "    Job ID: ${SLURM_JOB_ID:-local}"
echo " Node list: ${SLURM_JOB_NODELIST:-local}"
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

  if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v srun >/dev/null 2>&1; then
    srun --ntasks=1 \
      --cpus-per-task="${numCpuCores}" \
      --cpu-bind=cores \
      "${executable}" -c "${configPath}" >"${output}" 2>"${error}"
  else
    "${executable}" -c "${configPath}" >"${output}" 2>"${error}"
  fi

  endTime="$(date +%s.%N)"
  elapsed="$(awk -v start="${startTime}" -v end="${endTime}" 'BEGIN { printf "%.6f", end - start }')"

  echo "[$(date +"%d-%m-%Y %H:%M:%S")] Done ${config} in ${elapsed} seconds"
  echo
done

cd "${startDir}"

echo "End time: $(date +"%d-%m-%Y %H:%M:%S")"
echo "All FLAPIC manometer simulations completed."
