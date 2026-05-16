#!/usr/bin/env bash

set -u
set -o pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIC_BIN="${PIC_BIN:-${ROOT_DIR}/build-release/bin/PIC}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
LOG_DIR="${LOG_DIR:-${ROOT_DIR}/results/logs/ch-runs}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/videos_out}"

mkdir -p "${LOG_DIR}" "${OUT_ROOT}"

declare -a SIM_FAILURES=()
declare -a VIDEO_FAILURES=()
declare -a VIDEO_SKIPS=()

read_json_field() {
  local json_path="$1"
  local field="$2"
  "${PYTHON_BIN}" -c '
import json
import sys

with open(sys.argv[1], encoding="utf-8") as fh:
    data = json.load(fh)
value = data.get(sys.argv[2], "")
print(value if isinstance(value, str) else "")
' "${json_path}" "${field}"
}

case_key() {
  local path="$1"
  path="${path#${ROOT_DIR}/}"
  printf '%s\n' "${path%.json}" | tr '/ ' '__'
}

run_logged() {
  local log_path="$1"
  shift
  "$@" >"${log_path}" 2>&1
}

print_items() {
  if (($# == 0)); then
    echo "  aucun"
    return
  fi
  printf '  %s\n' "$@"
}

echo "== Build =="
if ! cmake -B "${ROOT_DIR}/build-release" -G Ninja -DCMAKE_BUILD_TYPE=Release \
    >"${LOG_DIR}/build-configure.log" 2>&1; then
  echo "[WARN] configuration CMake échouée; je continue avec le binaire existant s'il y en a un."
fi
if ! cmake --build "${ROOT_DIR}/build-release" >"${LOG_DIR}/build.log" 2>&1; then
  echo "[WARN] build échoué; je continue avec le binaire existant s'il y en a un."
fi

mapfile -t CONFIGS < <(find "${ROOT_DIR}/test"/CH* -type f -name '*.json' | sort)

echo "== Simulations =="
for config in "${CONFIGS[@]}"; do
  rel_config="${config#${ROOT_DIR}/}"
  key="$(case_key "${config}")"
  log_path="${LOG_DIR}/${key}.simulation.log"
  echo "[SIM] ${rel_config}"

  if [[ ! -x "${PIC_BIN}" ]]; then
    echo "  [FAIL] binaire introuvable: ${PIC_BIN}"
    SIM_FAILURES+=("${rel_config} (binaire introuvable)")
    continue
  fi

  if ! run_logged "${log_path}" "${PIC_BIN}" "${config}"; then
    echo "  [FAIL] voir ${log_path#${ROOT_DIR}/}"
    SIM_FAILURES+=("${rel_config}")
  fi
done

echo "== Vidéos =="
for config in "${CONFIGS[@]}"; do
  rel_config="${config#${ROOT_DIR}/}"
  rel_case="${rel_config#test/}"
  rel_case="${rel_case%.json}"
  key="$(case_key "${config}")"
  result_dir="$(read_json_field "${config}" "folder")"
  title="$(basename "${config}" .json)"

  if [[ -z "${result_dir}" ]]; then
    echo "[SKIP] ${rel_config}: champ folder absent"
    VIDEO_SKIPS+=("${rel_config} (folder absent)")
    continue
  fi

  abs_result_dir="${ROOT_DIR}/${result_dir}"
  output_mp4="${OUT_ROOT}/${rel_case}.mp4"
  mkdir -p "$(dirname "${output_mp4}")"
  case "${rel_config}" in
    test/CH4/*)
      input_pvd="${abs_result_dir}/normVelocity.pvd"
      renderer=("${PYTHON_BIN}" "${ROOT_DIR}/videos/videos.py" "${input_pvd}" --type field --out "${output_mp4}" --title "${title}")
      ;;
    test/CH6/*|test/CH7/*)
      input_pvd="${abs_result_dir}/particles.pvd"
      renderer=("${PYTHON_BIN}" "${ROOT_DIR}/videos/particler.py" "${input_pvd}" --out "${output_mp4}" --title "${title}")
      ;;
    *)
      echo "[SKIP] ${rel_config}: chapitre sans règle vidéo"
      VIDEO_SKIPS+=("${rel_config} (chapitre sans règle vidéo)")
      continue
      ;;
  esac

  if [[ ! -f "${input_pvd}" ]]; then
    echo "[SKIP] ${rel_config}: ${input_pvd#${ROOT_DIR}/} introuvable"
    VIDEO_SKIPS+=("${rel_config} (${input_pvd#${ROOT_DIR}/} introuvable)")
    continue
  fi

  log_path="${LOG_DIR}/${key}.video.log"
  echo "[VID] ${rel_config} -> ${output_mp4#${ROOT_DIR}/}"
  if ! run_logged "${log_path}" "${renderer[@]}"; then
    echo "  [FAIL] voir ${log_path#${ROOT_DIR}/}"
    VIDEO_FAILURES+=("${rel_config}")
  fi
done

echo "== Résumé =="
echo "Simulations échouées: ${#SIM_FAILURES[@]}"
print_items "${SIM_FAILURES[@]}"
echo "Vidéos échouées: ${#VIDEO_FAILURES[@]}"
print_items "${VIDEO_FAILURES[@]}"
echo "Vidéos ignorées: ${#VIDEO_SKIPS[@]}"
print_items "${VIDEO_SKIPS[@]}"
echo "Sortie vidéos: ${OUT_ROOT#${ROOT_DIR}/}"

exit 0
