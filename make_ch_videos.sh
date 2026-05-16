#!/usr/bin/env bash

set -u
set -o pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
OUT_ROOT="${OUT_ROOT:-${ROOT_DIR}/videos_out}"
LOG_DIR="${LOG_DIR:-${OUT_ROOT}/logs}"

mkdir -p "${OUT_ROOT}" "${LOG_DIR}"

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

print_items() {
  if (($# == 0)); then
    echo "  aucun"
    return
  fi
  printf '  %s\n' "$@"
}

mapfile -t CONFIGS < <(find "${ROOT_DIR}/test"/CH* -type f -name '*.json' | sort)

echo "== Génération des vidéos =="
for config in "${CONFIGS[@]}"; do
  rel_config="${config#${ROOT_DIR}/}"
  rel_case="${rel_config#test/}"
  rel_case="${rel_case%.json}"
  result_dir="$(read_json_field "${config}" "folder")"
  title="$(basename "${config}" .json)"
  key="$(case_key "${config}")"

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
      renderer=(
        "${PYTHON_BIN}" "${ROOT_DIR}/videos/videos.py" "${input_pvd}"
        --type field
        --out "${output_mp4}"
        --title "${title}"
      )
      ;;
    test/CH6/*|test/CH7/*)
      input_pvd="${abs_result_dir}/particles.pvd"
      renderer=(
        "${PYTHON_BIN}" "${ROOT_DIR}/videos/particler.py" "${input_pvd}"
        --out "${output_mp4}"
        --title "${title}"
      )
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
  if ! "${renderer[@]}" >"${log_path}" 2>&1; then
    echo "  [FAIL] voir ${log_path#${ROOT_DIR}/}"
    VIDEO_FAILURES+=("${rel_config}")
  fi
done

echo "== Résumé =="
echo "Vidéos échouées: ${#VIDEO_FAILURES[@]}"
print_items "${VIDEO_FAILURES[@]}"
echo "Vidéos ignorées: ${#VIDEO_SKIPS[@]}"
print_items "${VIDEO_SKIPS[@]}"
echo "Sortie: ${OUT_ROOT#${ROOT_DIR}/}"

exit 0
