#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PICM_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
WORKSPACE_ROOT="$(cd "${PICM_ROOT}/.." && pwd)"

PARTICLER="${PARTICLER:-${WORKSPACE_ROOT}/postProd/scripts/particler.py}"
OUT_DIR="${OUT_DIR:-${WORKSPACE_ROOT}/postProd/videos/tensions/section-8-8-3/personal}"
if [[ -z "${PYTHON:-}" && -x "${PICM_ROOT}/.venv/bin/python" ]]; then
  PYTHON="${PICM_ROOT}/.venv/bin/python"
else
  PYTHON="${PYTHON:-python3}"
fi
export MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp/matplotlib-picm}"

FPS="${FPS:-20}"
CMAP="${CMAP:-viridis}"
MODE="${MODE:-speed}"
WIDTH="${WIDTH:-1280}"
HEIGHT="${HEIGHT:-720}"
WORKERS="${WORKERS:-}"
SAMPLE="${SAMPLE:-1}"

ANGLES=(30 45 60 90 120 135 150)

if [[ ! -f "${PARTICLER}" ]]; then
  echo "particler.py not found: ${PARTICLER}" >&2
  exit 1
fi

if ! command -v ffmpeg >/dev/null 2>&1; then
  echo "ffmpeg not found in PATH; particler.py needs it to encode MP4 videos." >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"

for angle in "${ANGLES[@]}"; do
  config="${SCRIPT_DIR}/angle-${angle}.json"
  if [[ ! -f "${config}" ]]; then
    echo "Skipping angle ${angle}: missing ${config}" >&2
    continue
  fi

  config_info="$("${PYTHON}" -c 'import json, sys
cfg = json.load(open(sys.argv[1]))
print("{}\t{:.17g}\t{:.17g}".format(
    cfg["folder"],
    float(cfg["nx"]) * float(cfg["dx"]),
    float(cfg["ny"]) * float(cfg["dy"]),
))' "${config}")"
  IFS=$'\t' read -r folder domain_x domain_y <<< "${config_info}"

  if [[ "${folder}" = /* ]]; then
    result_dir="${folder}"
  else
    result_dir="${PICM_ROOT}/${folder}"
  fi
  pvd="${result_dir}/particles.pvd"
  out="${OUT_DIR}/contact-angle-$(printf "%03d" "${angle}")deg.mp4"
  xlim_min="${XLIM_MIN:-0}"
  xlim_max="${XLIM_MAX:-${domain_x}}"
  ylim_min="${YLIM_MIN:-0}"
  ylim_max="${YLIM_MAX:-${domain_y}}"

  if [[ ! -f "${pvd}" ]]; then
    echo "Skipping angle ${angle}: missing ${pvd}" >&2
    continue
  fi

  title="Contact angle ${angle} deg"
  args=(
    "${PARTICLER}"
    "${pvd}"
    "--out" "${out}"
    "--fps" "${FPS}"
    "--cmap" "${CMAP}"
    "--mode" "${MODE}"
    "--width" "${WIDTH}"
    "--height" "${HEIGHT}"
    "--sample" "${SAMPLE}"
    "--xlim" "${xlim_min}" "${xlim_max}"
    "--ylim" "${ylim_min}" "${ylim_max}"
    "--title" "${title}"
  )

  if [[ -n "${WORKERS}" ]]; then
    args+=("--workers" "${WORKERS}")
  fi

  echo
  echo "==> Rendering ${out}"
  "${PYTHON}" "${args[@]}"
done

echo
echo "Videos written to: ${OUT_DIR}"
