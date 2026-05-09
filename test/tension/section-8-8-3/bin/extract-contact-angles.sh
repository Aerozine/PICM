#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PICM_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
WORKSPACE_ROOT="$(cd "${PICM_ROOT}/.." && pwd)"

if [[ -z "${PYTHON:-}" && -x "${PICM_ROOT}/.venv/bin/python" ]]; then
  PYTHON="${PICM_ROOT}/.venv/bin/python"
else
  PYTHON="${PYTHON:-python3}"
fi

VALUE_SCRIPT="${VALUE_SCRIPT:-${WORKSPACE_ROOT}/postProd/scripts/value.py}"
CONTACT_SCRIPT="${CONTACT_SCRIPT:-${WORKSPACE_ROOT}/postProd/scripts/contact-angles.py}"
CSV_DIR="${CSV_DIR:-${WORKSPACE_ROOT}/postProd/results/outputs/tensions/contact-angles/personal}"
TABLE_OUT="${TABLE_OUT:-${CSV_DIR}/contact-angles-table.tex}"
SUMMARY_OUT="${SUMMARY_OUT:-${CSV_DIR}/contact-angles-summary.txt}"
FIELD_NAME="${FIELD_NAME:-label}"
CONTACT_BAND="${CONTACT_BAND:-1}"

if [[ -n "${ANGLES_OVERRIDE:-}" ]]; then
  read -r -a ANGLES <<< "${ANGLES_OVERRIDE}"
else
  ANGLES=(30 45 60 90 120 135 150)
fi

if [[ ! -f "${VALUE_SCRIPT}" ]]; then
  echo "value.py not found: ${VALUE_SCRIPT}" >&2
  exit 1
fi

if [[ ! -f "${CONTACT_SCRIPT}" ]]; then
  echo "contact-angles.py not found: ${CONTACT_SCRIPT}" >&2
  exit 1
fi

mkdir -p "${CSV_DIR}"

common_dx=""
common_dy=""
extracted=0
csv_files=()

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
    float(cfg["dx"]),
    float(cfg["dy"]),
))' "${config}")"
  IFS=$'\t' read -r folder dx dy <<< "${config_info}"

  if [[ -z "${common_dx}" ]]; then
    common_dx="${dx}"
    common_dy="${dy}"
  elif [[ "${dx}" != "${common_dx}" || "${dy}" != "${common_dy}" ]]; then
    echo "Warning: angle ${angle} has dx=${dx}, dy=${dy}; using dx=${common_dx}, dy=${common_dy} for the final table." >&2
  fi

  if [[ "${folder}" = /* ]]; then
    result_dir="${folder}"
  else
    result_dir="${PICM_ROOT}/${folder}"
  fi

  pvd="${result_dir}/${FIELD_NAME}.pvd"
  csv="${CSV_DIR}/contact-angle-${angle}.csv"

  if [[ ! -f "${pvd}" ]]; then
    echo "Skipping angle ${angle}: missing ${pvd}" >&2
    continue
  fi

  echo
  echo "==> Extracting final ${FIELD_NAME} field for ${angle} deg"
  "${PYTHON}" "${VALUE_SCRIPT}" "${pvd}" \
    --field "${FIELD_NAME}" \
    --mode field \
    --last \
    --out "${csv}"

  extracted=$((extracted + 1))
  csv_files+=("${csv}")
done

if [[ "${extracted}" -eq 0 ]]; then
  echo "No contact-angle CSV was extracted." >&2
  exit 1
fi

echo
echo "==> Measuring contact angles and writing LaTeX table"
"${PYTHON}" "${CONTACT_SCRIPT}" \
  --field "${FIELD_NAME}" \
  --dx "${common_dx}" \
  --dy "${common_dy}" \
  --contact-band "${CONTACT_BAND}" \
  --output "${TABLE_OUT}" \
  "${csv_files[@]}" \
  > "${SUMMARY_OUT}"

cat "${SUMMARY_OUT}"

echo
echo "CSV files: ${CSV_DIR}"
echo "LaTeX table: ${TABLE_OUT}"
echo "Summary: ${SUMMARY_OUT}"
