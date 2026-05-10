#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PICM_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
PIC_BIN="${PIC_BIN:-${PICM_ROOT}/build-release/bin/PIC}"

if [[ ! -x "${PIC_BIN}" ]]; then
  echo "PIC executable not found or not executable: ${PIC_BIN}" >&2
  echo "Build PIC first, or run with PIC_BIN=/path/to/PIC $0" >&2
  exit 1
fi

CONFIG_DIR="test/tension/section-8-8-3/personal"
CONFIGS=(
  "angle-30.json"
  "angle-45.json"
  "angle-60.json"
  "angle-90.json"
  "angle-120.json"
  "angle-135.json"
  "angle-150.json"
)

cd "${PICM_ROOT}"

for config in "${CONFIGS[@]}"; do
  config_path="${CONFIG_DIR}/${config}"
  echo
  echo "==> Running ${config_path}"
  "${PIC_BIN}" "${config_path}"
done

echo
echo "All personal contact-angle simulations completed."
