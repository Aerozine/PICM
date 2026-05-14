#!/usr/bin/env bash
set -uo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIC_BIN="${ROOT_DIR}/build-release/bin/PIC"

cd "${ROOT_DIR}" || exit 1

echo "== Build release =="
cmake -B build-release -G Ninja -DCMAKE_BUILD_TYPE=Release || exit 1
cmake --build build-release || exit 1

echo
echo "== Run test/CH* cases =="

shopt -s nullglob
CH_DIRS=(test/CH*)
shopt -u nullglob

if [ "${#CH_DIRS[@]}" -eq 0 ]; then
  echo "No test/CH* directory found."
  exit 0
fi

mapfile -t CONFIGS < <(find "${CH_DIRS[@]}" -type f -name '*.json' | sort)

if [ "${#CONFIGS[@]}" -eq 0 ]; then
  echo "No JSON config found under test/CH*."
  exit 0
fi

failed=0
failed_configs=()

for config in "${CONFIGS[@]}"; do
  echo
  echo "Running: ${config}"
  if ! "${PIC_BIN}" "${config}"; then
    failed=$((failed + 1))
    failed_configs+=("${config}")
    echo "FAILED: ${config}"
  fi
done

echo
if [ "${failed}" -eq 0 ]; then
  echo "All ${#CONFIGS[@]} CH cases completed."
  exit 0
fi

echo "${failed} / ${#CONFIGS[@]} CH cases failed:"
for config in "${failed_configs[@]}"; do
  echo "  - ${config}"
done
exit 1
