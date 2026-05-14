#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PICM_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
WORKSPACE_ROOT="$(cd "${PICM_ROOT}/.." && pwd)"

PIC="${PICM_ROOT}/build-release/bin/PIC"
PYTHON="${PYTHON:-python3}"
VALUE_PY="${WORKSPACE_ROOT}/postProd/scripts/value.py"
CONTACT_PY="${WORKSPACE_ROOT}/postProd/scripts/contact-angles.py"

OUT_DIR="${WORKSPACE_ROOT}/postProd/results/outputs/tensions/contact-angles/section-8-8-3"
TABLE_OUT="${OUT_DIR}/contact-angles-table.tex"
SUMMARY_OUT="${OUT_DIR}/contact-angles-summary.txt"
PLOT_OUT="${WORKSPACE_ROOT}/postProd/results/images/Tensions/contact-angles/section-8-8-3/contact-angles-validation"

ANGLES="30 45 60 90 120 135 150"

cd "${PICM_ROOT}"
mkdir -p "${OUT_DIR}"

for angle in ${ANGLES}; do
  echo
  echo "==> Simulation angle ${angle}"
  "${PIC}" "test/tension/section-8-8-3/angle-${angle}.json"
done

CSV_FILES=""

for angle in ${ANGLES}; do
  echo
  echo "==> Extraction CSV angle ${angle}"
  pvd="${PICM_ROOT}/results/tension/angle-${angle}/label.pvd"
  csv="${OUT_DIR}/contact-angle-${angle}.csv"

  "${PYTHON}" "${VALUE_PY}" "${pvd}" \
    --field label \
    --mode field \
    --last \
    --out "${csv}"

  CSV_FILES="${CSV_FILES} ${csv}"
done

echo
echo "==> Calcul des angles de contact"
"${PYTHON}" "${CONTACT_PY}" \
  --field label \
  --dx 1e-6 \
  --dy 1e-6 \
  --contact-band 1 \
  --output "${TABLE_OUT}" \
  --plot-output "${PLOT_OUT}" \
  ${CSV_FILES} \
  > "${SUMMARY_OUT}"

cat "${SUMMARY_OUT}"

echo
echo "Resultats des simulations :"
for angle in ${ANGLES}; do
  echo "  ${PICM_ROOT}/results/tension/angle-${angle}"
done

echo
echo "Resultats du post-traitement :"
echo "  CSV     : ${OUT_DIR}"
echo "  Table   : ${TABLE_OUT}"
echo "  Resume  : ${SUMMARY_OUT}"
echo "  Plot 1  : ${PLOT_OUT}-measured-vs-ideal.pdf / ${PLOT_OUT}-measured-vs-ideal.png"
echo "  Plot 2  : ${PLOT_OUT}-relative-error.pdf / ${PLOT_OUT}-relative-error.png"
