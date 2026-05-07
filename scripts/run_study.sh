#!/usr/bin/env bash
# run_study.sh — build (CPU), run all studies, generate plots
#
# Usage (from project root):
#   bash scripts/run_study.sh [config]

set -euo pipefail

CONFIG="${1:-test/PIC/free.json}"
BINARY="./build/bin/PIC"
BINARY_DBG="./build-debug/bin/PIC"
BUILD_DIR="build"
BUILD_DBG="build-debug"
SCRIPTS="scripts"
PYTHON="${PYTHON:-python3}"

PPC_MIN=1
PPC_MAX=3 #10
PPC_REPEATS=2 #25
METHODS="vanilla_pic" #,flip,apic,semilagrangian"
SOLVERS="red_black_gauss_seidel" #,cg,miccg0"

NCPU=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4)

# ── 1a. release build ─────────────────────────────────────────────────────────
echo "[build] release build in $BUILD_DIR/ …"
mkdir -p "$BUILD_DIR"
cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release -DUSE_GPU=OFF \
    > "$BUILD_DIR/cmake.log" 2>&1 \
    || { echo "[error] CMake failed:"; cat "$BUILD_DIR/cmake.log"; exit 1; }
cmake --build "$BUILD_DIR" -j"$NCPU" \
    || { echo "[error] release build failed"; exit 1; }
[[ -x "$BINARY" ]] || { echo "[error] binary not found: $BINARY"; exit 1; }
echo "[build] release done → $BINARY"

# ── 1b. debug build (for solver iteration counts) ─────────────────────────────
echo "[build] debug build in $BUILD_DBG/ …"
mkdir -p "$BUILD_DBG"
cmake -S . -B "$BUILD_DBG" -DCMAKE_BUILD_TYPE=Debug -DUSE_GPU=OFF \
    > "$BUILD_DBG/cmake.log" 2>&1 \
    || { echo "[error] CMake failed:"; cat "$BUILD_DBG/cmake.log"; exit 1; }
cmake --build "$BUILD_DBG" -j"$NCPU" \
    || { echo "[error] debug build failed"; exit 1; }
[[ -x "$BINARY_DBG" ]] || { echo "[error] debug binary not found: $BINARY_DBG"; exit 1; }
echo "[build] debug  done → $BINARY_DBG"
echo ""

# ── 2. run studies sequentially ──────────────────────────────────────────────
echo "[study] config: $CONFIG"
fail=0

echo ""
echo "[study] ppc sweep (ppc $PPC_MIN–$PPC_MAX, $PPC_REPEATS repeats) …"
$PYTHON "$SCRIPTS/study_ppc.py" "$BINARY" "$CONFIG" \
    --min-ppc "$PPC_MIN" --max-ppc "$PPC_MAX" --repeats "$PPC_REPEATS" \
    || { echo "[FAIL] ppc study"; fail=1; }

echo ""
echo "[study] method comparison ($METHODS) …"
$PYTHON "$SCRIPTS/compare_methods.py" "$BINARY" "$CONFIG" \
    --methods "$METHODS" \
    || { echo "[FAIL] method compare"; fail=1; }

echo ""
echo "[study] solver iteration study ($SOLVERS) — using debug binary …"
$PYTHON "$SCRIPTS/study_solver.py" "$BINARY_DBG" "$CONFIG" \
    --solvers "$SOLVERS" \
    || { echo "[FAIL] solver study"; fail=1; }

echo ""
echo "── outputs ──────────────────────────────────────────────────────────"

BASE=$(python3 -c "import json; cfg=json.load(open('$CONFIG')); print(cfg.get('folder','results'))")
for f in \
    "${BASE}/study_ppc/ke_vs_ppc.png" \
    "${BASE}/compare_methods/ke_comparison.png" \
    "${BASE}/study_solver/iters_vs_timestep.png"; do
    [[ -f "$f" ]] && echo "  ✓ $f" || echo "  ✗ $f"
done

exit $fail
