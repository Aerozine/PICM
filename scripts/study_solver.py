#!/usr/bin/env python3
"""
study_solver.py — compare pressure-solver iteration counts across solver types.

Requires a DEBUG build (compiled without -DNDEBUG) so that [DEBUG] lines are
emitted to stderr.  Release builds are silent and will produce empty plots.

Usage:
    python scripts/study_solver.py <binary> <config.json> [options]
    python scripts/study_solver.py <config.json> --no-run [options]

Positional:
    binary              Path to simulation binary (e.g. ./build/bin/PIC)
    config              Path to base JSON config file

Options:
    --solvers LIST      Comma-separated solver types to compare
                        (default: red_black_gauss_seidel,cg,miccg0)
    --repeats INT       Runs per solver for averaging (default: 3)
    --out DIR           Output directory (default: results/study_solver)
    --no-run            Skip simulations, re-plot from saved CSVs
    --nt INT            Override nt (number of timesteps, default: from config)

Each run folder: <out_dir>/<solver>_rep<R>/
Outputs: iters_vs_timestep.png, iters_vs_timestep.csv
"""

import argparse
import copy
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(it, **kw):
        return it
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from run import find_binary

CONVERGED_RE = re.compile(r"\[DEBUG\]\s+(\S+)\s+converged in\s+(\d+)\s+iter")
MAX_ITER_RE  = re.compile(r"\[DEBUG\]\s+(\S+).*reached maxIters\s*=\s*(\d+)")

ALL_SOLVERS = ["red_black_gauss_seidel", "cg", "miccg0", "jacobi", "gauss_seidel"]


def run_and_capture(binary: Path, cfg: dict, tmp_cfg: Path, out_folder: str,
                    solver_type: str, rep: int) -> list[int]:
    """Run simulation and return list of iteration counts (one per pressure solve)."""
    c = copy.deepcopy(cfg)
    c["solver"] = copy.deepcopy(cfg.get("solver", {}))
    c["solver"]["type"] = solver_type
    c["folder"] = out_folder
    c["write_particles"] = False
    c["write_u"] = False
    c["write_v"] = False
    c["write_p"] = False
    c["write_div"] = False
    c["write_norm_velocity"] = False
    c["write_smoke"] = False

    with open(tmp_cfg, "w") as f:
        json.dump(c, f)

    result = subprocess.run(
        [str(binary), str(tmp_cfg)],
        capture_output=True, text=True
    )

    counts = []
    for line in (result.stdout + result.stderr).splitlines():
        m = CONVERGED_RE.search(line)
        if m:
            counts.append(int(m.group(2)))
            continue
        m = MAX_ITER_RE.search(line)
        if m:
            counts.append(int(m.group(2)))

    if result.returncode != 0 and not counts:
        print(f"  [WARN] {solver_type} rep={rep} exit={result.returncode} — "
              f"no [DEBUG] lines found (release build?)")

    return counts


def load_or_run(binary, base_cfg, solver_type, repeats, out_dir, tmp_dir):
    all_counts = []  # list of lists (one list per repeat)
    for rep in range(repeats):
        run_folder = out_dir / f"{solver_type}_rep{rep}"
        run_folder.mkdir(parents=True, exist_ok=True)
        cache = run_folder / "iters.json"

        # Load from cache — but discard empty caches (left by release-build runs)
        if cache.exists():
            with open(cache) as f:
                counts = json.load(f)
            if counts:
                tqdm.write(f"  {solver_type} rep={rep}  {len(counts)} steps (cached)")
                all_counts.append(counts)
                continue
            cache.unlink()

        tmp_cfg = tmp_dir / f"{solver_type}_rep{rep}.json"
        tqdm.write(f"  running {solver_type} rep={rep} …")
        counts = run_and_capture(binary, base_cfg, tmp_cfg,
                                  str(run_folder), solver_type, rep)
        mean_str = f"{np.mean(counts):.1f}" if counts else "none"
        tqdm.write(f"  {solver_type} rep={rep}  {len(counts)} pressure solves  mean={mean_str}")
        if counts:  # only cache non-empty results
            with open(cache, "w") as f:
                json.dump(counts, f)

        all_counts.append(counts)
    return all_counts


def align_and_average(all_counts):
    """Align repeat runs to the shortest length, return mean and std arrays."""
    if not all_counts:
        return np.array([]), np.array([])
    min_len = min(len(c) for c in all_counts)
    arr = np.array([c[:min_len] for c in all_counts], dtype=float)
    return arr.mean(axis=0), arr.std(axis=0)


def main():
    parser = argparse.ArgumentParser(
        usage="%(prog)s <binary> <config.json> [options]\n"
              "       %(prog)s <config.json> --no-run [options]"
    )
    parser.add_argument("binary", help="Path to simulation binary, or config path when --no-run")
    parser.add_argument("config", nargs="?", default=None)
    parser.add_argument("--solvers", default="red_black_gauss_seidel,cg,miccg0")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--out", default=None)
    parser.add_argument("--no-run", action="store_true")
    parser.add_argument("--nt", type=int, default=None)
    args = parser.parse_args()

    if args.config is None:
        if not args.no_run:
            parser.error("binary and config are required unless --no-run is set")
        config_path = Path(args.binary).resolve()
        binary_path = None
    else:
        binary_path = Path(args.binary).resolve()
        config_path = Path(args.config).resolve()

    root = config_path.parents[2]
    with open(config_path) as f:
        base_cfg = json.load(f)

    if args.nt is not None:
        base_cfg["nt"] = args.nt

    base_folder = Path(base_cfg.get("folder", "results")).resolve()
    out_dir = Path(args.out).resolve() if args.out else base_folder / "study_solver"
    out_dir.mkdir(parents=True, exist_ok=True)

    solver_list = [s.strip() for s in args.solvers.split(",")]

    binary = binary_path
    if binary is not None:
        print(f"[binary] {binary}")

    results = {}
    for solver in tqdm(solver_list, unit="solver", desc="solvers"):
        all_counts = load_or_run(binary, base_cfg, solver,
                                  args.repeats, out_dir, out_dir)
        mean, std = align_and_average(all_counts)
        results[solver] = (mean, std)

    # ── CSV ──────────────────────────────────────────────────────────────────
    csv_path = out_dir / "iters_vs_timestep.csv"
    header = "timestep," + ",".join(
        f"{s}_mean,{s}_std" for s in solver_list
    )
    max_len = max((len(v[0]) for v in results.values()), default=0)
    with open(csv_path, "w") as f:
        f.write(header + "\n")
        for t in range(max_len):
            row = [str(t)]
            for s in solver_list:
                m, sd = results[s]
                row.append(f"{m[t]:.2f}" if t < len(m) else "")
                row.append(f"{sd[t]:.2f}" if t < len(sd) else "")
            f.write(",".join(row) + "\n")
    print(f"\n[done] CSV  → {csv_path}")

    # ── plot ─────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(11, 5))
    colors = plt.cm.tab10.colors
    for idx, solver in enumerate(solver_list):
        mean, std = results[solver]
        if len(mean) == 0:
            print(f"  [WARN] {solver}: no data to plot")
            continue
        xs = np.arange(len(mean))
        c = colors[idx % len(colors)]
        ax.plot(xs, mean, lw=1.5, color=c, label=solver)
        ax.fill_between(xs, mean - std, mean + std, alpha=0.2, color=c)

    ax.set_xlabel("Timestep index")
    ax.set_ylabel("Pressure-solver iterations")
    ax.set_title(f"Iterations per timestep by solver  (N={args.repeats} repeats, shaded = ±std)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    out_png = out_dir / "iters_vs_timestep.png"
    fig.savefig(str(out_png), dpi=150)
    plt.close(fig)
    print(f"[done] plot → {out_png}")

    # ── summary table ─────────────────────────────────────────────────────────
    print("\n── Summary (mean iterations per step) ──")
    for solver in solver_list:
        mean, _ = results[solver]
        if len(mean):
            print(f"  {solver:<28} avg={mean.mean():.1f}  max={mean.max():.0f}")
        else:
            print(f"  {solver:<28} no data")


if __name__ == "__main__":
    main()
