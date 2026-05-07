#!/usr/bin/env python3
"""
study_ppc.py — sweep ppcx/ppcy from 1 to 10, run N_REPEAT simulations per count
to build a kinetic energy curve with error bars (variation comes from the
random initial particle placement in each cell).

Usage:
    python scripts/study_ppc.py <binary> <config.json> [options]
    python scripts/study_ppc.py <config.json> --no-run [options]

Positional:
    binary            Path to simulation binary (e.g. ./build/bin/PIC)
    config            Path to base JSON config file

Options:
    --min-ppc INT     Minimum ppc value (default: 1)
    --max-ppc INT     Maximum ppc value (default: 10)
    --repeats INT     Runs per ppc value   (default: 5)
    --out DIR         Output directory for plots (default: results/study_ppc)
    --no-run          Skip simulations, plot from existing results

Each run gets its own output folder: <base_folder>/ppc<N>_rep<R>/
The KE is sampled at the final timestep of each run.
"""

import argparse
import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from run import find_binary, parse_pvd, kinetic_energy_from_vtp, kinetic_energy_from_vti


def run_simulation(binary: Path, config: dict, out_folder: str,
                   ppc: int, rep: int, tmp_dir: Path) -> Path:
    cfg = copy.deepcopy(config)
    cfg["ppcx"] = ppc
    cfg["ppcy"] = ppc
    cfg["folder"] = out_folder
    cfg["write_particles"] = True   # need VTP for KE
    cfg["write_u"] = False
    cfg["write_v"] = False
    cfg["write_p"] = False
    cfg["write_div"] = False
    cfg["write_norm_velocity"] = False
    cfg["write_smoke"] = False

    tmp_cfg = tmp_dir / f"ppc{ppc}_rep{rep}.json"
    with open(tmp_cfg, "w") as f:
        json.dump(cfg, f)

    result = subprocess.run([str(binary), str(tmp_cfg)], capture_output=True)
    if result.returncode != 0:
        print(f"  [FAIL] ppc={ppc} rep={rep}: {result.stderr.decode()[:200]}")
    return Path(out_folder)


def extract_final_ke(folder: Path, nx: int, ny: int,
                     density: float, dx: float, dy: float) -> float:
    # Prefer particles.pvd (VTP) over label.pvd (VTI) — label has no velocity data
    pvd_files = sorted(folder.glob("*.pvd"),
                       key=lambda p: (0 if "particle" in p.stem else 1, p.name))
    if not pvd_files:
        return float("nan")
    entries = parse_pvd(pvd_files[0])
    if not entries:
        return float("nan")
    t, fpath = entries[-1]
    if not fpath.exists():
        return float("nan")
    ext = fpath.suffix.lower()
    if ext == ".vtp":
        return kinetic_energy_from_vtp(fpath, density, dx, dy)
    elif ext == ".vti":
        return kinetic_energy_from_vti(fpath, "normVelocity",
                                        nx, ny, density, dx, dy)
    return float("nan")


def main():
    parser = argparse.ArgumentParser(
        usage="%(prog)s <binary> <config.json> [options]\n"
              "       %(prog)s <config.json> --no-run [options]"
    )
    parser.add_argument("binary", help="Path to simulation binary, or config path when --no-run")
    parser.add_argument("config", nargs="?", default=None,
                        help="Path to base JSON config (omit when passing config as first arg with --no-run)")
    parser.add_argument("--min-ppc", type=int, default=1)
    parser.add_argument("--max-ppc", type=int, default=10)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--out", default=None)
    parser.add_argument("--no-run", action="store_true")
    args = parser.parse_args()

    # Allow: script config.json --no-run  (binary omitted)
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

    nx       = int(base_cfg.get("nx", 64))
    ny       = int(base_cfg.get("ny", 64))
    dx       = float(base_cfg.get("dx", 1.0))
    dy       = float(base_cfg.get("dy", 1.0))
    density  = float(base_cfg.get("density", 1000.0))
    base_folder = Path(base_cfg.get("folder", "results")).resolve()

    out_dir = Path(args.out).resolve() if args.out else base_folder / "study_ppc"
    out_dir.mkdir(parents=True, exist_ok=True)

    binary = binary_path  # None when --no-run with single positional arg

    ppc_values = list(range(args.min_ppc, args.max_ppc + 1))
    ke_mean    = []
    ke_std     = []

    for ppc in ppc_values:
        ke_runs = []
        for rep in range(args.repeats):
            run_folder = out_dir / f"ppc{ppc}_rep{rep}"
            run_folder.mkdir(parents=True, exist_ok=True)

            if not args.no_run:
                print(f"  running ppc={ppc} rep={rep} …", flush=True)
                run_simulation(binary, base_cfg, str(run_folder),
                               ppc, rep, out_dir)

            ke = extract_final_ke(run_folder, nx, ny, density, dx, dy)
            ke_runs.append(ke)
            print(f"  ppc={ppc} rep={rep}  KE={ke:.4g}")

        ke_runs = [v for v in ke_runs if not np.isnan(v)]
        ke_mean.append(np.mean(ke_runs) if ke_runs else float("nan"))
        ke_std.append(np.std(ke_runs)  if ke_runs else float("nan"))

    # ── plot ─────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.errorbar(ppc_values, ke_mean, yerr=ke_std,
                fmt="o-", capsize=5, lw=1.5, label="mean ± std")
    ax.set_xlabel("ppcx = ppcy  (particles per cell per direction)")
    ax.set_ylabel("Final kinetic energy (proportional, a.u.)")
    ax.set_title(f"Kinetic energy vs PPC  (N={args.repeats} runs per point)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    out_png = out_dir / "ke_vs_ppc.png"
    fig.savefig(str(out_png), dpi=150)
    plt.close(fig)
    print(f"\n[done] plot saved → {out_png}")

    # save CSV
    csv_path = out_dir / "ke_vs_ppc.csv"
    with open(csv_path, "w") as f:
        f.write("ppc,ke_mean,ke_std\n")
        for p, m, s in zip(ppc_values, ke_mean, ke_std):
            f.write(f"{p},{m},{s}\n")
    print(f"[done] data  saved → {csv_path}")


if __name__ == "__main__":
    main()
