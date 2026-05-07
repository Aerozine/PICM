#!/usr/bin/env python3
"""
compare_methods.py — run the same scenario with multiple advection methods and
compare kinetic energy time series on a single plot.

Usage:
    python scripts/compare_methods.py <binary> <config.json> [options]
    python scripts/compare_methods.py <config.json> --no-run [options]

Positional:
    binary              Path to simulation binary (e.g. ./build/bin/PIC)
    config              Path to base JSON config file

Options:
    --methods LIST      Comma-separated method names to compare
                        (default: vanilla_pic,flip,apic,semilagrangian)
    --out DIR           Output directory (default: results/compare_methods)
    --no-run            Skip simulations, re-plot from existing PVD/VTP files
    --nt INT            Override nt (number of timesteps)

Output: ke_comparison.png, ke_comparison.csv
"""

import argparse
import copy
import json
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
from run import (find_binary, parse_pvd,
                 kinetic_energy_from_vtp, kinetic_energy_from_vti)

ALL_METHODS = ["vanilla_pic", "flip", "apic", "semilagrangian", "gpic"]


def run_method(binary: Path, base_cfg: dict, method: str,
               run_folder: Path, tmp_cfg: Path) -> bool:
    cfg = copy.deepcopy(base_cfg)
    cfg["method"]           = method
    cfg["folder"]           = str(run_folder)
    cfg["write_particles"]  = True   # VTP for KE
    cfg["write_u"]          = False
    cfg["write_v"]          = False
    cfg["write_p"]          = False
    cfg["write_div"]        = False
    cfg["write_norm_velocity"] = False
    cfg["write_smoke"]      = False

    with open(tmp_cfg, "w") as f:
        json.dump(cfg, f)

    result = subprocess.run([str(binary), str(tmp_cfg)], capture_output=True)
    if result.returncode != 0:
        err = result.stderr.decode(errors="replace")[:300]
        print(f"  [FAIL] method={method}: {err}")
        return False
    return True


def extract_ke_series(folder: Path, nx: int, ny: int,
                      density: float, dx: float, dy: float) -> tuple[list, list]:
    # Prefer particles.pvd (VTP) — label.pvd has no velocity data
    pvd_files = sorted(folder.glob("*.pvd"),
                       key=lambda p: (0 if "particle" in p.stem else 1, p.name))
    if not pvd_files:
        return [], []
    entries = parse_pvd(pvd_files[0])
    if not entries:
        return [], []

    times, ke_vals = [], []
    for t, fpath in entries:
        if not fpath.exists():
            continue
        ext = fpath.suffix.lower()
        if ext == ".vtp":
            ke = kinetic_energy_from_vtp(fpath, density, dx, dy)
        elif ext == ".vti":
            ke = kinetic_energy_from_vti(fpath, "normVelocity",
                                          nx, ny, density, dx, dy)
        else:
            continue
        times.append(t)
        ke_vals.append(ke)
    return times, ke_vals


def main():
    parser = argparse.ArgumentParser(
        usage="%(prog)s <binary> <config.json> [options]\n"
              "       %(prog)s <config.json> --no-run [options]"
    )
    parser.add_argument("binary", help="Path to simulation binary, or config path when --no-run")
    parser.add_argument("config", nargs="?", default=None)
    parser.add_argument("--methods",
                        default="vanilla_pic,flip,apic,semilagrangian")
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

    nx      = int(base_cfg.get("nx", 64))
    ny      = int(base_cfg.get("ny", 64))
    dx      = float(base_cfg.get("dx", 1.0))
    dy      = float(base_cfg.get("dy", 1.0))
    density = float(base_cfg.get("density", 1000.0))

    base_folder = Path(base_cfg.get("folder", "results")).resolve()
    out_dir = Path(args.out).resolve() if args.out else base_folder / "compare_methods"
    out_dir.mkdir(parents=True, exist_ok=True)

    method_list = [m.strip() for m in args.methods.split(",")]

    binary = binary_path
    if binary is not None:
        print(f"[binary] {binary}")

    results = {}  # method → (times, ke_vals)

    for method in tqdm(method_list, unit="method", desc="methods"):
        run_folder = out_dir / method
        run_folder.mkdir(parents=True, exist_ok=True)

        if not args.no_run:
            tmp_cfg = out_dir / f"_tmp_{method}.json"
            ok = run_method(binary, base_cfg, method, run_folder, tmp_cfg)
            if not ok:
                tqdm.write(f"  [SKIP] {method} failed")
                continue
            tmp_cfg.unlink(missing_ok=True)

        times, ke_vals = extract_ke_series(run_folder, nx, ny, density, dx, dy)
        if not times:
            tqdm.write(f"  [WARN] {method}: no output found in {run_folder}")
            continue
        results[method] = (times, ke_vals)
        print(f"  {method}: {len(times)} timesteps, "
              f"final KE={ke_vals[-1]:.4g}")

    if not results:
        print("[ERROR] No data to plot.")
        sys.exit(1)

    # ── CSV ──────────────────────────────────────────────────────────────────
    csv_path = out_dir / "ke_comparison.csv"
    with open(csv_path, "w") as f:
        f.write("method,time,ke\n")
        for method, (times, ke_vals) in results.items():
            for t, ke in zip(times, ke_vals):
                f.write(f"{method},{t},{ke}\n")
    print(f"\n[done] CSV  → {csv_path}")

    # ── plot ─────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 5))
    colors = plt.cm.tab10.colors
    for idx, (method, (times, ke_vals)) in enumerate(results.items()):
        ax.plot(times, ke_vals, lw=1.5, color=colors[idx % len(colors)],
                label=method)

    ax.set_xlabel("Simulation time [s]")
    ax.set_ylabel("Kinetic energy (proportional, a.u.)")
    ax.set_title("Kinetic energy over time — method comparison")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    out_png = out_dir / "ke_comparison.png"
    fig.savefig(str(out_png), dpi=150)
    plt.close(fig)
    print(f"[done] plot → {out_png}")


if __name__ == "__main__":
    main()
