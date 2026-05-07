#!/usr/bin/env python3
"""
run.py — launch a simulation and plot results including kinetic energy.

Usage:
    python scripts/run.py <config.json> [options]

Options:
    --bin PATH      Path to simulation binary (default: auto-detect)
    --no-run        Skip the simulation, only plot existing results
    --fields NAME   Comma-separated VTI field names to plot (default: p,normVelocity)
    --last          Plot only the last timestep (default: all)
    --out DIR       Output directory for plots (default: <folder>/plots)
    --ke            Plot kinetic energy time series (requires VTP particle files)
    --density FLOAT Fluid density kg/m³ used for kinetic energy (default: from config)
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── helpers ─────────────────────────────────────────────────────────────────

def find_binary(root: Path) -> Path:
    candidates = [
        root / "build-gpu" / "bin" / "PIC",
        root / "build" / "bin" / "PIC",
        root / "build-gpu" / "PIC",
        root / "build" / "PIC",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(
        f"Simulation binary not found. Tried: {[str(c) for c in candidates]}"
    )


def load_config(config_path: Path) -> dict:
    with open(config_path) as f:
        return json.load(f)


# ── VTI reader ───────────────────────────────────────────────────────────────

def read_field_2d(vti_path: Path, field: str, nx: int, ny: int):
    try:
        import pyvista as pv
        mesh = pv.read(str(vti_path))
        store = mesh.point_data if field in mesh.point_data else mesh.cell_data
        arr = store[field]
        dims = mesh.dimensions
        return arr.reshape(dims[1] - 1, dims[0] - 1) if arr.size == (dims[0]-1)*(dims[1]-1) \
               else arr.reshape(dims[1], dims[0])
    except ImportError:
        import value as vmod
        return np.array(vmod.extract_field(str(vti_path), field)).reshape(ny, nx)


# ── VTP reader for particles ─────────────────────────────────────────────────

def read_vtp(vtp_path: Path):
    """Returns (x, y, speed) arrays from a particle VTP file."""
    import pyvista as pv
    mesh = pv.read(str(vtp_path))
    pts = np.array(mesh.points)
    x, y = pts[:, 0], pts[:, 1]
    if "normVelocity" in mesh.point_data:
        speed = np.array(mesh.point_data["normVelocity"])
    else:
        speed = np.zeros(len(x))
    return x, y, speed


def kinetic_energy_from_vtp(vtp_path: Path, density: float, dx: float, dy: float) -> float:
    """
    Kinetic energy from a particle VTP file.
    KE = Σ  ½ * (density * dx * dy / n_per_cell) * v²
    Approximated as ½ * density * dx * dy * mean(v²) * n_particles
    when all particles have equal mass.
    """
    import pyvista as pv
    mesh = pv.read(str(vtp_path))
    if mesh.n_points == 0:
        return 0.0
    if "normVelocity" in mesh.point_data:
        speed = np.array(mesh.point_data["normVelocity"])
    else:
        pts = np.array(mesh.points)
        speed = np.zeros(pts.shape[0])
    # Each particle represents density * cell_area / ppc fluid
    # We don't store ppc here so use total: KE = ½ * m_total/N * Σ v²
    # m_total ≈ density * Σ(fluid cells) * dx*dy  — not available here.
    # Best approximation: KE proportional to Σ v² (same mass per particle)
    return float(0.5 * np.sum(speed ** 2))


def kinetic_energy_from_vti(vti_path: Path, field: str, nx: int, ny: int,
                             density: float, dx: float, dy: float) -> float:
    """KE from grid normVelocity: ½ * density * dx*dy * Σ v²"""
    try:
        arr = read_field_2d(vti_path, field, nx, ny)
        return float(0.5 * density * dx * dy * np.sum(arr ** 2))
    except Exception:
        return float("nan")


# ── PVD parser ────────────────────────────────────────────────────────────────

def parse_pvd(pvd_path: Path) -> list:
    import xml.etree.ElementTree as ET
    tree = ET.parse(pvd_path)
    root = tree.getroot()
    entries = []
    for ds in root.iter("DataSet"):
        t = float(ds.get("timestep", 0))
        f = pvd_path.parent / ds.get("file", "")
        entries.append((t, f))
    return sorted(entries, key=lambda x: x[0])


# ── plotting ──────────────────────────────────────────────────────────────────

def plot_field(arr, title, out_path, cmap="viridis", dx=1.0, dy=1.0):
    ny, nx = arr.shape
    fig, ax = plt.subplots(figsize=(8, 8 * ny / nx))
    im = ax.imshow(arr, origin="lower", cmap=cmap,
                   extent=[0, nx * dx, 0, ny * dy], aspect="equal")
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(str(out_path), dpi=120)
    plt.close(fig)


def plot_particles(x, y, xmax, ymax, title, out_path):
    fig, ax = plt.subplots(figsize=(8, 8 * ymax / xmax))
    ax.scatter(x, y, s=0.5, c="steelblue", alpha=0.5, linewidths=0)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(str(out_path), dpi=120)
    plt.close(fig)


def plot_kinetic_energy(times, ke_values, out_path, label=""):
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(times, ke_values, lw=1.5, label=label or "KE")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Kinetic energy [J/m] (proportional)")
    ax.set_title("Kinetic energy over time")
    ax.grid(True, alpha=0.3)
    if label:
        ax.legend()
    fig.tight_layout()
    fig.savefig(str(out_path), dpi=120)
    plt.close(fig)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Run simulation and plot results.")
    parser.add_argument("config", help="Path to JSON config file")
    parser.add_argument("--bin", default=None, help="Path to simulation binary")
    parser.add_argument("--no-run", action="store_true")
    parser.add_argument("--fields", default="p,normVelocity")
    parser.add_argument("--last", action="store_true",
                        help="Plot only the last timestep")
    parser.add_argument("--out", default=None)
    parser.add_argument("--ke", action="store_true",
                        help="Plot kinetic energy time series")
    parser.add_argument("--density", type=float, default=None)
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    root = config_path.parents[2]
    cfg = load_config(config_path)

    nx = int(cfg.get("nx", 64))
    ny = int(cfg.get("ny", 64))
    dx = float(cfg.get("dx", 1.0))
    dy = float(cfg.get("dy", 1.0))
    density = args.density or float(cfg.get("density", 1000.0))
    result_folder = Path(cfg.get("folder", "results")).resolve()
    field_names = [f.strip() for f in args.fields.split(",")]

    if not args.no_run:
        binary = Path(args.bin).resolve() if args.bin else find_binary(root)
        print(f"[run] binary : {binary}")
        print(f"[run] config : {config_path}")
        ret = subprocess.run([str(binary), str(config_path)])
        if ret.returncode != 0:
            print(f"[run] simulation failed (exit {ret.returncode})", file=sys.stderr)
            sys.exit(ret.returncode)
        print("[run] simulation complete")

    out_dir = Path(args.out).resolve() if args.out else result_folder / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    pvd_files = sorted(result_folder.glob("*.pvd"))
    if not pvd_files:
        print(f"[plot] No .pvd files found in {result_folder}")
        return

    print(f"[plot] found {len(pvd_files)} PVD file(s)")

    for pvd_path in pvd_files:
        entries = parse_pvd(pvd_path)
        if not entries:
            continue
        ext = entries[0][1].suffix.lower() if entries else ""
        to_plot = [entries[-1]] if args.last else entries

        # ── kinetic energy time series ──────────────────────────────────────
        if args.ke:
            times, ke_vals = [], []
            for t, fpath in entries:
                if not fpath.exists():
                    continue
                if ext == ".vtp":
                    ke = kinetic_energy_from_vtp(fpath, density, dx, dy)
                elif ext == ".vti" and "normVelocity" in field_names:
                    ke = kinetic_energy_from_vti(fpath, "normVelocity", nx, ny,
                                                  density, dx, dy)
                else:
                    continue
                times.append(t)
                ke_vals.append(ke)
            if times:
                ke_png = out_dir / f"{pvd_path.stem}_kinetic_energy.png"
                plot_kinetic_energy(times, ke_vals, ke_png)
                print(f"[plot] KE saved → {ke_png.name}")

        # ── field / particle frames ─────────────────────────────────────────
        for t, fpath in to_plot:
            if not fpath.exists():
                continue
            tag = f"t={t:.4g}"
            stem = fpath.stem

            if ext == ".vti":
                for field in field_names:
                    try:
                        arr = read_field_2d(fpath, field, nx, ny)
                    except (KeyError, Exception) as e:
                        print(f"[plot] {fpath.name}/{field}: {e}")
                        continue
                    cmap = "RdBu_r" if "p" in field.lower() else "viridis"
                    out_png = out_dir / f"{stem}_{field}.png"
                    plot_field(arr, f"{field}  {tag}", out_png,
                               cmap=cmap, dx=dx, dy=dy)
                    print(f"[plot] saved {out_png.name}")

            elif ext == ".vtp":
                try:
                    x, y, _ = read_vtp(fpath)
                    out_png = out_dir / f"{stem}_particles.png"
                    plot_particles(x, y, nx * dx, ny * dy,
                                   f"particles  {tag}", out_png)
                    print(f"[plot] saved {out_png.name}")
                except Exception as e:
                    print(f"[plot] {fpath.name}: {e}")

    print(f"[plot] all plots saved to {out_dir}")


if __name__ == "__main__":
    main()
