#!/usr/bin/env python3
"""
sim2mp4.py  –  PVD+VTP (particules) ou PVD+VTI (champ scalaire) → MP4
Détecte automatiquement le type de fichier, ou forçage avec --type.

Usage:
    python sim2mp4.py <file.pvd> [options]

Options communes:
    --out  PATH             Fichier MP4 de sortie         (défaut: <pvd_stem>.mp4)
    --fps  INT              Images par seconde             (défaut: 30)
    --cmap NAME             Colormap matplotlib            (défaut: viridis)
    --vmin FLOAT            Minimum de l'échelle couleur  (défaut: auto)
    --vmax FLOAT            Maximum de l'échelle couleur  (défaut: auto)
    --dpi  INT              DPI figure                     (défaut: 150)
    --title TEXT            Titre affiché                 (défaut: stem du pvd)
    --type particles|field  Forcer le mode               (défaut: auto-détection)
    --sample INT            Garder 1 frame sur N          (défaut: 1)

Options particules (VTP):
    --width INT             Largeur vidéo en pixels        (défaut: 1280)
    --height INT            Hauteur vidéo en pixels        (défaut: 720)
    --mode speed|density    Couleur par vitesse ou densité (défaut: speed)
    --vmax-percentile FLOAT Percentile pour vmax auto      (défaut: 99.5)
    --particle-size FLOAT   Rayon particule en pixels      (défaut: 1)
    --particle-alpha FLOAT  Opacité particules 0..1        (défaut: 1.0)
    --particle-style NAME   point|square|disk|soft         (défaut: point)
    --gamma FLOAT           Gamma densité                  (défaut: 0.75)
    --margin FLOAT          Marge auto-fit                 (défaut: 0.08)
    --panel-height INT      Hauteur bande colorbar px      (défaut: auto)
    --panel-position top|bottom                            (défaut: bottom)
    --no-title              Cacher le titre
    --xlim FLOAT FLOAT      Limites physiques x
    --ylim FLOAT FLOAT      Limites physiques y
    --quality INT           CRF x264                      (défaut: 14)
    --preset NAME           Preset x264                   (défaut: slower)
    --pixfmt NAME           yuv420p|yuv444p               (défaut: yuv420p)
    --workers INT           Nombre de workers             (défaut: cpu_count)

Options champ scalaire (VTI):
    --global-scale          Calculer vmin/vmax sur toutes les frames
"""

import argparse
import os
import struct
import subprocess
import sys
import zlib
import xml.etree.ElementTree as ET
from collections import deque
from concurrent.futures import Future, ProcessPoolExecutor
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

    class tqdm:
        """Minimal tqdm fallback."""
        def __init__(self, iterable=None, total=None, desc="", **kw):
            self._it = iterable
            self._n = 0
            self._total = total or (len(iterable) if iterable is not None else 0)
            self._desc = desc

        def __iter__(self):
            for item in self._it:
                self._n += 1
                print(f"\r  {self._desc} {self._n}/{self._total}", end="", flush=True)
                yield item
            print()

        def __enter__(self): return self
        def __exit__(self, *a): print()
        def update(self, n=1): self._n += n
        def set_postfix(self, **kw): pass
        def close(self): print()
        @staticmethod
        def write(msg): print(msg)


# ═══════════════════════════════════════════════════════════════════════════════
#  COMMUN – PVD & XML
# ═══════════════════════════════════════════════════════════════════════════════

def parse_pvd(pvd_path: Path) -> list[Path]:
    tree = ET.parse(pvd_path)
    base = pvd_path.parent
    return [base / ds.get("file")
            for ds in tree.getroot().iter("DataSet") if ds.get("file")]


def _extract_xml(raw: bytes) -> bytes:
    idx = raw.find(b"<AppendedData")
    if idx == -1:
        return raw
    us = raw.find(b"_", idx)
    if us == -1:
        return raw
    return raw[:us] + b"\n  </AppendedData>\n</VTKFile>"


def _bin_start(raw: bytes) -> int:
    m = raw.find(b"  _")
    return (m + 3) if m != -1 else (raw.find(b"_") + 1)


def _decode(raw: bytes, bstart: int, offset: int,
            compressed: bool, dtype) -> np.ndarray | None:
    chunk = raw[bstart + offset:]
    try:
        if compressed:
            _, _, _, csz = struct.unpack_from("<IIII", chunk)
            return np.frombuffer(zlib.decompress(chunk[16: 16 + csz]), dtype=dtype)
        (nb,) = struct.unpack_from("<I", chunk)
        return np.frombuffer(chunk[4: 4 + nb], dtype=dtype)
    except Exception as e:
        tqdm.write(f"  [WARN] decode @ offset {offset}: {e}")
        return None


def detect_type(paths: list[Path]) -> str:
    """Retourne 'particles' (VTP) ou 'field' (VTI) selon l'extension."""
    for p in paths:
        suf = p.suffix.lower()
        if suf == ".vtp":
            return "particles"
        if suf == ".vti":
            return "field"
    return "field"


# ═══════════════════════════════════════════════════════════════════════════════
#  MODE PARTICULES (VTP)
# ═══════════════════════════════════════════════════════════════════════════════

def _cmap(name: str):
    return matplotlib.colormaps[name]


def _load_vtp(path: Path):
    raw = path.read_bytes()
    try:
        root = ET.fromstring(_extract_xml(raw))
    except ET.ParseError as e:
        tqdm.write(f"  [WARN] {path.name}: {e}")
        return None, None, None

    compressed = "compressor" in root.attrib
    piece = root.find(".//Piece")
    if piece is None:
        return None, None, None

    n = int(piece.get("NumberOfPoints", 0))
    if n == 0:
        return (np.empty(0, np.float32),) * 3

    bs = _bin_start(raw)
    pts_da = root.find(".//Points/DataArray")
    if pts_da is None:
        return None, None, None

    dt = np.float32 if pts_da.get("type", "Float32") == "Float32" else np.float64
    pts = _decode(raw, bs, int(pts_da.get("offset", 0)), compressed, dt)
    if pts is None or pts.size < n * 3:
        return None, None, None
    pts = pts[:n * 3].reshape(n, 3)

    spd = None
    for da in root.iter("DataArray"):
        if da.get("Name") == "normVelocity":
            dt2 = np.float32 if da.get("type", "Float32") == "Float32" else np.float64
            spd = _decode(raw, bs, int(da.get("offset", 0)), compressed, dt2)
            if spd is not None:
                spd = spd[:n].astype(np.float32)
            break
    if spd is None:
        spd = np.zeros(n, np.float32)

    return pts[:, 0].astype(np.float32), pts[:, 1].astype(np.float32), spd


def _scan_vtp(paths, n_scan=12, vmax_percentile=99.5):
    if not paths:
        return (0.0, 1.0), (0.0, 1.0), 0.0, 1.0

    indices = sorted(set(
        [0] + list(range(0, len(paths), max(1, len(paths) // n_scan))) + [len(paths) - 1]
    ))
    xs, ys, ss = [], [], []

    with tqdm(indices, desc="  scanning", unit="frame", leave=False) as pbar:
        for i in pbar:
            x, y, s = _load_vtp(paths[i])
            if x is not None and len(x) > 0:
                xs.append(x); ys.append(y); ss.append(s)

    if not xs:
        return (0.0, 1.0), (0.0, 1.0), 0.0, 1.0

    xa = np.concatenate(xs); ya = np.concatenate(ys)
    sa = np.concatenate(ss); sa = sa[np.isfinite(sa)]

    if sa.size == 0:
        vmin, vmax = 0.0, 1.0
    else:
        vmin = float(np.nanmin(sa))
        pct = float(np.clip(vmax_percentile, 50.0, 100.0))
        vmax = float(np.nanpercentile(sa, pct))
        if abs(vmax - vmin) < 1e-30:
            vmax = float(np.nanmax(sa))
        if abs(vmax - vmin) < 1e-30:
            vmax = vmin + 1.0

    return ((float(np.nanmin(xa)), float(np.nanmax(xa))),
            (float(np.nanmin(ya)), float(np.nanmax(ya))),
            vmin, vmax)


def _expand_limits(lo, hi, fallback_span=1.0):
    if not np.isfinite(lo) or not np.isfinite(hi):
        return (0.0, fallback_span)
    if abs(hi - lo) < 1e-30:
        c = 0.5 * (lo + hi); half = 0.5 * fallback_span
        return (c - half, c + half)
    return (lo, hi)


def _fit_limits_to_canvas(xlim, ylim, W, H, margin=0.08):
    xmin, xmax = _expand_limits(float(xlim[0]), float(xlim[1]))
    ymin, ymax = _expand_limits(float(ylim[0]), float(ylim[1]))
    xmid = 0.5 * (xmin + xmax); ymid = 0.5 * (ymin + ymax)
    xspan = (xmax - xmin) * (1.0 + 2.0 * max(0.0, margin))
    yspan = (ymax - ymin) * (1.0 + 2.0 * max(0.0, margin))
    canvas_aspect = W / max(1, H)
    data_aspect = xspan / max(yspan, 1e-30)
    if data_aspect > canvas_aspect:
        yspan = xspan / canvas_aspect
    else:
        xspan = yspan * canvas_aspect
    return ((xmid - 0.5 * xspan, xmid + 0.5 * xspan),
            (ymid - 0.5 * yspan, ymid + 0.5 * yspan))


def _particle_offsets(radius, style):
    style = str(style or "point").lower()
    r = max(0.0, float(radius))
    if style == "point" or r <= 0.5:
        return [(0, 0, 1.0)]
    ri = int(np.ceil(r))
    offsets = []
    for dy in range(-ri, ri + 1):
        for dx in range(-ri, ri + 1):
            if style == "square":
                offsets.append((dx, dy, 1.0)); continue
            d = (dx * dx + dy * dy) ** 0.5
            if style == "disk":
                if d <= r:
                    offsets.append((dx, dy, 1.0))
            else:
                if d <= r + 0.5:
                    w = float(np.clip(r + 0.5 - d, 0.0, 1.0))
                    if w > 0:
                        offsets.append((dx, dy, w))
    return offsets or [(0, 0, 1.0)]


def _rasterise(x, y, speed, W, H, xlim, ylim, cmap_name, vmin, vmax,
               mode, particle_radius, particle_alpha, particle_style, gamma):
    bg = np.full((H, W, 3), 255, dtype=np.float32)
    if x is None or len(x) == 0:
        return bg.astype(np.uint8).tobytes()

    cm = _cmap(cmap_name)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    px = ((x - xlim[0]) / (xlim[1] - xlim[0]) * (W - 1)).astype(np.int32)
    py = ((ylim[1] - y) / (ylim[1] - ylim[0]) * (H - 1)).astype(np.int32)
    mask = (px >= 0) & (px < W) & (py >= 0) & (py < H) & np.isfinite(speed)
    px, py, speed = px[mask], py[mask], speed[mask]
    if len(px) == 0:
        return bg.astype(np.uint8).tobytes()

    counts = np.zeros((H, W), dtype=np.float32)
    weights = np.zeros((H, W), dtype=np.float32)
    for dx, dy, w in _particle_offsets(particle_radius, particle_style):
        xx = px + dx; yy = py + dy
        m = (xx >= 0) & (xx < W) & (yy >= 0) & (yy < H)
        if not np.any(m):
            continue
        np.add.at(counts, (yy[m], xx[m]), w)
        if mode == "speed":
            np.add.at(weights, (yy[m], xx[m]), speed[m] * w)

    occupied = counts > 0
    if not np.any(occupied):
        return bg.astype(np.uint8).tobytes()

    if mode == "density":
        data = np.log1p(counts)
        data = data / (float(data.max()) + 1e-10)
        data = np.power(data, max(1e-6, gamma))
        rgb = (cm(data)[..., :3] * 255.0).astype(np.float32)
    else:
        mean_spd = np.zeros_like(weights)
        mean_spd[occupied] = weights[occupied] / counts[occupied]
        rgb = (cm(norm(mean_spd))[..., :3] * 255.0).astype(np.float32)

    alpha_base = float(np.clip(particle_alpha, 0.0, 1.0))
    if alpha_base >= 0.999:
        bg[occupied] = rgb[occupied]
    else:
        alpha = alpha_base * np.clip(counts / max(1.0, float(counts.max())), 0.50, 1.0)
        alpha = alpha[..., None]
        bg[occupied] = (1.0 - alpha[occupied]) * bg[occupied] + alpha[occupied] * rgb[occupied]

    return np.clip(bg, 0, 255).astype(np.uint8).tobytes()


def _worker_vtp(args):
    (path, W, H, xlim, ylim, cmap_name, vmin, vmax, mode,
     particle_radius, particle_alpha, particle_style, gamma) = args
    x, y, spd = _load_vtp(Path(path))
    return _rasterise(x, y, spd, W, H, xlim, ylim, cmap_name, vmin, vmax,
                      mode, particle_radius, particle_alpha, particle_style, gamma)


def _make_panel(W, panel_h, dpi, title, cmap_name, vmin, vmax, label,
                show_title=True):
    panel_h = int(max(80, panel_h))
    fig = plt.figure(figsize=(W / dpi, panel_h / dpi), dpi=dpi)
    fig.patch.set_facecolor("white")
    left, width = 0.13, 0.74
    bar_bottom = 0.24 if show_title else 0.34
    cax = fig.add_axes([left, bar_bottom, width, 0.28])
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    sm = plt.cm.ScalarMappable(cmap=_cmap(cmap_name), norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_label(label, color="black", fontsize=12, labelpad=8)
    cbar.outline.set_edgecolor("black")
    cbar.outline.set_linewidth(0.8)
    cbar.ax.tick_params(axis="x", colors="black", labelsize=10, length=4, width=0.8)
    cbar.ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    plt.setp(cbar.ax.xaxis.get_ticklabels(), rotation=0, ha="center")
    if show_title and title:
        fig.text(0.5, 0.78, title, ha="center", va="center",
                 fontsize=14, color="black", fontweight="bold")
    fig.canvas.draw()
    rgba_buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
    rgba_buf = rgba_buf.reshape(fig.canvas.get_width_height()[::-1] + (4,))
    panel = rgba_buf[..., :3].copy()
    plt.close(fig)
    return panel


def _open_ffmpeg(out_path, W, H, fps, quality, preset, pixfmt):
    return subprocess.Popen([
        "ffmpeg", "-y",
        "-f", "rawvideo", "-vcodec", "rawvideo",
        "-s", f"{W}x{H}", "-pix_fmt", "rgb24",
        "-framerate", str(fps), "-loglevel", "error",
        "-i", "pipe:",
        "-vcodec", "libx264", "-preset", str(preset),
        "-crf", str(quality), "-tune", "animation",
        "-pix_fmt", str(pixfmt),
        "-movflags", "+faststart",
        str(out_path),
    ], stdin=subprocess.PIPE)


def build_mp4_particles(vtp_paths, out_path, fps, cmap_name,
                        vmin, vmax, vmax_percentile, W, H, dpi, title, mode,
                        xlim, ylim, margin, n_workers, prefetch,
                        particle_size, particle_alpha, particle_style,
                        panel_height, panel_position, no_title,
                        quality, preset, pixfmt, gamma):

    W += W % 2; H += H % 2

    print("┌─ Étape 1/3  Scan des frames (domaine & échelle)")
    xlim_s, ylim_s, vmin_s, vmax_s = _scan_vtp(vtp_paths, vmax_percentile=vmax_percentile)

    xlim_auto = xlim is None; ylim_auto = ylim is None
    if xlim is None: xlim = xlim_s
    if ylim is None: ylim = ylim_s
    if vmin is None: vmin = vmin_s
    if vmax is None: vmax = vmax_s
    if abs(vmax - vmin) < 1e-30: vmax = vmin + 1.0

    print("├─ Étape 2/3  Génération du panneau colorbar")
    panel_h = int(panel_height) if panel_height else max(110, int(H * 0.14))
    panel_h += panel_h % 2
    frame_H = H - panel_h
    if frame_H < 100:
        sys.exit("Panneau trop grand. Augmentez --height ou réduisez --panel-height.")
    frame_H += frame_H % 2
    total_H = frame_H + panel_h

    if xlim_auto or ylim_auto:
        fit_xlim, fit_ylim = _fit_limits_to_canvas(xlim, ylim, W, frame_H, margin=margin)
        if xlim_auto: xlim = fit_xlim
        if ylim_auto: ylim = fit_ylim

    particle_radius = max(0.0, float(particle_size)) if particle_size is not None else 1.0
    label = "mean ‖v‖" if mode == "speed" else "log(1 + particle count)"
    panel = _make_panel(W, panel_h, dpi, title, cmap_name, vmin, vmax, label,
                        show_title=not no_title)
    panel = panel[:panel_h, :W, :]

    empty_frame = np.full((frame_H, W, 3), 255, np.uint8)
    if panel_position == "top":
        template = np.concatenate([panel, empty_frame], axis=0)
        data_slice = slice(panel_h, panel_h + frame_H)
    else:
        template = np.concatenate([empty_frame, panel], axis=0)
        data_slice = slice(0, frame_H)

    print(f"│   canvas {W}×{total_H}  frame={W}×{frame_H}  panel={panel_h}px")
    print(f"│   x={xlim}  y={ylim}  vitesse=[{vmin:.3g}, {vmax:.3g}]")
    print(f"│   rayon={particle_radius:g}px  style={particle_style}  alpha={particle_alpha:g}  cmap={cmap_name}")
    print(f"└─ Étape 3/3  Encodage de {len(vtp_paths)} frames ({n_workers} workers)")

    worker_args = [
        (str(p), W, frame_H, xlim, ylim, cmap_name, vmin, vmax, mode,
         particle_radius, particle_alpha, particle_style, gamma)
        for p in vtp_paths
    ]

    proc = _open_ffmpeg(out_path, W, total_H, fps, quality, preset, pixfmt)
    if proc.stdin is None:
        sys.exit("Impossible d'ouvrir stdin ffmpeg.")

    outer = tqdm(total=len(worker_args), desc="  encodage", unit="frame")
    inner = tqdm(total=prefetch, desc="  en attente", unit="frame", leave=False)

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        pending: deque[tuple[int, Future]] = deque()

        def submit_next(idx):
            if idx < len(worker_args):
                fut = pool.submit(_worker_vtp, worker_args[idx])
                pending.append((idx, fut))
                inner.update(1)

        for i in range(min(prefetch, len(worker_args))):
            submit_next(i)
        next_submit = prefetch

        for _ in range(len(worker_args)):
            _, fut = pending.popleft()
            inner.update(-1)
            rgb_bytes = fut.result()
            submit_next(next_submit); next_submit += 1
            frame = template.copy()
            frame[data_slice] = np.frombuffer(rgb_bytes, dtype=np.uint8).reshape(frame_H, W, 3)
            proc.stdin.write(frame.tobytes())
            outer.set_postfix(queue=len(pending)); outer.update(1)

    outer.close(); inner.close()
    proc.stdin.close(); proc.wait()
    if proc.returncode != 0:
        sys.exit(f"ffmpeg a quitté avec le code {proc.returncode}")
    print(f"Sauvegardé: {out_path}")


# ═══════════════════════════════════════════════════════════════════════════════
#  MODE CHAMP SCALAIRE (VTI) – pipeline matplotlib → ffmpeg
# ═══════════════════════════════════════════════════════════════════════════════

def load_vti(vti_path: Path) -> np.ndarray | None:
    raw = vti_path.read_bytes()
    try:
        root = ET.fromstring(_extract_xml(raw))
    except ET.ParseError as exc:
        print(f"  [WARN] XML parse error in {vti_path.name}: {exc}")
        return None

    compressed = "compressor" in root.attrib
    img = root.find(".//ImageData")
    if img is None:
        print(f"  [WARN] Pas de <ImageData> dans {vti_path.name}"); return None

    extent = list(map(int, img.get("WholeExtent", "0 1 0 1 0 0").split()))
    nx, ny = extent[1], extent[3]

    da = root.find(".//DataArray")
    if da is None:
        print(f"  [WARN] Pas de <DataArray> dans {vti_path.name}"); return None

    dtype = np.float32 if da.get("type", "Float32") == "Float32" else np.float64
    bs = _bin_start(raw)
    arr = _decode(raw, bs, int(da.get("offset", 0)), compressed, dtype)
    if arr is None:
        return None

    expected = nx * ny
    if arr.size < expected:
        print(f"  [WARN] {vti_path.name}: attendu {expected} valeurs, obtenu {arr.size}")
        return None

    return arr[:expected].reshape(ny, nx)


def build_mp4_field(vti_paths: list[Path], out_path: Path,
                    fps: int, cmap: str,
                    vmin: float | None, vmax: float | None,
                    global_scale: bool,
                    title: str, dpi: int, quality: int = 18,
                    preset: str = "slow", pixfmt: str = "yuv420p"):

    print(f"Chargement de {len(vti_paths)} frames VTI…")
    frames = []
    for i, p in enumerate(vti_paths):
        if not p.exists():
            print(f"  [SKIP] {p.name} – fichier introuvable"); continue
        arr = load_vti(p)
        if arr is not None:
            frames.append(arr)
        if (i + 1) % 20 == 0 or (i + 1) == len(vti_paths):
            print(f"  Chargé {i+1}/{len(vti_paths)}", end="\r", flush=True)
    print(f"\n{len(frames)} frames valides.")

    if not frames:
        print("Aucune frame à encoder."); return

    if vmin is None: vmin = float(min(f.min() for f in frames))
    if vmax is None: vmax = float(max(f.max() for f in frames))
    if abs(vmax - vmin) < 1e-30: vmax = vmin + 1.0

    ny, nx = frames[0].shape
    aspect = nx / ny
    fig_w = max(6.0, 8.0 * aspect)
    fig, ax = plt.subplots(figsize=(fig_w, 5), dpi=dpi)
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")
    ax.set_title(title, color="white", fontsize=11)
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor("white")

    im = ax.imshow(frames[0], origin="lower", cmap=cmap,
                   vmin=vmin, vmax=vmax, aspect="auto", interpolation="nearest")
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white")
    step_text = ax.text(0.01, 0.97, "", transform=ax.transAxes,
                        color="white", fontsize=9, va="top")

    fig.canvas.draw()
    W_fig, H_fig = fig.canvas.get_width_height()
    W_fig += W_fig % 2; H_fig += H_fig % 2

    proc = _open_ffmpeg(out_path, W_fig, H_fig, fps, quality, preset, pixfmt)
    if proc.stdin is None:
        sys.exit("Impossible d'ouvrir stdin ffmpeg.")

    print(f"Encodage {len(frames)} frames → {out_path.name}  ({W_fig}×{H_fig})")
    for i, arr in enumerate(tqdm(frames, desc="  encodage", unit="frame")):
        im.set_data(arr)
        step_text.set_text(f"frame {i}/{len(frames)-1}")
        if not global_scale:
            lo, hi = float(arr.min()), float(arr.max())
            if abs(hi - lo) > 1e-30:
                im.set_clim(lo, hi)
        fig.canvas.draw()
        buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
        buf = buf.reshape(fig.canvas.get_width_height()[::-1] + (4,))
        rgb = buf[:H_fig, :W_fig, :3].copy()
        proc.stdin.write(rgb.tobytes())

    plt.close(fig)
    proc.stdin.close(); proc.wait()
    if proc.returncode != 0:
        sys.exit(f"ffmpeg a quitté avec le code {proc.returncode}")
    print(f"Sauvegardé: {out_path}  ({len(frames)} frames @ {fps} fps)")


# ═══════════════════════════════════════════════════════════════════════════════
#  CLI
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description="PVD+VTP (particules) ou PVD+VTI (champ scalaire) → MP4"
    )
    ap.add_argument("pvd", help="Chemin vers le fichier .pvd")

    # ── options communes ──────────────────────────────────────────────────────
    ap.add_argument("--out",    default=None)
    ap.add_argument("--fps",    type=int,   default=30)
    ap.add_argument("--cmap",   default="viridis")
    ap.add_argument("--vmin",   type=float, default=None)
    ap.add_argument("--vmax",   type=float, default=None)
    ap.add_argument("--dpi",    type=int,   default=150)
    ap.add_argument("--title",  default=None)
    ap.add_argument("--sample", type=int,   default=1,
                    help="Garder 1 frame sur N")
    ap.add_argument("--type", choices=["particles", "field"], default=None,
                    help="Forcer le mode (défaut: auto-détection par extension)")

    # ── options particules (VTP) ──────────────────────────────────────────────
    ap.add_argument("--width",  type=int, default=1280)
    ap.add_argument("--height", type=int, default=720)
    ap.add_argument("--mode",   choices=["speed", "density"], default="speed")
    ap.add_argument("--vmax-percentile", type=float, default=99.5)
    ap.add_argument("--particle-size",   type=float, default=1.0)
    ap.add_argument("--particle-alpha",  type=float, default=1.0)
    ap.add_argument("--particle-style",
                    choices=["point", "square", "disk", "soft"], default="point")
    ap.add_argument("--gamma",  type=float, default=0.75)
    ap.add_argument("--margin", type=float, default=0.08)
    ap.add_argument("--panel-height",    type=int, default=None)
    ap.add_argument("--panel-position",  choices=["top", "bottom"], default="bottom")
    ap.add_argument("--no-title",        action="store_true")
    ap.add_argument("--xlim",  type=float, nargs=2, default=None, metavar=("XMIN", "XMAX"))
    ap.add_argument("--ylim",  type=float, nargs=2, default=None, metavar=("YMIN", "YMAX"))
    ap.add_argument("--quality", type=int, default=14)
    ap.add_argument("--preset",  default="slower")
    ap.add_argument("--pixfmt",  choices=["yuv444p", "yuv420p"], default="yuv420p")
    ap.add_argument("--workers", type=int, default=os.cpu_count())

    # ── options champ scalaire (VTI) ──────────────────────────────────────────
    ap.add_argument("--global-scale", action="store_true",
                    help="Calculer vmin/vmax sur toutes les frames (VTI)")

    args = ap.parse_args()

    pvd_path = Path(args.pvd)
    if not pvd_path.exists():
        sys.exit(f"Fichier introuvable: {pvd_path}")

    out_path = Path(args.out) if args.out else pvd_path.with_suffix(".mp4")
    title = args.title or pvd_path.stem

    all_paths = parse_pvd(pvd_path)
    if not all_paths:
        sys.exit("Aucune entrée DataSet trouvée dans le PVD.")

    if args.sample > 1:
        all_paths = all_paths[::args.sample]

    # Auto-détection du type
    sim_type = args.type or detect_type(all_paths)
    print(f"sim2mp4  →  {out_path.name}  [mode: {sim_type}]")

    if sim_type == "particles":
        existing = [p for p in all_paths if p.exists()]
        if not existing:
            sys.exit("Aucun fichier VTP trouvé.")
        print(f"  {len(existing)} frames · {args.workers} workers")
        n_workers = max(1, int(args.workers or 1))
        prefetch = n_workers * 2
        build_mp4_particles(
            existing, out_path,
            fps=args.fps, cmap_name=args.cmap,
            vmin=args.vmin, vmax=args.vmax,
            vmax_percentile=args.vmax_percentile,
            W=args.width, H=args.height, dpi=args.dpi,
            title=title, mode=args.mode,
            xlim=tuple(args.xlim) if args.xlim else None,
            ylim=tuple(args.ylim) if args.ylim else None,
            margin=args.margin, n_workers=n_workers, prefetch=prefetch,
            particle_size=args.particle_size, particle_alpha=args.particle_alpha,
            particle_style=args.particle_style,
            panel_height=args.panel_height, panel_position=args.panel_position,
            no_title=args.no_title,
            quality=args.quality, preset=args.preset, pixfmt=args.pixfmt,
            gamma=args.gamma,
        )

    else:  # field / VTI
        print(f"  {len(all_paths)} frames")
        build_mp4_field(
            all_paths, out_path,
            fps=args.fps, cmap=args.cmap,
            vmin=args.vmin, vmax=args.vmax,
            global_scale=args.global_scale,
            title=title, dpi=args.dpi,
            quality=args.quality, preset=args.preset, pixfmt=args.pixfmt,
        )


if __name__ == "__main__":
    main()