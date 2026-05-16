#!/usr/bin/env python3
"""
Usage:
    python particler.py <file.pvd> --title ... [options]

All options are in the CONFIG dict below.
"""

# ══════════════════════════════════════════════════════════════════════════════
# CONFIG  –  edit here
# ══════════════════════════════════════════════════════════════════════════════

CONFIG = dict(

    # ── Output ────────────────────────────────────────────────────────────────
    out = "cross.mp4",
    fps     = 30,      # Frames per second
    width   = 1280,    # Output width  in pixels
    height  = 800,     # Output height in pixels

    # ── Domain (None → auto from data) ───────────────────────────────────────
    xlim    = None,    # (xmin, xmax) or None
    ylim    = None,    # (ymin, ymax) or None

    # ── Rendering mode ────────────────────────────────────────────────────────
    mode    = "speed", # "speed" = mean velocity per pixel | "density" = log count

    # ── Colormap ──────────────────────────────────────────────────────────────
    cmap    = "viridis",  # any matplotlib colormap
    vmin    = None,    # color-scale min (None → auto)
    vmax    = None,    # color-scale max (None → auto)
    vmax_percentile = 99.5,  # robust auto upper bound; ignores rare spikes

    # ── Colorbar ──────────────────────────────────────────────────────────────
    colorbar_label  = r"$||\boldsymbol{v}||$",  # label on the colorbar
    colorbar_unit   = "",           # unit appended in brackets, e.g. "m/s"
    colorbar_position = "left",     # "left" or "top"
    colorbar_size   = 0.12,         # fraction of body width/height used by the panel
    colorbar_dpi    = 120,
    show_colorbar   = True,   # False → no colorbar panel; frame uses full height

    # ── Background ────────────────────────────────────────────────────────────
    bg_color = (255, 255, 255),  # white background (R, G, B)

    # ── Padding ───────────────────────────────────────────────────────────────
    # Extra background space around the particle area (pixels).
    padding_top    = 0,
    padding_bottom = 15,

    # ── Top header ────────────────────────────────────────────────────────────
    # Clean title/status band above the animation. Use --title "PIC", "FLIP",
    # or "APIC" to label comparison videos.
    show_title     = True,
    header_height  = 58,      # pixels; set 0 to disable the reserved band
    header_bg      = (255, 255, 255),
    header_line    = (226, 229, 234),
    title_fontsize = 20,
    title_color    = "#111827",
    title_weight   = "bold",

    # ── Particle size ─────────────────────────────────────────────────────────
    # Each particle is drawn as a square of (2*particle_radius+1) pixels.
    # 0 = 1 pixel (single dot), 1 = 3×3, 2 = 5×5, 3 = 7×7, etc.
    particle_radius = 2,  # integer ≥ 0

    # ── Time-step label ───────────────────────────────────────────────────────
    show_timestep     = True,
    # The time index is drawn in the top header, aligned right when a title is present.
    timestep_label    = "Time step",
    timestep_fontsize = 11,
    timestep_color    = "#4b5563",
    timestep_weight   = "normal",

    # ── Title ─────────────────────────────────────────────────────────────────
    title   = None,   # None → pvd file stem

    # ── Parallelism ───────────────────────────────────────────────────────────
    workers  = None,  # None → all CPUs
    prefetch = None,  # None → workers*2
    sample   = 1,     # keep 1-in-N frames

)

# ══════════════════════════════════════════════════════════════════════════════

import argparse
import os
import struct
import subprocess
import sys
import zlib
import xml.etree.ElementTree as ET
from concurrent.futures import ProcessPoolExecutor, Future
from collections import deque
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-picm")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator

try:
    from tqdm import tqdm
except ImportError:
    class tqdm:
        def __init__(self, iterable=None, total=None, desc="", **_kwargs):
            self.iterable = iterable
            self.total = total if total is not None else len(iterable or [])
            self.count = 0
            self.desc = desc

        def __iter__(self):
            for item in self.iterable:
                yield item

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            print()

        def update(self, n=1):
            self.count += n
            if self.total:
                print(f"\r{self.desc} {self.count}/{self.total}", end="", flush=True)

        def set_postfix(self, **_kwargs):
            return None

        def close(self):
            print()

        @staticmethod
        def write(message):
            print(message)


def _cmap(name):
    return matplotlib.colormaps[name]


def parse_pvd(pvd_path: Path):
    tree = ET.parse(pvd_path)
    base = pvd_path.parent
    return [base / ds.get("file")
            for ds in tree.getroot().iter("DataSet") if ds.get("file")]


def _extract_xml(raw):
    idx = raw.find(b"<AppendedData")
    if idx == -1:
        return raw
    us = raw.find(b"_", idx)
    return raw[:us] + b"\n  </AppendedData>\n</VTKFile>"


def _bin_start(raw):
    m = raw.find(b"  _")
    return (m + 3) if m != -1 else (raw.find(b"_") + 1)


def _decode(raw, bstart, offset, compressed, dtype):
    chunk = raw[bstart + offset:]
    try:
        if compressed:
            _, _, _, csz = struct.unpack_from("<IIII", chunk)
            return np.frombuffer(zlib.decompress(chunk[16: 16 + csz]), dtype=dtype)
        else:
            (nb,) = struct.unpack_from("<I", chunk)
            return np.frombuffer(chunk[4: 4 + nb], dtype=dtype)
    except Exception as e:
        tqdm.write(f"  [WARN] decode @ offset {offset}: {e}")
        return None


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
        return np.empty(0, np.float32), np.empty(0, np.float32), np.empty(0, np.float32)

    bs = _bin_start(raw)

    pts_da = root.find(".//Points/DataArray")
    if pts_da is None:
        return None, None, None
    dt  = np.float32 if pts_da.get("type", "Float32") == "Float32" else np.float64
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


def _rasterise(x, y, speed, W, H, xlim, ylim, cmap_name, vmin, vmax, mode, bg_color, particle_radius):
    """
    Pure accumulate → mean → colormap. No blur, no blending, no post-processing.
    Pixels with at least one particle get the colormap colour.
    Pixels with no particle keep the background colour.
    """
    img = np.full((H, W, 3), bg_color, dtype=np.uint8)
    if x is None or len(x) == 0:
        return img.tobytes()

    cm   = _cmap(cmap_name)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)

    # Physical coordinates → pixel indices (row 0 = top)
    px = np.round((x - xlim[0]) / (xlim[1] - xlim[0]) * (W - 1)).astype(np.int32)
    py = np.round((ylim[1] - y) / (ylim[1] - ylim[0]) * (H - 1)).astype(np.int32)

    valid = (px >= 0) & (px < W) & (py >= 0) & (py < H) & np.isfinite(speed)
    px, py, speed = px[valid], py[valid], speed[valid]
    if len(px) == 0:
        return img.tobytes()

    count  = np.zeros((H, W), dtype=np.float32)
    values = np.zeros((H, W), dtype=np.float32)

    if particle_radius <= 0:
        # Single-pixel particles (fastest path)
        np.add.at(count,  (py, px), 1.0)
        np.add.at(values, (py, px), speed)
    else:
        # Square kernel: each particle covers a (2r+1)×(2r+1) block of pixels
        r = int(particle_radius)
        for dy in range(-r, r + 1):
            for dx in range(-r, r + 1):
                qx = np.clip(px + dx, 0, W - 1)
                qy = np.clip(py + dy, 0, H - 1)
                np.add.at(count,  (qy, qx), 1.0)
                np.add.at(values, (qy, qx), speed)

    occupied = count > 0

    if mode == "density":
        field = np.zeros((H, W), dtype=np.float32)
        field[occupied] = np.log1p(count[occupied])
        mx = field.max()
        if mx > 0:
            field /= mx
        rgba = cm(field)
    else:
        mean_spd = np.zeros((H, W), dtype=np.float32)
        mean_spd[occupied] = values[occupied] / count[occupied]
        rgba = cm(norm(mean_spd))

    img[occupied] = (rgba[..., :3][occupied] * 255).astype(np.uint8)
    return img.tobytes()


def _worker(args):
    path, W, H, xlim, ylim, cmap_name, vmin, vmax, mode, bg_color, particle_radius = args
    x, y, spd = _load_vtp(Path(path))
    return _rasterise(x, y, spd, W, H, xlim, ylim, cmap_name, vmin, vmax, mode, bg_color, particle_radius)


def _fit_limits_to_canvas(xlim, ylim, W, H):
    xmin, xmax = map(float, xlim)
    ymin, ymax = map(float, ylim)
    xspan, yspan = xmax - xmin, ymax - ymin
    if xspan <= 0 or yspan <= 0:
        return xlim, ylim
    xmid, ymid = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)
    if xspan / yspan > W / H:
        h2 = xspan / (W / H) / 2
        ymin, ymax = ymid - h2, ymid + h2
    else:
        w2 = yspan * (W / H) / 2
        xmin, xmax = xmid - w2, xmid + w2
    return (xmin, xmax), (ymin, ymax)


def _scan(paths, n_scan=12, vmax_percentile=99.5):
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
    xa = np.concatenate(xs)
    ya = np.concatenate(ys)
    sa = np.concatenate(ss)
    sa = sa[np.isfinite(sa)]
    if sa.size == 0:
        return ((float(xa.min()), float(xa.max())),
                (float(ya.min()), float(ya.max())),
                0.0, 1.0)
    vmin = 0.0 if float(sa.min()) >= 0.0 else float(np.nanpercentile(sa, 0.5))
    pct = float(np.clip(vmax_percentile, 50.0, 100.0))
    vmax = float(np.nanpercentile(sa, pct))
    if abs(vmax - vmin) < 1e-30:
        vmax = float(np.nanmax(sa))
    if abs(vmax - vmin) < 1e-30:
        vmax = vmin + 1.0
    return ((float(xa.min()), float(xa.max())),
            (float(ya.min()), float(ya.max())),
            vmin, vmax)


def _make_colorbar_panel(W, H, dpi, title, cmap_name, vmin, vmax, label, unit, orientation):
    full_label = f"{label} [{unit}]" if unit else label
    cm   = _cmap(cmap_name)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)

    fig = plt.figure(figsize=(W / dpi, H / dpi), dpi=dpi)
    fig.patch.set_facecolor("white")
    if title and orientation == "horizontal":
        fig.text(0.5, 0.82, title, ha="center", va="center",
                 fontsize=14, color="black", fontweight="bold")

    if orientation == "vertical":
        cax = fig.add_axes([0.30, 0.12, 0.20, 0.76])
    else:
        # Keep enough room below the bar so tick labels do not get clipped.
        cax = fig.add_axes([0.12, 0.42, 0.76, 0.22])

    sm   = plt.cm.ScalarMappable(cmap=cm, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation=orientation)
    cbar.set_label(full_label, color="black", fontsize=11, labelpad=6)
    cbar.outline.set_edgecolor("black")
    cbar.outline.set_linewidth(0.8)
    tick_axis = "y" if orientation == "vertical" else "x"
    cbar.ax.tick_params(axis=tick_axis, colors="black", labelsize=10, length=4, width=0.8)
    locator_axis = cbar.ax.yaxis if orientation == "vertical" else cbar.ax.xaxis
    locator_axis.set_major_locator(MaxNLocator(nbins=6))

    fig.canvas.draw()
    rgba  = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
    rgba  = rgba.reshape(fig.canvas.get_width_height()[::-1] + (4,))
    panel = rgba[..., :3].copy()
    plt.close(fig)

    if panel.shape[:2] != (H, W):
        try:
            from PIL import Image
            panel = np.array(Image.fromarray(panel).resize((W, H), Image.LANCZOS))
        except ImportError:
            fixed = np.full((H, W, 3), 255, np.uint8)
            h, w = min(panel.shape[0], H), min(panel.shape[1], W)
            fixed[:h, :w] = panel[:h, :w]
            panel = fixed
    return panel


def _rgb255(color):
    if isinstance(color, str):
        return tuple(int(round(c * 255)) for c in mcolors.to_rgb(color))
    vals = tuple(color)
    if max(vals) <= 1:
        return tuple(int(round(c * 255)) for c in vals)
    return tuple(int(c) for c in vals)


def _make_top_header(W, H, title, step_idx, cfg):
    bg = _rgb255(cfg["header_bg"])
    line = _rgb255(cfg["header_line"])
    header = np.full((H, W, 3), bg, dtype=np.uint8)
    if H <= 0:
        return header

    title = title if cfg["show_title"] else ""
    time_text = ""
    if cfg["show_timestep"]:
        time_text = f"{cfg['timestep_label']} {step_idx:04d}"

    if not title and not time_text:
        header[-1:, :, :] = line
        return header

    dpi = 120
    fig = plt.figure(figsize=(W / dpi, H / dpi), dpi=dpi)
    fig.patch.set_facecolor(np.array(bg) / 255.0)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    if title:
        ax.text(
            0.5, 0.56, title,
            transform=ax.transAxes,
            fontsize=cfg["title_fontsize"],
            color=cfg["title_color"],
            fontweight=cfg["title_weight"],
            va="center", ha="center",
        )

    if time_text:
        x = 0.975 if title else 0.5
        ha = "right" if title else "center"
        ax.text(
            x, 0.54, time_text,
            transform=ax.transAxes,
            fontsize=cfg["timestep_fontsize"],
            color=cfg["timestep_color"],
            fontweight=cfg["timestep_weight"],
            va="center", ha=ha,
        )

    fig.canvas.draw()
    rgba = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
    rgba = rgba.reshape(fig.canvas.get_width_height()[::-1] + (4,))
    plt.close(fig)

    rendered = rgba[..., :3].copy()
    if rendered.shape[:2] != (H, W):
        try:
            from PIL import Image
            rendered = np.array(Image.fromarray(rendered).resize((W, H), Image.LANCZOS))
        except ImportError:
            fixed = header.copy()
            h, w = min(rendered.shape[0], H), min(rendered.shape[1], W)
            fixed[:h, :w] = rendered[:h, :w]
            rendered = fixed

    rendered[-1:, :, :] = line
    return rendered.astype(np.uint8)


def _open_ffmpeg(out_path, W, H, fps):
    return subprocess.Popen([
        "ffmpeg", "-y",
        "-f", "rawvideo", "-vcodec", "rawvideo",
        "-s", f"{W}x{H}", "-pix_fmt", "rgb24",
        "-framerate", str(fps), "-loglevel", "error",
        "-i", "pipe:",
        "-vcodec", "libx264", "-preset", "slow",
        "-crf", "16", "-pix_fmt", "yuv420p",
        "-movflags", "+faststart",
        str(out_path),
    ], stdin=subprocess.PIPE)


def build_mp4(vtp_paths, out_path, cfg):
    W = cfg["width"]  + cfg["width"]  % 2
    H = cfg["height"] + cfg["height"] % 2
    pad_t = max(0, int(cfg["padding_top"]))
    pad_b = max(0, int(cfg["padding_bottom"]))
    header_h = max(0, int(cfg["header_height"]))
    if not cfg["show_title"] and not cfg["show_timestep"]:
        header_h = 0
    header_h += header_h % 2
    body_H = H - header_h
    if body_H < 120:
        sys.exit("header_height too large — increase height or reduce header_height")

    print("┌─ Step 1/3  Scanning frames for domain & scale")
    xlim_s, ylim_s, vmin_s, vmax_s = _scan(vtp_paths, vmax_percentile=cfg["vmax_percentile"])
    xlim = cfg["xlim"] or xlim_s
    ylim = cfg["ylim"] or ylim_s
    vmin = cfg["vmin"] if cfg["vmin"] is not None else vmin_s
    vmax = cfg["vmax"] if cfg["vmax"] is not None else vmax_s
    if abs(vmax - vmin) < 1e-30:
        vmax = vmin + 1.0
    print(f"│   x={xlim}  y={ylim}  speed=[{vmin:.3g}, {vmax:.3g}]")

    if cfg["show_colorbar"]:
        print("├─ Step 2/3  Building colorbar panel")
        colorbar_position = cfg["colorbar_position"]
        if colorbar_position == "left":
            panel_w = max(110, int(W * cfg["colorbar_size"]))
            panel_w += panel_w % 2
            frame_W = W - panel_w
            if frame_W < 100:
                sys.exit("colorbar_size too large — increase width or reduce colorbar_size")
            frame_H = body_H
            inner_H = frame_H - pad_t - pad_b
            if inner_H < 10:
                sys.exit("vertical padding too large — increase height or reduce padding")
            total_H = header_h + frame_H
            panel = _make_colorbar_panel(
                panel_w, frame_H, cfg["colorbar_dpi"], None,
                cfg["cmap"], vmin, vmax,
                cfg["colorbar_label"], cfg["colorbar_unit"], "vertical",
            )
            xlim, ylim = _fit_limits_to_canvas(xlim, ylim, frame_W, inner_H)
            template = np.full((total_H, W, 3), cfg["bg_color"], dtype=np.uint8)
            template[header_h: header_h + frame_H, :panel_w] = panel
            print(
                f"│   canvas {W}×{total_H}  header_h={header_h}  colorbar_w={panel_w}  "
                f"padding_top={pad_t}  padding_bottom={pad_b}  x={xlim}  y={ylim}"
            )
        elif colorbar_position == "top":
            panel_h = max(100, int(body_H * cfg["colorbar_size"]))
            panel_h += panel_h % 2
            frame_W = W
            frame_H = body_H - panel_h
            if frame_H < 100:
                sys.exit("colorbar_size too large — increase height or reduce colorbar_size")
            frame_H += frame_H % 2
            inner_H = frame_H - pad_t - pad_b
            if inner_H < 10:
                sys.exit("vertical padding too large — increase height or reduce padding")
            total_H = header_h + panel_h + frame_H
            panel = _make_colorbar_panel(
                W, panel_h, cfg["colorbar_dpi"], cfg["title"] if header_h == 0 else None,
                cfg["cmap"], vmin, vmax,
                cfg["colorbar_label"], cfg["colorbar_unit"], "horizontal",
            )
            xlim, ylim = _fit_limits_to_canvas(xlim, ylim, frame_W, inner_H)
            template = np.full((total_H, W, 3), cfg["bg_color"], dtype=np.uint8)
            template[header_h: header_h + panel_h] = panel
            print(
                f"│   canvas {W}×{total_H}  header_h={header_h}  colorbar_h={panel_h}  "
                f"padding_top={pad_t}  padding_bottom={pad_b}  x={xlim}  y={ylim}"
            )
        else:
            sys.exit('colorbar_position must be "left" or "top"')
    else:
        print("├─ Step 2/3  No colorbar panel (show_colorbar=False)")
        panel_h = 0
        panel_w = 0
        frame_W = W
        frame_H = body_H
        total_H = header_h + frame_H
        inner_H = frame_H - pad_t - pad_b
        if inner_H < 10:
            sys.exit("vertical padding too large — increase height or reduce padding")
        xlim, ylim = _fit_limits_to_canvas(xlim, ylim, W, inner_H)
        template = np.full((total_H, W, 3), cfg["bg_color"], dtype=np.uint8)
        print(
            f"│   canvas {W}×{total_H}  (no panel)  header_h={header_h}  "
            f"padding_top={pad_t}  padding_bottom={pad_b}  x={xlim}  y={ylim}"
        )

    n_workers = cfg["workers"] or os.cpu_count()
    prefetch  = cfg["prefetch"] or n_workers * 2
    print(f"└─ Step 3/3  Encoding {len(vtp_paths)} frames ({n_workers} workers)")

    worker_args = [
        (str(p), frame_W, inner_H, xlim, ylim,
         cfg["cmap"], vmin, vmax, cfg["mode"], cfg["bg_color"], cfg["particle_radius"])
        for p in vtp_paths
    ]

    proc  = _open_ffmpeg(out_path, W, total_H, cfg["fps"])
    outer = tqdm(total=len(worker_args), desc="  encoding", unit="frame", position=0)
    inner = tqdm(total=prefetch, desc="  queued  ", unit="frame", position=1, leave=False)

    def write_frame(step_idx, rgb_bytes):
        frame = template.copy()
        if header_h > 0:
            frame[:header_h] = _make_top_header(W, header_h, cfg["title"], step_idx, cfg)

        row0 = header_h + (panel_h if cfg["show_colorbar"] and cfg["colorbar_position"] == "top" else 0) + pad_t
        col0 = panel_w if cfg["show_colorbar"] and cfg["colorbar_position"] == "left" else 0
        frame[row0: row0 + inner_H, col0: col0 + frame_W] = (
            np.frombuffer(rgb_bytes, dtype=np.uint8).reshape(inner_H, frame_W, 3)
        )

        proc.stdin.write(frame.tobytes())

    if n_workers == 1:
        inner.close()
        for step_idx, args in enumerate(worker_args):
            write_frame(step_idx, _worker(args))
            outer.update(1)
        outer.close()
        proc.stdin.close()
        proc.wait()
        if proc.returncode != 0:
            sys.exit(f"ffmpeg exited with code {proc.returncode}")
        print(f"Saved: {out_path}")
        return

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        pending: deque[tuple[int, Future]] = deque()

        def submit_next(idx):
            if idx < len(worker_args):
                pending.append((idx, pool.submit(_worker, worker_args[idx])))
                inner.update(1)

        for i in range(min(prefetch, len(worker_args))):
            submit_next(i)
        next_submit = prefetch

        for step_idx in range(len(worker_args)):
            _, fut = pending.popleft()
            inner.update(-1)
            rgb_bytes = fut.result()
            submit_next(next_submit); next_submit += 1

            write_frame(step_idx, rgb_bytes)
            outer.set_postfix(queue=len(pending))
            outer.update(1)

    outer.close(); inner.close()
    proc.stdin.close(); proc.wait()
    if proc.returncode != 0:
        sys.exit(f"ffmpeg exited with code {proc.returncode}")
    print(f"Saved: {out_path}")


def main():
    ap = argparse.ArgumentParser(description="PVD+VTP particles → MP4")
    ap.add_argument("pvd")
    ap.add_argument("--out",      default=None)
    ap.add_argument("--fps",      type=int,   default=None)
    ap.add_argument("--cmap",     default=None)
    ap.add_argument("--vmin",     type=float, default=None)
    ap.add_argument("--vmax",     type=float, default=None)
    ap.add_argument("--width",    type=int,   default=None)
    ap.add_argument("--height",   type=int,   default=None)
    ap.add_argument("--title",    default=None)
    ap.add_argument("--no-title", action="store_true",
                    help="Hide the top title while keeping the time index if enabled")
    ap.add_argument("--header-height", type=int, default=None,
                    help="Top header height in pixels (0 disables the reserved band)")
    ap.add_argument("--title-fontsize", type=int, default=None)
    ap.add_argument("--no-timestep", action="store_true",
                    help="Hide the time index in the top header")
    ap.add_argument("--timestep-label", default=None)
    ap.add_argument("--mode",     choices=["speed", "density"], default=None)
    ap.add_argument("--workers",  type=int,   default=None)
    ap.add_argument("--prefetch", type=int,   default=None)
    ap.add_argument("--sample",   type=int,   default=None)
    ap.add_argument("--xlim",           type=float, nargs=2, default=None)
    ap.add_argument("--ylim",           type=float, nargs=2, default=None)
    ap.add_argument("--padding-top",    type=int,   default=None,
                    help="Background padding above the particle area, in pixels")
    ap.add_argument("--padding-bottom", type=int,   default=None,
                    help="Background padding below the particle area, in pixels")
    ap.add_argument("--particle-radius", type=int,   default=None,
                    help="Particle display radius in pixels (0=1px, 1=3x3, 2=5x5, …)")
    args = ap.parse_args()

    pvd_path = Path(args.pvd)
    if not pvd_path.exists():
        sys.exit(f"Not found: {pvd_path}")

    cfg = dict(CONFIG)
    if args.out      is not None: cfg["out"]      = args.out
    if args.fps      is not None: cfg["fps"]      = args.fps
    if args.cmap     is not None: cfg["cmap"]     = args.cmap
    if args.vmin     is not None: cfg["vmin"]     = args.vmin
    if args.vmax     is not None: cfg["vmax"]     = args.vmax
    if args.width    is not None: cfg["width"]    = args.width
    if args.height   is not None: cfg["height"]   = args.height
    if args.title    is not None: cfg["title"]    = args.title
    if args.no_title: cfg["show_title"] = False
    if args.header_height is not None: cfg["header_height"] = args.header_height
    if args.title_fontsize is not None: cfg["title_fontsize"] = args.title_fontsize
    if args.no_timestep: cfg["show_timestep"] = False
    if args.timestep_label is not None: cfg["timestep_label"] = args.timestep_label
    if args.mode     is not None: cfg["mode"]     = args.mode
    if args.workers  is not None: cfg["workers"]  = args.workers
    if args.prefetch is not None: cfg["prefetch"] = args.prefetch
    if args.sample   is not None: cfg["sample"]   = args.sample
    if args.xlim            is not None: cfg["xlim"]            = tuple(args.xlim)
    if args.ylim            is not None: cfg["ylim"]            = tuple(args.ylim)
    if args.padding_top     is not None: cfg["padding_top"]     = args.padding_top
    if args.padding_bottom  is not None: cfg["padding_bottom"]  = args.padding_bottom
    if args.particle_radius is not None: cfg["particle_radius"] = args.particle_radius

    if cfg["title"] is None:
        cfg["title"] = pvd_path.stem
    out_path = Path(cfg["out"]) if cfg["out"] else pvd_path.with_suffix(".mp4")

    all_paths = [p for p in parse_pvd(pvd_path) if p.exists()]
    if cfg["sample"] and cfg["sample"] > 1:
        all_paths = all_paths[::cfg["sample"]]
    if not all_paths:
        sys.exit("No VTP files found.")

    cfg["workers"] = max(1, cfg["workers"] or os.cpu_count())
    print(f"particler → {out_path.name}  |  {len(all_paths)} frames  |  mode={cfg['mode']}  |  cmap={cfg['cmap']}")
    build_mp4(all_paths, out_path, cfg)


if __name__ == "__main__":
    main()
