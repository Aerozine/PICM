import struct
import sys
import zlib
import csv
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

# =============================================================================
# CONFIGURATION — modifiez ces valeurs avant de lancer le script
# =============================================================================

PVD_PATH   = Path("../data/SL/von-karman/reference/normVelocity.pvd")
FIELD_NAME = "normVelocity"
OUTPUT_CSV = Path("normVel-sl-vk-ref.csv")

# --- Mode d'extraction ---
# "point" : extrait la valeur en un seul point (i, j) — CSV léger, 1 colonne de valeur
# "field" : extrait tout le champ à chaque pas de temps — CSV large (nx*ny colonnes)
#           nécessaire pour calculer E_inf, E2, etc. dans picPlot.py
MODE = "point"   # "point"  ou  "field"

# Utilisé uniquement si MODE = "point"
I_INDEX = 125     # indice colonne  (axe x)
J_INDEX = 25    # indice ligne    (axe y)

# =============================================================================


def parse_pvd(pvd_path: Path) -> list:
    tree = ET.parse(pvd_path)
    root = tree.getroot()
    base = pvd_path.parent
    entries = []
    for ds in root.iter("DataSet"):
        t = ds.get("timestep", ds.get("time", None))
        f = ds.get("file")
        if f:
            entries.append((float(t) if t is not None else None, base / f))
    return entries


def _extract_xml(raw: bytes) -> bytes:
    tag = b"<AppendedData"
    idx = raw.find(tag)
    if idx == -1:
        return raw
    underscore = raw.find(b"_", idx)
    if underscore == -1:
        return raw
    return raw[:underscore] + b"\n  </AppendedData>\n</VTKFile>"


def _read_appended_data(raw: bytes, compressed: bool, dtype: np.dtype) -> np.ndarray | None:
    marker = raw.find(b"  _")
    data_start = marker + 3 if marker != -1 else raw.find(b"_") + 1
    chunk = raw[data_start:]
    try:
        if compressed:
            _num_blocks, uncomp_size, _last, comp_size = struct.unpack_from("<IIII", chunk, 0)
            payload = zlib.decompress(chunk[16 : 16 + comp_size])
        else:
            (raw_bytes,) = struct.unpack_from("<I", chunk, 0)
            payload = chunk[4 : 4 + raw_bytes]
    except Exception as exc:
        print(f"  [WARN] Binary decode failed: {exc}")
        return None
    return np.frombuffer(payload, dtype=dtype)


def load_vti(vti_path: Path) -> np.ndarray | None:
    """Parse un fichier VTI → tableau 2-D numpy shape (ny, nx), ou None."""
    raw = vti_path.read_bytes()
    try:
        root = ET.fromstring(_extract_xml(raw))
    except ET.ParseError as exc:
        print(f"  [WARN] XML parse error in {vti_path.name}: {exc}")
        return None

    compressed = "compressor" in root.attrib

    img = root.find(".//ImageData")
    if img is None:
        print(f"  [WARN] No <ImageData> in {vti_path.name}")
        return None

    extent = list(map(int, img.get("WholeExtent", "0 1 0 1 0 0").split()))
    nx, ny = extent[1], extent[3]

    da = root.find(".//DataArray")
    if da is None:
        print(f"  [WARN] No <DataArray> in {vti_path.name}")
        return None

    dtype = np.float32 if da.get("type", "Float32") == "Float32" else np.float64
    arr = _read_appended_data(raw, compressed, dtype)
    if arr is None:
        return None

    expected = nx * ny
    if arr.size < expected:
        print(f"  [WARN] {vti_path.name}: expected {expected}, got {arr.size}")
        return None

    return arr[:expected].reshape(ny, nx)


# =============================================================================
# MODE POINT — un seul point (i, j) par timestep
# =============================================================================

def main_point():
    print(f"[MODE POINT]  PVD: {PVD_PATH}")
    print(f"  Point (i={I_INDEX}, j={J_INDEX})  →  arr[j={J_INDEX}, i={I_INDEX}]")
    print(f"  Sortie: {OUTPUT_CSV}\n")

    entries = parse_pvd(PVD_PATH)
    if not entries:
        print("[ERREUR] Aucun fichier VTI trouvé.")
        sys.exit(1)

    rows = []
    for step_idx, (timestep, vti_path) in enumerate(entries):
        print(f"  [{step_idx+1}/{len(entries)}] {vti_path.name}", end="  ")
        arr = load_vti(vti_path)
        if arr is None:
            print("→ ignoré (erreur)")
            continue

        ny, nx = arr.shape
        if J_INDEX >= ny or I_INDEX >= nx:
            print(f"→ ignoré (hors limites : ({ny},{nx}))")
            continue

        value = float(arr[J_INDEX, I_INDEX])
        print(f"→ {value:.6g}")
        rows.append({
            "step":     step_idx,
            "timestep": timestep if timestep is not None else "",
            "file":     vti_path.name,
            f"{FIELD_NAME}[j={J_INDEX},i={I_INDEX}]": value,
        })

    if not rows:
        print("\n[ERREUR] Aucune valeur extraite.")
        sys.exit(1)

    _write_csv(rows, OUTPUT_CSV)


# =============================================================================
# MODE FIELD — tout le champ aplati par timestep
# CSV : step | timestep | file | u[0,0] | u[0,1] | ... | u[ny-1,nx-1]
# Les colonnes sont nommées  {FIELD_NAME}[j=J,i=I]
# =============================================================================

def main_field():
    print(f"[MODE FIELD]  PVD: {PVD_PATH}")
    print(f"  Sortie: {OUTPUT_CSV}\n")

    entries = parse_pvd(PVD_PATH)
    if not entries:
        print("[ERREUR] Aucun fichier VTI trouvé.")
        sys.exit(1)

    rows = []
    field_cols = None   # noms des colonnes u[j=...,i=...], construits au 1er timestep

    for step_idx, (timestep, vti_path) in enumerate(entries):
        print(f"  [{step_idx+1}/{len(entries)}] {vti_path.name}", end="  ")
        arr = load_vti(vti_path)
        if arr is None:
            print("→ ignoré (erreur)")
            continue

        ny, nx = arr.shape

        # Construire les noms de colonnes une seule fois
        if field_cols is None:
            field_cols = [
                f"{FIELD_NAME}[j={j},i={i}]"
                for j in range(ny)
                for i in range(nx)
            ]
            print(f"  (champ {ny}×{nx} = {ny*nx} valeurs / timestep)")
        
        flat = arr.ravel().tolist()
        row = {
            "step":     step_idx,
            "timestep": timestep if timestep is not None else "",
            "file":     vti_path.name,
        }
        row.update(dict(zip(field_cols, flat)))
        rows.append(row)
        print(f"→ ok")

    if not rows:
        print("\n[ERREUR] Aucune valeur extraite.")
        sys.exit(1)

    _write_csv(rows, OUTPUT_CSV)


# =============================================================================
# Utilitaire CSV
# =============================================================================

def _write_csv(rows: list, path: Path):
    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n✓ CSV écrit : {path}  ({len(rows)} lignes, {len(fieldnames)} colonnes)")


# =============================================================================

def main():
    if MODE == "point":
        main_point()
    elif MODE == "field":
        main_field()
    else:
        print(f"[ERREUR] MODE inconnu : '{MODE}'. Choisir 'point' ou 'field'.")
        sys.exit(1)


if __name__ == "__main__":
    main()
