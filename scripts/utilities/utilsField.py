"""
utilsField.py — Plot à partir d'un CSV mode "field"
===================================================

Version avec labels associés directement aux fichiers.

Tu peux fournir :
    labels = {
        "normVel-cfl-1e-2-field.csv": r"$CFL=10^{-2}$",
        "normVel-cfl-1e-1-field.csv": r"$CFL=10^{-1}$",
    }

Les clés acceptées sont :
    - le chemin complet du fichier ;
    - le nom du fichier avec extension ;
    - le stem du fichier sans extension.
"""

import glob
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utilities.utilsStyle import COLORS, style_ax, add_hline, finalize


# ============================================================
#  FORMULAS
# ============================================================

def _field_values(df: pd.DataFrame, field: str) -> np.ndarray:
    """Extrait toutes les colonnes field[j=*,i=*] → tableau (T, N)."""
    cols = [c for c in df.columns if c.startswith(f"{field}[")]
    if not cols:
        raise ValueError(
            f"Aucune colonne '{field}[j=...,i=...]' trouvée.\n"
            f"Colonnes disponibles : {list(df.columns)}"
        )
    return df[cols].to_numpy(dtype=float)


def E_inf_formula(field: str, u_ref: float = 1.0):
    """E_inf(t) = max_i |u_i(t) - u_ref|."""
    def _f(df):
        U = _field_values(df, field)
        return np.max(np.abs(U - u_ref), axis=1)
    return _f


def E2_formula(field: str, u_ref: float = 1.0):
    """E2(t) = sqrt(mean_i((u_i(t) - u_ref)^2))."""
    def _f(df):
        U = _field_values(df, field)
        return np.sqrt(np.mean((U - u_ref) ** 2, axis=1))
    return _f


def point_formula(field: str, i: int, j: int):
    """Extrait u[j,i](t) depuis un CSV field."""
    col = f"{field}[j={j},i={i}]"
    def _f(df):
        if col not in df.columns:
            raise ValueError(f"Colonne '{col}' absente du DataFrame.")
        return df[col].to_numpy(dtype=float)
    return _f


# ============================================================
#  DÉTECTION AUTOMATIQUE DU CHAMP
# ============================================================

def _detect_field_name(df: pd.DataFrame) -> str | None:
    for col in df.columns:
        m = re.match(r"^(.+)\[j=\d+,i=\d+\]$", col)
        if m:
            return m.group(1)
    return None


def _auto_formula(df: pd.DataFrame):
    field = _detect_field_name(df)
    if field is None:
        raise ValueError(
            "Impossible de détecter le champ dans ce CSV.\n"
            "Vérifiez que le fichier est bien en MODE='field'.\n"
            f"Colonnes trouvées : {list(df.columns)}"
        )
    print(f"  [auto] champ '{field}' détecté → E_inf par défaut")
    return E_inf_formula(field)


# ============================================================
#  HELPERS INTERNES
# ============================================================

def _natural_sort_key(path):
    """Tri alphabétique stable par nom de fichier."""
    return os.path.basename(str(path))



def _label_key_candidates(fpath):
    fpath = Path(fpath)
    return {
        str(fpath),
        fpath.as_posix(),
        fpath.name,
        fpath.stem,
    }


def _order_files_from_label_dict(csv_files, labels):
    """
    Si labels est un dictionnaire {fichier: label}, impose l'ordre des courbes
    donné par l'ordre d'insertion du dictionnaire.

    Les clés peuvent être :
        - le chemin complet ;
        - le chemin POSIX ;
        - le nom du fichier avec extension ;
        - le stem du fichier sans extension.
    """
    csv_files = [Path(f) for f in csv_files]

    if not isinstance(labels, dict):
        return sorted(csv_files, key=_natural_sort_key)

    remaining = list(csv_files)
    ordered = []

    for key in labels.keys():
        key = str(key)
        matches = [f for f in remaining if key in _label_key_candidates(f)]

        if len(matches) == 0:
            available = "\n".join(f"  - {f.name}" for f in csv_files)
            raise KeyError(
                f"La clé de légende ne correspond à aucun fichier : {key}\n\n"
                f"Fichiers disponibles :\n{available}"
            )

        if len(matches) > 1:
            found = "\n".join(f"  - {f}" for f in matches)
            raise KeyError(
                f"La clé de légende est ambiguë : {key}\n\n"
                f"Elle correspond à plusieurs fichiers :\n{found}"
            )

        match = matches[0]
        ordered.append(match)
        remaining.remove(match)

    if remaining:
        missing = "\n".join(f"  - {f.name}" for f in remaining)
        raise KeyError(
            "Certains fichiers trouvés par le pattern n'ont pas de label associé :\n"
            f"{missing}\n\n"
            "Ajoute ces fichiers dans le dictionnaire `labels`, ou rends ton pattern plus restrictif."
        )

    return ordered


def _label_from_file(fpath: Path, labels):
    """
    Retourne le label associé à un fichier.

    labels peut être :
        - None : label = stem du fichier ;
        - list/tuple : association positionnelle ;
        - dict : association directe au fichier.

    Pour un dict, les clés possibles sont :
        str(fpath), fpath.name, fpath.stem.
    """
    if labels is None:
        return fpath.stem

    if not isinstance(labels, dict):
        raise TypeError(
            "labels doit être None, une liste, ou un dictionnaire. "
            "Pour associer directement un label à chaque fichier, utilise un dictionnaire."
        )

    candidates = [
        str(fpath),
        fpath.as_posix(),
        fpath.name,
        fpath.stem,
    ]

    for key in candidates:
        if key in labels:
            return labels[key]

    raise KeyError(
        "Aucun label trouvé pour le fichier :\n"
        f"  {fpath}\n\n"
        "Clés essayées :\n"
        + "\n".join(f"  - {key}" for key in candidates)
        + "\n\n"
        "Ajoute une de ces clés dans ton dictionnaire `labels`."
    )


def _resolve_labels(csv_files, labels):
    """Résout les labels, soit positionnels, soit par dictionnaire fichier → label."""
    if labels is None:
        return [Path(f).stem for f in csv_files]

    if isinstance(labels, dict):
        return [_label_from_file(Path(f), labels) for f in csv_files]

    if len(labels) != len(csv_files):
        raise ValueError(
            "Le nombre de labels doit être égal au nombre de fichiers CSV.\n"
            f"Nombre de fichiers : {len(csv_files)}\n"
            f"Nombre de labels   : {len(labels)}\n"
            "Alternative recommandée : utilise un dictionnaire {nom_fichier: label}."
        )

    return list(labels)


def _print_mapping(csv_files, labels):
    print("\nAssociation fichiers → légendes :")
    for fpath, label in zip(csv_files, labels):
        print(f"  {Path(fpath).name:50s} -> {label}")
    print()


# ============================================================
#  FONCTION PRINCIPALE
# ============================================================

def plot_field(
    csv_files,
    formula=None,
    x_col="timestep",
    x_formula=None,
    labels=None,
    skip_first=False,
    xlabel="",
    ylabel="",
    title="",
    hlines=None,
    xlim=None,
    ylim=None,
    xscale="linear",
    yscale="linear",
    legend_title=None,
    legend_loc="best",
    colors=None,
    figsize=None,
    save_path=None,
    show=True,
    print_mapping=False,
):
    """
    Trace une ou plusieurs courbes à partir de CSV mode "field".

    labels peut être :
        - None : noms de fichiers ;
        - liste : labels dans l'ordre de csv_files ;
        - dict : {nom_fichier: label}, recommandé.
    """
    if isinstance(csv_files, (str, Path)):
        csv_files = [Path(csv_files)]
    else:
        csv_files = [Path(f) for f in csv_files]

    if not csv_files:
        raise FileNotFoundError("Aucun fichier CSV fourni.")

    resolved_labels = _resolve_labels(csv_files, labels)

    if print_mapping:
        _print_mapping(csv_files, resolved_labels)

    _colors = colors if colors is not None else COLORS
    _resolved_formula = None

    fig, ax = plt.subplots(figsize=figsize or plt.rcParams["figure.figsize"])

    for k, (fpath, label) in enumerate(zip(csv_files, resolved_labels)):
        df = pd.read_csv(fpath)
        if skip_first:
            df = df.iloc[1:].reset_index(drop=True)

        if formula is None and _resolved_formula is None:
            _resolved_formula = _auto_formula(df)
        current_formula = _resolved_formula if formula is None else formula

        x = np.asarray(x_formula(df) if x_formula is not None else df[x_col], dtype=float)
        y = np.asarray(current_formula(df), dtype=float)

        ax.plot(x, y, color=_colors[k % len(_colors)], label=label)

    if hlines:
        for yval, hlbl in hlines:
            add_hline(ax, yval, label=hlbl)

    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    style_ax(ax, xlabel=xlabel, ylabel=ylabel, title=title)
    ax.legend(title=legend_title, loc=legend_loc)

    finalize(fig, save_path=save_path, show=show)
    return fig, ax


# ============================================================
#  RACCOURCI ÉTUDE GÉNÉRIQUE
# ============================================================

def plot_field_study(
    pattern="*-field.csv",
    formula=None,
    x_col="timestep",
    x_formula=None,
    labels=None,
    skip_first=False,
    xlabel=r"Temps $t$",
    ylabel=r"$E$",
    title="",
    hlines=None,
    xlim=None,
    ylim=None,
    xscale="linear",
    yscale="linear",
    legend_title=None,
    legend_loc="best",
    colors=None,
    figsize=None,
    save_path=None,
    show=True,
    print_mapping=False,
):
    """
    Trace tous les CSV field trouvés par `pattern`.

    Important : la légende peut être donnée directement par fichier :
        labels = {
            "fichierA-field.csv": "label A",
            "fichierB-field.csv": "label B",
        }
    """
    csv_files = _order_files_from_label_dict(glob.glob(pattern), labels)

    if not csv_files:
        raise FileNotFoundError(f"Aucun fichier trouvé : {pattern}")

    return plot_field(
        csv_files=csv_files,
        formula=formula,
        x_col=x_col,
        x_formula=x_formula,
        labels=labels,
        skip_first=skip_first,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        hlines=hlines,
        xlim=xlim,
        ylim=ylim,
        xscale=xscale,
        yscale=yscale,
        legend_title=legend_title,
        legend_loc=legend_loc,
        colors=colors,
        figsize=figsize,
        save_path=save_path,
        show=show,
        print_mapping=print_mapping,
    )
