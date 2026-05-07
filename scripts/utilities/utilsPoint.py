"""
utilsPoint.py — Plot à partir d'un CSV mode "point"
====================================================

Version avec labels associés directement aux fichiers.

Tu peux fournir :
    labels = {
        "normVel-cfl-1e-2.csv": r"$CFL=10^{-2}$",
        "normVel-cfl-1e-1.csv": r"$CFL=10^{-1}$",
    }

Les clés acceptées sont :
    - le chemin complet du fichier ;
    - le nom du fichier avec extension ;
    - le stem du fichier sans extension.

Exemple :
    labels = {
        "../results/.../normVel-cfl-1e-2.csv": r"$CFL=10^{-2}$",
        "normVel-cfl-1e-1.csv": r"$CFL=10^{-1}$",
        "normVel-cfl-1": r"$CFL=1$",
    }
"""

import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from utilities.utilsStyle import COLORS, style_ax, add_hline, finalize


# ============================================================
#  HELPERS INTERNES
# ============================================================

def _default_formula(df: pd.DataFrame) -> np.ndarray:
    """Dernière colonne du CSV — convient pour tous les CSV mode 'point'."""
    return df.iloc[:, -1].to_numpy(dtype=float)


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
        print(f"  {Path(fpath).name:45s} -> {label}")
    print()


# ============================================================
#  FONCTION PRINCIPALE
# ============================================================

def plot_point(
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
    Trace une ou plusieurs courbes à partir de CSV mode "point".

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

    _formula = formula if formula is not None else _default_formula
    _colors = colors if colors is not None else COLORS

    fig, ax = plt.subplots(figsize=figsize or plt.rcParams["figure.figsize"])

    for k, (fpath, label) in enumerate(zip(csv_files, resolved_labels)):
        df = pd.read_csv(fpath)
        if skip_first:
            df = df.iloc[1:].reset_index(drop=True)

        x = np.asarray(x_formula(df) if x_formula is not None else df[x_col], dtype=float)
        y = np.asarray(_formula(df), dtype=float)

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

def plot_point_study(
    pattern="*.csv",
    formula=None,
    x_col="timestep",
    x_formula=None,
    labels=None,
    skip_first=False,
    xlabel=r"Temps $t$",
    ylabel=r"$u$",
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
    Trace tous les CSV trouvés par `pattern`.

    Si `labels` est un dictionnaire, l’ordre des courbes et de la légende
    suit exactement l’ordre des clés du dictionnaire.

    Important : la légende peut être donnée directement par fichier :
        labels = {
            "fichierA.csv": "label A",
            "fichierB.csv": "label B",
        }
    """
    csv_files = _order_files_from_label_dict(glob.glob(pattern), labels)

    if not csv_files:
        raise FileNotFoundError(f"Aucun fichier trouvé : {pattern}")

    return plot_point(
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
