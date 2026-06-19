# -*- coding: utf-8 -*-
"""Regenera las figuras del ARTÍCULO JIC SIN títulos embebidos.

Sigue la convención RIC/JIC de formato de figuras (el título va en el PIE de
figura, no dentro de la imagen, y los multipaneles se rotulan con letras A/B/C/D
y se describen en el pie):

  Fig 1  superposicion_limpia.png       <- Chignolin (sin título; solo clave de color + RMSD)
  Fig 2  red_catalitica_heatmap_limpia.png <- heatmap ocupancia H-bond (sin título)
  Fig 3  convergencia_limpia.png         <- 2x2 con paneles A/B/C/D (sin títulos de panel)

NO modifica los scripts del póster (que sí llevan título): escribe copias nuevas.

Uso:  conda run -n jic-folding python src/figuras_articulo_jic_limpias.py
Salida: results/figures/*_limpia.png  (300 dpi)
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Patch
import numpy as np

PROY = "."
AN = os.path.join(PROY, "EXPLORACION", "DIRECCION_B_leishmania",
                  "results", "system", "prod", "analysis")
OUT = os.path.join(PROY, "results", "figures")
DPI = 300
plt.rcParams.update({"font.family": "serif", "font.size": 11})


# --------------------------------------------------------------------------- #
# FIGURA 1 — superposición Chignolin (sin título; clave de color + RMSD al pie)
# --------------------------------------------------------------------------- #
def _recortar_blanco(img, umbral=0.96, margen=0.015):
    """Recorta el margen blanco alrededor de la molécula para que llene el
    recuadro (sugerencia de revisión: usar el espacio en blanco de la figura).
    """
    rgb = img[..., :3]
    tinta = (rgb < umbral).any(axis=2)
    if img.shape[-1] == 4:                 # respeta el canal alfa si existe
        tinta &= img[..., 3] > 0.05
    ys, xs = np.where(tinta)
    if len(xs) == 0:
        return img
    h, w = img.shape[0], img.shape[1]
    py, px = int((ys.max() - ys.min()) * margen), int((xs.max() - xs.min()) * margen)
    y0, y1 = max(0, ys.min() - py), min(h, ys.max() + py)
    x0, x1 = max(0, xs.min() - px), min(w, xs.max() + px)
    return img[y0:y1, x0:x1]


def figura_superposicion():
    raw = os.path.join(OUT, "superposition_raw.png")
    img = _recortar_blanco(mpimg.imread(raw))   # molécula sin margen blanco
    fig = plt.figure(figsize=(6.0, 5.0))
    # Aprovechamos el espacio en blanco inferior de la figura para dar más
    # tamaño a las etiquetas (leyenda de color y RMSD), según sugerencia de
    # revisión: reservamos un 17 % inferior y subimos el tamaño de letra.
    ax = fig.add_axes([0.0, 0.15, 1.0, 0.83])
    ax.imshow(img)
    ax.axis("off")
    handles = [
        Patch(facecolor="#999999", label="Experimental (PDB 5AWL)"),
        Patch(facecolor="#3a6fb0", label="AlphaFold2 (pLDDT 94)"),
        Patch(facecolor="#e8820e", label="MD (mejor frame)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, 0.065), columnspacing=1.2,
               handletextpad=0.4)
    fig.text(0.5, 0.02,
             "RMSD Cα:  experimental–AlphaFold2 0.90 Å  ·  experimental–MD 0.15 Å  ·  "
             "AlphaFold2–MD 0.97 Å",
             ha="center", fontsize=8.5)
    out = os.path.join(OUT, "superposicion_limpia.png")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print("[OK] Fig 1 ->", out)


# --------------------------------------------------------------------------- #
# FIGURA 2 — heatmap ocupancia de puentes de H de la red catalítica (sin título)
#   Datos finales en red_catalitica_resumen.txt (100 ns, tetrámero holo).
# --------------------------------------------------------------------------- #
def figura_heatmap():
    filas = [
        "Ser111·N2 (sustrato)",
        "Tyr193·N8 (catalítica)",
        "Lys197·ribosa-NADP",
        "Arg17·2'-P NADPH",
        "NADPH·N1 (cofactor-sust.)",
    ]
    cols = ["sitio A", "sitio B", "sitio C", "sitio D"]
    M = np.array([
        [97, 85, 99, 92],
        [99, 100, 100, 98],
        [100, 100, 100, 100],
        [100, 100, 83, 94],
        [100, 100, 100, 98],
    ], dtype=float)

    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    im = ax.imshow(M, cmap="YlGn", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(cols)), labels=cols)
    ax.set_yticks(range(len(filas)), labels=filas)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                    color="white" if v >= 70 else "black", fontsize=10)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("ocupancia (%)")
    ax.set_xticks(np.arange(-0.5, len(cols), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(filas), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", length=0)
    fig.tight_layout()
    out = os.path.join(OUT, "red_catalitica_heatmap_limpia.png")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print("[OK] Fig 2 ->", out)


# --------------------------------------------------------------------------- #
# FIGURA 3 — convergencia 2x2 con paneles A/B/C/D (sin títulos de panel)
#   Misma lógica que EXPLORACION/.../src/convergencia_bloques.py, sin títulos.
# --------------------------------------------------------------------------- #
DT_NS = 0.05
BLOCK_NS = 20.0
CHAINS = "ABCD"


def _running_mean(y):
    return np.cumsum(y) / np.arange(1, len(y) + 1)


def _blocks(t_ns, y, blk=BLOCK_NS):
    out = []
    edges = np.arange(0, t_ns[-1] + 1e-9, blk)
    for e in edges:
        m = (t_ns >= e) & (t_ns < e + blk)
        if m.sum() >= 2:
            out.append((e + blk / 2, float(y[m].mean()), float(y[m].std())))
    return out


def _panel_label(ax, letra):
    ax.text(0.02, 0.97, letra, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top", ha="left")


def figura_convergencia():
    r = np.genfromtxt(os.path.join(AN, "rmsd_all.dat"), skip_header=1)
    rmsd_t = (r[:, 0] - 1) * DT_NS
    rmsd = r[:, 1]
    d = np.genfromtxt(os.path.join(AN, "dist_cat_mdtraj.dat"), skip_header=1)
    dt_ns = d[:, 0]
    dcat = {CHAINS[k]: d[:, k + 1] for k in range(4)}

    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))

    # A — RMSD Cα: crudo + media acumulada
    ax[0, 0].plot(rmsd_t, rmsd, color="#9ecae1", lw=0.6, label="RMSD por frame")
    ax[0, 0].plot(rmsd_t, _running_mean(rmsd), color="#08519c", lw=2, label="media acumulada")
    ax[0, 0].set_xlabel("tiempo (ns)"); ax[0, 0].set_ylabel("RMSD Cα (Å)")
    ax[0, 0].legend(fontsize=8); ax[0, 0].grid(alpha=0.3)
    _panel_label(ax[0, 0], "A")

    # B — RMSD por bloques de 20 ns
    bl = _blocks(rmsd_t, rmsd)
    bx = [b[0] for b in bl]; bm = [b[1] for b in bl]; bs = [b[2] for b in bl]
    ax[0, 1].bar(bx, bm, width=BLOCK_NS * 0.8, yerr=bs, capsize=4, color="#6baed6")
    ax[0, 1].set_xlabel("tiempo (ns)"); ax[0, 1].set_ylabel("RMSD Cα medio (Å)")
    ax[0, 1].grid(axis="y", alpha=0.3)
    _panel_label(ax[0, 1], "B")

    # C — distancia catalítica: media acumulada por sitio
    colors = {"A": "#e41a1c", "B": "#377eb8", "C": "#4daf4a", "D": "#984ea3"}
    for s in CHAINS:
        ax[1, 0].plot(dt_ns, _running_mean(dcat[s]), color=colors[s], lw=1.8, label=f"sitio {s}")
    ax[1, 0].set_xlabel("tiempo (ns)"); ax[1, 0].set_ylabel("distancia C4N–C6 (Å)")
    ax[1, 0].set_ylim(3.2, 5.2); ax[1, 0].legend(fontsize=8, ncol=2); ax[1, 0].grid(alpha=0.3)
    _panel_label(ax[1, 0], "C")

    # D — distancia catalítica por bloques de 20 ns, agrupada por sitio
    nb = len(_blocks(dt_ns, dcat["A"]))
    x = np.arange(nb); w = 0.2
    for j, s in enumerate(CHAINS):
        bl = _blocks(dt_ns, dcat[s])
        ax[1, 1].bar(x + j * w, [b[1] for b in bl], w, yerr=[b[2] for b in bl],
                     capsize=2, color=colors[s], label=f"sitio {s}")
    ax[1, 1].set_xticks(x + 1.5 * w)
    ax[1, 1].set_xticklabels(
        [f"{int(b[0]-BLOCK_NS/2)}-{int(b[0]+BLOCK_NS/2)}" for b in _blocks(dt_ns, dcat['A'])],
        fontsize=8)
    ax[1, 1].set_xlabel("bloque (ns)"); ax[1, 1].set_ylabel("distancia media (Å)")
    ax[1, 1].set_ylim(3.2, 5.2); ax[1, 1].legend(fontsize=8, ncol=2); ax[1, 1].grid(axis="y", alpha=0.3)
    _panel_label(ax[1, 1], "D")

    fig.tight_layout()
    out = os.path.join(OUT, "convergencia_limpia.png")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print("[OK] Fig 3 ->", out)


def main():
    os.makedirs(OUT, exist_ok=True)
    figura_superposicion()
    figura_heatmap()
    figura_convergencia()
    print("\nOK -> figuras limpias (sin título embebido) en", OUT)


if __name__ == "__main__":
    main()
