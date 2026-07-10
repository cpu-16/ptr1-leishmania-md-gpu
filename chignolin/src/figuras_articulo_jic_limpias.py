# -*- coding: utf-8 -*-
"""Regenera las figuras del ARTÍCULO JIC SIN títulos embebidos.

Atiende las sugerencias del asesor sobre formato de figuras (la convención
RIC/JIC es que el título va en el PIE de figura, no dentro de la imagen, y que
los multipaneles se rotulan con letras A/B/C/D y se describen en el pie):

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
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

PROY = "."
AN = os.path.join(PROY, "EXPLORACION", "DIRECCION_B_leishmania",
                  "results", "system", "prod", "analysis")
OUT = os.path.join(PROY, "results", "figures")
DPI = 300
# Tipografia mas grande: estas figuras se colocan grandes en el poster A0 y a
# 1-2 m el texto interno debe leerse (antes ~8pt -> ilegible).
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
    "font.size": 15, "axes.labelsize": 15, "axes.titlesize": 15,
    "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 13,
})

# Paleta de MARCA del poster, daltonico-segura (luminancias bien separadas):
# navy / teal / ambar / gris. Se usa con redundancia no-cromatica (estilos de
# linea y hachurado) para que la distincion no dependa solo del color.
PAL = {"navy": "#1F2D4E", "teal": "#2A9D8F", "amber": "#E9A23B", "grey": "#6B6B6B"}
SITE_COLOR = {"A": "#1F2D4E", "B": "#2A9D8F", "C": "#E9A23B", "D": "#6B6B6B"}
SITE_LS = {"A": "-", "B": "--", "C": "-.", "D": ":"}
SITE_HATCH = {"A": "", "B": "//", "C": "..", "D": "xx"}
# Rampa secuencial navy -> teal -> ambar (de marca, monotona en luminancia).
CMAP_BRAND = LinearSegmentedColormap.from_list(
    "brand", ["#1F2D4E", "#2A6E78", "#3FA08C", "#E9A23B"])
# Rampa para HEATMAPS: navy -> teal -> teal claro, SIN ambar, para que una figura
# de valores altos no inunde de ambar y le robe el rol de acento escaso al ambar.
CMAP_HEAT = LinearSegmentedColormap.from_list(
    "heat", ["#1F2D4E", "#2A9D8F", "#BFE3DC"])


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
    # Figura mas ancha para que la leyenda + linea de RMSD QUEPAN sin recortarse
    # (con la tipografia grande, el texto largo se salia del lienzo de 6").
    fig = plt.figure(figsize=(6.9, 5.2))
    ax = fig.add_axes([0.0, 0.16, 1.0, 0.82])
    ax.imshow(img)
    ax.axis("off")
    handles = [
        Patch(facecolor="#999999", label="Experimental (5AWL)"),
        Patch(facecolor="#1F2D4E", label="AlphaFold2 (pLDDT 94)"),
        Patch(facecolor="#E9A23B", label="MD (representativa)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=11, bbox_to_anchor=(0.5, 0.075), columnspacing=0.9,
               handletextpad=0.35)
    # Comparación homogénea: medoide del ensemble de MD frente al modelo único de
    # AlphaFold. Reportar aquí el mejor frame (0.15 Å) los haría incomparables.
    fig.text(0.5, 0.025,
             "RMSD Cα:  exp–AF2 0.90 Å  ·  exp–MD 0.85 Å  ·  AF2–MD 0.41 Å",
             ha="center", fontsize=9.5)
    out = os.path.join(OUT, "superposicion_limpia.png")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print("[OK] Fig 1 ->", out)


# --------------------------------------------------------------------------- #
# FIGURA 2 — heatmap ocupancia de puentes de H de la red catalítica (sin título)
#   Datos finales en red_catalitica_resumen.txt (100 ns, tetrámero holo).
# --------------------------------------------------------------------------- #
def figura_heatmap():
    # Numeración UniProt de L. panamensis (A0A088SA10), la misma de `system.pdb`.
    # ⚠️ Hasta el 10 jul 2026 estas filas decían Ser111/Tyr193/Lys197/Arg17: eran los
    # índices 0-based de mdtraj, uno menos que el número de residuo. Los residuos
    # físicos siempre fueron los correctos; sólo el rótulo estaba corrido.
    filas = [
        "Ser112·N2 (sustrato)",
        "Tyr194·N8 (catalítica)",
        "Lys198·ribosa-NADP",
        "Arg18·2'-P NADPH",
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

    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    # Rampa navy->teal (CMAP_HEAT), SIN ambar: como los valores son altos (83-100%)
    # casi todas las celdas caerian en el extremo calido y el ambar inundaria la
    # figura. vmin=80 reparte el rango 80-100 sobre navy->teal claro.
    im = ax.imshow(M, cmap=CMAP_HEAT, vmin=80, vmax=100, aspect="auto")
    ax.set_xticks(range(len(cols)), labels=cols)
    ax.set_yticks(range(len(filas)), labels=filas)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            # Celdas oscuras (navy/teal, valores ~80-92) -> texto blanco;
            # celdas claras (teal palido, ~93-100) -> texto negro.
            ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                    color="white" if v < 93 else "black", fontsize=15)
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
    ax.text(0.02, 0.97, letra, transform=ax.transAxes, fontsize=17,
            fontweight="bold", va="top", ha="left")


def figura_convergencia():
    r = np.genfromtxt(os.path.join(AN, "rmsd_all.dat"), skip_header=1)
    rmsd_t = (r[:, 0] - 1) * DT_NS
    rmsd = r[:, 1]
    d = np.genfromtxt(os.path.join(AN, "dist_cat_mdtraj.dat"), skip_header=1)
    dt_ns = d[:, 0]
    dcat = {CHAINS[k]: d[:, k + 1] for k in range(4)}

    fig, ax = plt.subplots(2, 2, figsize=(9.6, 6.6))

    # A — RMSD Cα: crudo + media acumulada (teal claro + navy, paleta de marca)
    ax[0, 0].plot(rmsd_t, rmsd, color="#A9D6CE", lw=0.9, label="RMSD por frame")
    ax[0, 0].plot(rmsd_t, _running_mean(rmsd), color=PAL["navy"], lw=2.6, label="media acumulada")
    ax[0, 0].set_xlabel("tiempo (ns)"); ax[0, 0].set_ylabel("RMSD Cα (Å)")
    ax[0, 0].legend(); ax[0, 0].grid(alpha=0.3)
    _panel_label(ax[0, 0], "A")

    # B — RMSD por bloques de 20 ns (teal)
    bl = _blocks(rmsd_t, rmsd)
    bx = [b[0] for b in bl]; bm = [b[1] for b in bl]; bs = [b[2] for b in bl]
    ax[0, 1].bar(bx, bm, width=BLOCK_NS * 0.8, yerr=bs, capsize=4, color=PAL["teal"])
    ax[0, 1].set_xlabel("tiempo (ns)"); ax[0, 1].set_ylabel("RMSD Cα medio (Å)")
    ax[0, 1].grid(axis="y", alpha=0.3)
    _panel_label(ax[0, 1], "B")

    # C — distancia catalítica por sitio: paleta de marca + ESTILO DE LINEA
    # distinto por sitio (redundancia no-cromatica, daltonico-segura).
    for s in CHAINS:
        ax[1, 0].plot(dt_ns, _running_mean(dcat[s]), color=SITE_COLOR[s],
                      ls=SITE_LS[s], lw=2.4, label=f"sitio {s}")
    ax[1, 0].set_xlabel("tiempo (ns)"); ax[1, 0].set_ylabel("distancia C4N–C6 (Å)")
    ax[1, 0].set_ylim(3.2, 5.2); ax[1, 0].legend(ncol=2); ax[1, 0].grid(alpha=0.3)
    _panel_label(ax[1, 0], "C")

    # D — distancia por bloques, agrupada por sitio: paleta de marca + HACHURADO
    # distinto por sitio (redundancia no-cromatica).
    nb = len(_blocks(dt_ns, dcat["A"]))
    x = np.arange(nb); w = 0.2
    for j, s in enumerate(CHAINS):
        bl = _blocks(dt_ns, dcat[s])
        ax[1, 1].bar(x + j * w, [b[1] for b in bl], w, yerr=[b[2] for b in bl],
                     capsize=2, color=SITE_COLOR[s], hatch=SITE_HATCH[s],
                     edgecolor="white", linewidth=0.4, label=f"sitio {s}")
    ax[1, 1].set_xticks(x + 1.5 * w)
    ax[1, 1].set_xticklabels(
        [f"{int(b[0]-BLOCK_NS/2)}-{int(b[0]+BLOCK_NS/2)}" for b in _blocks(dt_ns, dcat['A'])])
    ax[1, 1].set_xlabel("bloque (ns)"); ax[1, 1].set_ylabel("distancia media (Å)")
    ax[1, 1].set_ylim(3.2, 5.2); ax[1, 1].legend(ncol=2); ax[1, 1].grid(axis="y", alpha=0.3)
    _panel_label(ax[1, 1], "D")

    fig.tight_layout()
    out = os.path.join(OUT, "convergencia_limpia.png")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print("[OK] Fig 3 ->", out)


def figura_convergencia_poster():
    """Version POSTER (1x2, fuentes grandes, legible a 1-2 m): RMSD Cα en el tiempo
    + distancia catalitica por sitio. El 2x2 detallado (figura_convergencia) queda
    para el articulo; aqui se evita el texto microscopico de 4 paneles en un slot
    pequeno (conflicto A2 'heroe >=4x apoyo' vs C2 legibilidad)."""
    r = np.genfromtxt(os.path.join(AN, "rmsd_all.dat"), skip_header=1)
    rmsd_t = (r[:, 0] - 1) * DT_NS
    rmsd = r[:, 1]
    d = np.genfromtxt(os.path.join(AN, "dist_cat_mdtraj.dat"), skip_header=1)
    dt_ns = d[:, 0]
    dcat = {CHAINS[k]: d[:, k + 1] for k in range(4)}

    fig, ax = plt.subplots(1, 2, figsize=(7.8, 3.6))
    ax[0].plot(rmsd_t, rmsd, color="#A9D6CE", lw=1.0, label="RMSD por frame")
    ax[0].plot(rmsd_t, _running_mean(rmsd), color=PAL["navy"], lw=2.8, label="media acumulada")
    ax[0].set_xlabel("tiempo (ns)"); ax[0].set_ylabel("RMSD Cα (Å)")
    ax[0].grid(alpha=0.3)
    # leyendas ARRIBA de cada panel (fuera de los ejes) para no tapar las curvas
    ax[0].legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False,
                 fontsize=11, columnspacing=1.0, handlelength=1.6, handletextpad=0.5)
    for s in CHAINS:
        ax[1].plot(dt_ns, _running_mean(dcat[s]), color=SITE_COLOR[s], ls=SITE_LS[s],
                   lw=2.4, label=f"sitio {s}")
    ax[1].set_xlabel("tiempo (ns)"); ax[1].set_ylabel("distancia C4N–C6 (Å)")
    ax[1].set_ylim(3.2, 5.2); ax[1].grid(alpha=0.3)
    ax[1].legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4, frameon=False,
                 fontsize=10, columnspacing=0.8, handlelength=1.3, handletextpad=0.35)
    fig.subplots_adjust(top=0.82, bottom=0.17, left=0.10, right=0.985, wspace=0.34)
    out = os.path.join(OUT, "convergencia_poster.png")
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    print("[OK] Fig 3b (poster) ->", out)


def main():
    os.makedirs(OUT, exist_ok=True)
    figura_superposicion()
    figura_heatmap()
    figura_convergencia()
    figura_convergencia_poster()
    print("\nOK -> figuras limpias (sin título embebido) en", OUT)


if __name__ == "__main__":
    main()
