# -*- coding: utf-8 -*-
"""Genera el banner de cabecera del repositorio (assets/banner.png).

Identidad visual coherente con el póster/artículo del proyecto:
navy + teal + ámbar sobre blanco, tipografía DejaVu Serif/Sans.
Elemento distintivo: el mensaje central codificado por color
(teal = lo accesible/logrado; ámbar = lo que aún no alcanza).
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
TETRA = os.path.join(REPO, "figures", "ptr1", "ptr1_tetramero.png")
OUT = os.path.join(HERE, "banner.png")

NAVY = "#1F2D4E"; TEAL = "#2A9D8F"; AMBER = "#E9A23B"; INK = "#232323"; GREY = "#6B7280"
W, H = 1280, 430

fig = plt.figure(figsize=(W / 100, H / 100), dpi=150)
fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off"); ax.invert_yaxis()

# Acento lateral teal (firma sobria)
ax.add_patch(plt.Rectangle((0, 0), 10, H, color=TEAL))

LX = 54
ax.text(LX, 48, "JIC 2026   ·   BIOLOGÍA ESTRUCTURAL COMPUTACIONAL ACCESIBLE",
        color=TEAL, fontsize=11, fontweight="bold", family="DejaVu Sans", va="center")

ax.text(LX, 100, "De AlphaFold a la dinámica molecular",
        color=NAVY, fontsize=24, fontweight="bold", family="DejaVu Serif", va="center")
ax.text(LX, 140, "accesible en una GPU de consumo",
        color=NAVY, fontsize=24, fontweight="bold", family="DejaVu Serif", va="center")

ax.add_patch(plt.Rectangle((LX, 174), 380, 2.2, color=TEAL))

# Mensaje central: el color codifica el significado
ax.text(LX, 212, "Recuperar la estructura de una proteína es accesible.",
        color=TEAL, fontsize=14.5, fontweight="bold", family="DejaVu Sans", va="center")
ax.text(LX, 240, "Medir su dinámica completa, todavía no.",
        color=AMBER, fontsize=14.5, fontweight="bold", family="DejaVu Sans", va="center")

ax.text(LX, 286, "Control estructural con Chignolin   ·   Aplicación a la PTR1 de Leishmania panamensis",
        color=INK, fontsize=11.5, family="DejaVu Sans", style="italic", va="center")

ax.text(LX, 366, "OpenMM   ·   AlphaFold 2/3   ·   Markov State Model   ·   1× RTX 4060   ·   software libre",
        color=GREY, fontsize=10.5, fontweight="bold", family="DejaVu Sans", va="center")

# Figura-héroe: el tetrámero holo a la derecha (franja angosta, sin invadir el texto)
if os.path.exists(TETRA):
    axi = fig.add_axes([0.725, 0.05, 0.26, 0.90])
    axi.imshow(mpimg.imread(TETRA)); axi.axis("off")

fig.savefig(OUT, dpi=150, facecolor="white")
print("Banner ->", OUT)
