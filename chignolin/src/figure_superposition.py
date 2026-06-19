"""Figura de superposición: las 3 estructuras de Chignolin coinciden.

Carga la estructura experimental (5AWL), la de AlphaFold (rank_001) y el mejor
frame de la MD, las superpone y renderiza un PNG (headless).

Uso: conda run -n jic-folding pymol -cq src/figure_superposition.py
"""
import os
from pymol import cmd

os.makedirs("results/figures", exist_ok=True)

cmd.load("data/5awl.pdb", "exp")
cmd.load("data/chignolin_af.pdb", "af")
cmd.load("data/md_best.pdb", "md")

cmd.remove("not polymer")      # quitar aguas/iones cristalográficos
cmd.remove("hydrogens")

# superponer AF y MD sobre la experimental
cmd.align("af", "exp")
cmd.align("md", "exp")

cmd.hide("everything")
cmd.show("cartoon")
cmd.show("sticks")
cmd.set("cartoon_transparency", 0.25)
cmd.set("stick_radius", 0.15)
cmd.color("grey60", "exp")     # experimental
cmd.color("marine", "af")      # AlphaFold
cmd.color("orange", "md")      # MD
cmd.bg_color("white")
cmd.set("ray_opaque_background", 1)
cmd.orient()
cmd.ray(1400, 1000)
cmd.png("results/figures/superposition_raw.png", dpi=150)
print("OK -> results/figures/superposition_raw.png")
