"""Renderiza la trayectoria de Chignolin a imágenes PNG (sin abrir ventana).

Se ejecuta con PyMOL en modo headless:
    conda run -n jic-folding pymol -cq src/render_movie.py

Deja los frames en results/movie/. Luego se pueden unir en un GIF con ffmpeg
o ImageMagick (ver README / instrucciones).
"""
import glob
import os

from pymol import cmd

TOP = "results/prep/system.pdb"
OUTDIR = "results/movie"
N_FRAMES = 60          # nº aprox. de fotogramas del clip


def main():
    dcds = sorted(glob.glob("results/trajectories/round0_*.dcd"))
    if not dcds:
        raise SystemExit("No hay trayectorias todavía.")
    dcd = dcds[0]
    os.makedirs(OUTDIR, exist_ok=True)

    cmd.load(TOP, "chig")
    cmd.load_traj(dcd, "chig", 1)
    cmd.remove("not polymer")   # quitar agua e iones (robusto)
    cmd.remove("hydrogens")
    cmd.dss()
    cmd.hide("everything")
    cmd.show("cartoon")
    cmd.show("sticks")
    cmd.color("green", "elem C")
    cmd.color("red", "elem O")
    cmd.color("blue", "elem N")
    cmd.set("cartoon_transparency", 0.3)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.intra_fit("chig")
    cmd.orient()

    n = cmd.count_states("chig")
    step = max(1, n // N_FRAMES)
    i = 0
    for s in range(1, n + 1, step):
        cmd.set("state", s)
        cmd.ray(640, 480)
        cmd.png(os.path.join(OUTDIR, f"frame_{i:04d}.png"), dpi=100)
        i += 1
    print(f"Renderizados {i} frames (de {n} estados) en {OUTDIR}/")


main()
