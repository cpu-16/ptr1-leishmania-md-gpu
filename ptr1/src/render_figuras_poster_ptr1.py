# -*- coding: utf-8 -*-
"""Re-render en ALTA RESOLUCIÓN de las 3 figuras estructurales de PTR1 para el póster JIC.

Genera, con PyMOL (ray-traced, fondo blanco, antialias):
  1) ptr1_tetramero.png      -> homotetrámero holo (4 subunidades + NADPH + sustrato)
  2) ptr1_sitio_activo.png   -> sitio activo (cadena A): Tyr114 específica de especie
                                sobre el sustrato + tríada catalítica Ser112/Tyr194/Lys198
  3) ptr1_superpos_1e92.png  -> superposición del modelo MD vs cristal 1E92 (RMSD 1.05 Å)

Paleta del póster: proteína en TEAL (#2A9D8F); cofactor NADPH en ÁMBAR (#E9A23B);
sustrato HBI en magenta para contraste; cristal de referencia en gris.

⚠️ NUMERACIÓN (aclarado el 10 jul 2026). El PDB de vista `md_final_holo_chains.pdb` estaba
renumerado en −1, así que ahí la Tyr de especie es `resi 113` y la tríada `resi 111/193/197`.
En `system.pdb` (la topología que se simuló) y en UniProt de *L. panamensis* esos mismos
residuos son **Tyr114** y la tríada **Ser112–Tyr194–Lys198**, y la arginina es **Arg18**.
Los residuos físicos siempre fueron los correctos; lo que estaba corrido era el rótulo.

Uso:  pymol -cq src/render_figuras_poster_ptr1.py
"""
import os
from pymol import cmd

ROOT = "."
B = f"{ROOT}/EXPLORACION/DIRECCION_B_leishmania"
OUT = f"{ROOT}/results/figures/poster"
os.makedirs(OUT, exist_ok=True)

def base_settings():
    # Paleta del póster (redefinida tras cada reinitialize, que borra colores custom)
    cmd.set_color("teal_p",  [0.165, 0.616, 0.561])   # #2A9D8F
    cmd.set_color("amber_p", [0.914, 0.635, 0.231])   # #E9A23B
    cmd.bg_color("white")
    cmd.set("ray_shadows", 0)
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("ray_trace_mode", 0)


# ===================== 1) TETRÁMERO HOLO =====================
cmd.reinitialize()
base_settings()
cmd.load(f"{B}/results/system/prod/view/md_final_holo_chains.pdb", "ptr1")
cmd.hide("everything")
cmd.show("cartoon", "polymer")
# 4 subunidades en familia fría que armoniza con navy/teal
cmd.color("teal_p",    "chain A and polymer")
cmd.color("palecyan",  "chain B and polymer")
cmd.color("lightblue", "chain C and polymer")
cmd.color("gray70",    "chain D and polymer")
# Ligandos por resname (sin restringir cadena, así siempre aparecen)
cmd.show("sticks", "resn NPH or resn HBI")
cmd.color("amber_p", "resn NPH")     # NADPH (cofactor)
cmd.color("hotpink", "resn HBI")     # sustrato (7,8-dihidrobiopterina)
cmd.set("stick_radius", 0.30, "resn NPH or resn HBI")
cmd.util.cnc("resn NPH or resn HBI")
cmd.orient("polymer")
cmd.zoom("polymer", 3)
cmd.ray(2600, 2200)
cmd.png(f"{OUT}/ptr1_tetramero.png", dpi=300)
print("[OK] 1/3 ptr1_tetramero.png  | ligandos:", cmd.count_atoms("resn NPH or resn HBI"))


# ===================== 2) SITIO ACTIVO (Tyr114) =====================
cmd.reinitialize()
base_settings()
cmd.load(f"{B}/results/system/prod/view/md_final_holo_chains.pdb", "pan")
cmd.hide("everything")
cmd.show("cartoon", "pan and chain A and polymer")
cmd.set("cartoon_transparency", 0.80, "pan and chain A")
cmd.color("teal_p", "pan and chain A")
# Sustrato + cofactor en el sitio A
cmd.show("sticks", "pan and chain A and (resn HBI or resn NPH)")
cmd.color("hotpink", "pan and chain A and resn HBI")
cmd.color("amber_p", "pan and chain A and resn NPH")
# Tríada catalítica Ser112/Tyr194/Lys198 en UniProt. Los `resi 111+193+197` de abajo son
# los del PDB de vista, que estaba renumerado en -1. Ver el aviso del docstring.
# Se selecciona por número Y por nombre de residuo: si el PDB trae otra numeración
# (p. ej. `system.pdb`, donde 111/193/197 son ALA/MET/ALA), la selección sale vacía y el
# script aborta, en vez de colorear el residuo equivocado sin decir nada.
triada = ("pan and chain A and ((resi 111 and resn SER) or (resi 193 and resn TYR) "
          "or (resi 197 and resn LYS))")
if cmd.count_atoms(triada) == 0:
    raise SystemExit(
        "\n🔴 La tríada no se encontró con `resi 111+193+197`.\n"
        "   Este PDB usa otra numeración. En `system.pdb` (UniProt) es 112/194/198.\n"
        "   Ver el aviso de NUMERACIÓN en el docstring.\n"
    )
cmd.show("sticks", triada)
cmd.color("slate", triada)
# Tyr114 específica de especie (resi 113): destacada en azul intenso (protagonista)
tyr114 = "pan and chain A and resi 113"
cmd.show("sticks", tyr114)
cmd.color("marine", tyr114)
cmd.util.cnc(f"({triada}) or ({tyr114}) or (pan and chain A and (resn HBI or resn NPH))")
# H-bond específico de especie: Tyr114-OH ... O10 del sustrato (en rojo, sin etiqueta)
cmd.distance("hb", "pan and chain A and resi 113 and name OH",
             "pan and chain A and resn HBI and name O10")
cmd.color("red", "hb")
cmd.hide("labels", "hb")
cmd.set("dash_width", 5)
# Sin etiquetas quemadas: las anotaciones se ponen como texto en el póster (más legibles)
cmd.orient("pan and chain A and resn HBI")
cmd.zoom("pan and chain A and (resn HBI or resi 111+113+193+197)", 3.5)
cmd.turn("y", 15)
cmd.ray(2600, 2000)
cmd.png(f"{OUT}/ptr1_sitio_activo.png", dpi=300)
print("[OK] 2/3 ptr1_sitio_activo.png | HBI en A:",
      cmd.count_atoms("pan and chain A and resn HBI"))


# ===================== 3) SUPERPOSICIÓN vs cristal 1E92 =====================
cmd.reinitialize()
base_settings()
cmd.load(f"{B}/results/system/prod/view/md_final_holo.pdb", "md")
cmd.load(f"{B}/data/ptr1_Lmajor_1E92.pdb", "xtal")
rms = cmd.align("md and polymer", "xtal and polymer")
cmd.hide("everything")
cmd.show("cartoon", "polymer")
cmd.color("teal_p", "md and polymer")     # modelo MD (panamensis)
cmd.color("gray70", "xtal and polymer")   # cristal 1E92 (major)
cmd.orient("md")
cmd.ray(2600, 2000)
cmd.png(f"{OUT}/ptr1_superpos_1e92.png", dpi=300)
print(f"[OK] 3/3 ptr1_superpos_1e92.png | RMSD align = {rms[0]:.2f} A sobre {rms[1]} atomos")

print("\nFiguras del póster en:", OUT)
