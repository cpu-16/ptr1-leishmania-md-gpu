# -*- coding: utf-8 -*-
"""Valida AF3 (de novo) por MONÓMERO (una subunidad) contra nuestro modelo MD y el
cristal 1E92, y mide si el NADPH cae en el mismo bolsillo. Alinear el tetrámero
completo falla por la simetría de las 4 cadenas idénticas; el monómero es la
comparación correcta del plegamiento."""
from pymol import cmd
import numpy as np, collections

B = "."
cmd.load(B + "/alphafold_server/extracted/fold_ptr1_lpanamensis_tetramero_nadph_model_0.cif", "af3")
cmd.load(B + "/EXPLORACION/DIRECCION_B_leishmania/results/system/prod/view/md_final_holo.pdb", "mdl")
cmd.load(B + "/EXPLORACION/DIRECCION_B_leishmania/data/ptr1_Lmajor_1E92.pdb", "xtal")


def com(sel):
    m = cmd.get_model(sel)
    return np.array([a.coord for a in m.atom]).mean(0) if m.atom else None


def coms(sel):
    g = collections.defaultdict(list)
    for a in cmd.get_model(sel).atom:
        g[(a.chain, a.resi)].append(a.coord)
    return [np.array(v).mean(0) for v in g.values()]


print("=== MONOMERO: AF3 cadena A vs NUESTRO MODELO MD (subunidad A, resi 1-288) ===")
r = cmd.align("af3 and chain A and polymer", "mdl and resi 1-288 and polymer")
print("   RMSD monomero = %.2f A (%d atomos)" % (r[0], r[1]))
e = com("af3 and chain E and resn NDP")            # NDP unido a la proteina cadena A
dmin = min(np.linalg.norm(e - x) for x in coms("mdl and resn NPH"))
print("   NADPH AF3(de novo) vs nuestro NADPH(trasplantado): centroide a %.1f A" % dmin)
print("   -> si <~3-4 A, AF3 coloco el cofactor en el MISMO bolsillo")

print("=== MONOMERO: AF3 cadena A vs CRISTAL 1E92 cadena A (major) ===")
r2 = cmd.align("af3 and chain A and polymer", "xtal and chain A and polymer")
print("   RMSD monomero = %.2f A (%d atomos)" % (r2[0], r2[1]))
e2 = com("af3 and chain E and resn NDP")
nap = coms("xtal and chain A and resn NAP")
if nap:
    print("   NADPH AF3 vs NADP+ cristal: centroide a %.1f A" % min(np.linalg.norm(e2 - x) for x in nap))

# figura monomero af3 vs cristal (ya alineadas)
cmd.hide("everything"); cmd.bg_color("white")
cmd.set("ray_shadows", 0); cmd.set("ray_opaque_background", 1)
cmd.show("cartoon", "((af3 and chain A) or (xtal and chain A)) and polymer")
cmd.color("marine", "af3 and chain A"); cmd.color("grey70", "xtal and chain A")
cmd.show("sticks", "af3 and chain E and resn NDP"); cmd.color("yellow", "af3 and chain E")
cmd.show("sticks", "xtal and chain A and resn NAP"); cmd.color("orange", "xtal and chain A and resn NAP")
cmd.orient("af3 and chain A"); cmd.ray(1200, 900)
cmd.png(B + "/alphafold_server/af3_vs_cristal_monomero.png", dpi=140)
print("[OK] figura -> alphafold_server/af3_vs_cristal_monomero.png")
