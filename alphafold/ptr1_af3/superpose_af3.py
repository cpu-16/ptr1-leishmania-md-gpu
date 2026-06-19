# -*- coding: utf-8 -*-
"""Valida el modelo AF3 (de novo) contra nuestro modelo MD y el cristal 1E92:
RMSD del esqueleto y posición del NADPH (de novo vs trasplantado)."""
from pymol import cmd
import numpy as np, collections

B = "."
af3 = B + "/alphafold_server/extracted/fold_ptr1_lpanamensis_tetramero_nadph_model_0.cif"
mdl = B + "/EXPLORACION/DIRECCION_B_leishmania/results/system/prod/view/md_final_holo.pdb"
xtal = B + "/EXPLORACION/DIRECCION_B_leishmania/data/ptr1_Lmajor_1E92.pdb"

cmd.load(af3, "af3"); cmd.load(mdl, "mdl"); cmd.load(xtal, "xtal")

def coms(obj, resn):
    g = collections.defaultdict(list)
    for a in cmd.get_model(f"{obj} and resn {resn}").atom:
        g[(a.chain, a.resi)].append(a.coord)
    return [np.array(v).mean(0) for v in g.values()]

print("=== 1) AF3 (de novo) vs NUESTRO MODELO MD (panamensis) ===")
r1 = cmd.align("af3 and polymer", "mdl and polymer")
print(f"   RMSD esqueleto = {r1[0]:.2f} A  ({r1[1]} atomos alineados)")

# posicion del NADPH: de novo (AF3, NDP) vs trasplantado (nuestro, NPH), ya alineadas las proteinas
ndp = coms("af3", "NDP"); nph = coms("mdl", "NPH")
print(f"   NADPH de novo (AF3): {len(ndp)} | NADPH trasplantado (modelo): {len(nph)}")
dmins = []
for c in nph:
    d = min(np.linalg.norm(c - x) for x in ndp)
    dmins.append(d)
print(f"   Distancia centroide NADPH AF3<->modelo por sitio (A): {[round(d,1) for d in sorted(dmins)]}")
print(f"   -> media {np.mean(dmins):.1f} A  (si <~3 A: AF3 colocó el cofactor en el MISMO sitio)")

print("\n=== 2) AF3 (panamensis) vs CRISTAL 1E92 (major) ===")
r2 = cmd.align("af3 and polymer", "xtal and polymer")
print(f"   RMSD esqueleto = {r2[0]:.2f} A  ({r2[1]} atomos alineados)")
# NADPH de AF3 (NDP) vs NADP+ del cristal (NAP)
nap = coms("xtal", "NAP")
dmins2 = [min(np.linalg.norm(c - x) for x in nap) for c in coms("af3", "NDP")]
print(f"   Distancia NADPH(AF3)<->NADP+(cristal) por sitio (A): {[round(d,1) for d in sorted(dmins2)]}")
print(f"   -> media {np.mean(dmins2):.1f} A")

# figura de superposicion AF3 vs cristal (ya alineadas)
cmd.hide("everything"); cmd.bg_color("white"); cmd.set("ray_shadows", 0)
cmd.set("ray_opaque_background", 1)
cmd.show("cartoon", "(af3 or xtal) and polymer")
cmd.color("marine", "af3 and polymer"); cmd.color("grey70", "xtal and polymer")
cmd.show("sticks", "af3 and resn NDP"); cmd.color("yellow", "af3 and resn NDP")
cmd.show("sticks", "xtal and resn NAP"); cmd.color("orange", "xtal and resn NAP")
cmd.orient("af3"); cmd.ray(1200, 900)
cmd.png(B + "/alphafold_server/af3_vs_cristal.png", dpi=140)
print("\n[OK] figura -> alphafold_server/af3_vs_cristal.png")
