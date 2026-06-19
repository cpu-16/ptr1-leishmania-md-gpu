# -*- coding: utf-8 -*-
"""
Arma el TETRÁMERO HOLO de PTR1 de *Leishmania panamensis* por superposición
del monómero de AlphaFold DB sobre cada cadena del cristal 1E92 (*L. major*).

Estrategia (paso 1 de la hoja de ruta de PTR1):
  1. Cargar la plantilla 1E92 (tetrámero A/B/C/D con NAP=NADP+ y HBI=dihidrobiopterina).
  2. Por cada cadena A/B/C/D: copiar el monómero AF y superponerlo (structure-based)
     sobre esa cadena de 1E92, luego renombrar su cadena.
  3. Trasplantar los ligandos NAP + HBI desde 1E92 (quedan en su sitio relativo al
     tetrámero, porque los monómeros AF se movieron al marco de 1E92).
  4. Guardar el tetrámero holo: 4 x proteína (L. panamensis, full-length) + 4 NAP + 4 HBI.

NOTA: el NAP queda como NADP+ en esta etapa. La conversión a NADPH es el paso
siguiente, separado, porque requiere reconstruir la topología del cofactor
(H en C4 de la nicotinamida, carga, geometría) para la parametrización.

Aguas cristalográficas y EDO se descartan; la solvatación posterior repone el agua.
(Conservar las aguas estructurales del sitio activo queda como mejora futura.)

Uso:  pymol -cq src/build_tetramer.py
"""
import os
from pymol import cmd

# Ruta base fija de la Rama B (PyMOL no resuelve __file__ de forma fiable en modo -cq).
BASE = "."
TEMPLATE = os.path.join(BASE, "data", "ptr1_Lmajor_1E92.pdb")
AF_MONO  = os.path.join(BASE, "data", "ptr1_Lpanamensis_alphafold.pdb")
OUTDIR   = os.path.join(BASE, "results", "tetramer")
OUT_PDB  = os.path.join(OUTDIR, "ptr1_Lpanamensis_tetramer_holo.pdb")
OUT_LOG  = os.path.join(OUTDIR, "build_tetramer.log")

CHAINS = ["A", "B", "C", "D"]

os.makedirs(OUTDIR, exist_ok=True)
log_lines = []

def log(msg):
    print(msg)
    log_lines.append(str(msg))

# ---------------------------------------------------------------------------
cmd.reinitialize()
cmd.load(TEMPLATE, "tmpl")
cmd.load(AF_MONO, "afmono")

log("# Armado del tetrámero PTR1 L. panamensis (superposición sobre 1E92)")
log(f"Plantilla : {TEMPLATE}")
log(f"Monómero  : {AF_MONO}")
log(f"Átomos plantilla : {cmd.count_atoms('tmpl')}")
log(f"Átomos monómero  : {cmd.count_atoms('afmono')}")
log("")
log("RMSD de superposición por cadena (super, sobre proteína):")

# 1+2) copiar y superponer el monómero sobre cada cadena de la plantilla
prot_objs = []
for ch in CHAINS:
    obj = f"prot_{ch}"
    cmd.create(obj, "afmono")
    # super: mueve 'obj' (mobile) sobre la cadena ch de la plantilla (target).
    # Restringimos el target a proteína de esa cadena; super alinea por estructura.
    res = cmd.super(f"{obj}", f"tmpl and chain {ch} and polymer.protein")
    rmsd, natoms = res[0], res[1]
    # renombrar la cadena de la copia a la cadena destino
    cmd.alter(obj, f"chain='{ch}'")
    cmd.alter(obj, "segi=''")
    cmd.sort(obj)
    prot_objs.append(obj)
    log(f"  cadena {ch}: RMSD = {rmsd:6.3f} A sobre {natoms} atomos alineados")

# 3) ligandos desde la plantilla (en su marco, que ahora es el del tetrámero AF)
cmd.create("ligs", "tmpl and (resn NAP+HBI)")
n_nap = cmd.count_atoms("ligs and resn NAP")
n_hbi = cmd.count_atoms("ligs and resn HBI")
log("")
log(f"Ligandos trasplantados: NAP (NADP+) = {n_nap} atomos, HBI (dihidrobiopterina) = {n_hbi} atomos")

# 4) combinar y guardar
sel_final = " or ".join(prot_objs + ["ligs"])
cmd.create("tetramer", sel_final)

n_prot = cmd.count_atoms("tetramer and polymer.protein")
n_total = cmd.count_atoms("tetramer")
log("")
log(f"Tetrámero final: {n_total} atomos totales ({n_prot} proteína, {n_nap+n_hbi} ligando)")

cmd.set("pdb_use_ter_records", 1)
cmd.save(OUT_PDB, "tetramer")
log(f"\nGuardado: {OUT_PDB}")

with open(OUT_LOG, "w") as fh:
    fh.write("\n".join(log_lines) + "\n")
print(f"\nLog: {OUT_LOG}")
