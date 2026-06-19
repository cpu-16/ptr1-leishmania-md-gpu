# -*- coding: utf-8 -*-
"""
Validación biológica del tetrámero holo de PTR1 (paso 1, post-revisión de Codex).

Responde "¿es 1E92 una plantilla legítima para L. panamensis?" con números:
  - Identidad de secuencia global L. panamensis (AF) vs L. major (1E92 cadena A),
    emparejando residuos por Cα más cercano en el marco ya superpuesto.
  - Conservación del SITIO ACTIVO: residuos a < 8 Å del cofactor (NAP) y del
    sustrato (HBI), y qué sustituciones hay entre especies.

Uso (env jic-folding):  python src/validate_tetramer.py
"""
import os
import numpy as np

BASE = "."
TMPL = os.path.join(BASE, "data", "ptr1_Lmajor_1E92.pdb")
TET  = os.path.join(BASE, "results", "tetramer", "ptr1_Lpanamensis_tetramer_holo.pdb")
OUT_LOG = os.path.join(BASE, "results", "tetramer", "validate_tetramer.log")

THREE2ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

log_lines = []
def log(m):
    print(m); log_lines.append(str(m))

def parse(path, want_chain=None):
    res = {}; lig = {}
    for ln in open(path):
        rec = ln[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        ch = ln[21]
        if want_chain and ch != want_chain:
            continue
        resn = ln[17:20].strip(); resi = ln[22:27].strip()
        name = ln[12:16].strip(); el = ln[76:78].strip()
        x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
        if rec == "ATOM" and resn in THREE2ONE:
            if name == "CA":
                res[(ch, resi)] = (resn, np.array([x, y, z]))
        elif rec == "HETATM" and resn in ("NAP", "HBI") and el != "H":
            lig.setdefault((resn, ch), []).append([x, y, z])
    return res, lig

tmpl_res, tmpl_lig = parse(TMPL, want_chain="A")
af_res, _ = parse(TET, want_chain="A")

af_items = list(af_res.items())
af_xyz = np.array([v[1] for _, v in af_items])

pairs = []
for (ch, resi), (resn, ca) in tmpl_res.items():
    d = np.sqrt(((af_xyz - ca) ** 2).sum(1)); j = d.argmin()
    pairs.append((resn, af_items[j][1][0], d[j], resi))

good = [(a, b) for a, b, dist, _ in pairs if dist < 2.0]
ident = sum(1 for a, b in good if a == b)
log("# Validación biológica del tetrámero PTR1 L. panamensis (plantilla 1E92)")
log(f"Residuos 1E92 cadena A emparejados (Cα<2Å): {len(good)} de {len(tmpl_res)} modelados")
log(f"IDENTIDAD GLOBAL (sobre emparejados): {ident}/{len(good)} = {100*ident/len(good):.1f}%")

def near_lig(ligname, cutoff=8.0):
    L = np.array(tmpl_lig[(ligname, 'A')])
    out = []
    for (resn, af_resn, dist, resi) in pairs:
        ca = tmpl_res[('A', resi)][1]
        dmin = np.sqrt(((L - ca) ** 2).sum(1)).min()
        if dmin < cutoff:
            ri = int(resi) if resi.lstrip('-').isdigit() else resi
            out.append((ri, resn, af_resn, dmin, dist))
    return sorted(out, key=lambda t: t[3])

for lg in ("NAP", "HBI"):
    rows = near_lig(lg)
    cons = sum(1 for r in rows if r[1] == r[2])
    log(f"\n== Sitio de {lg} (Cα < 8 Å del ligando), cadena A: {len(rows)} residuos ==")
    log(f"   conservados L.major→L.panamensis: {cons}/{len(rows)} = {100*cons/len(rows):.0f}%")
    diffs = [r for r in rows if r[1] != r[2]]
    if diffs:
        log("   sustituciones: " + ", ".join(
            f"{THREE2ONE.get(r[1], r[1])}{r[0]}{THREE2ONE.get(r[2], r[2])}({r[3]:.1f}Å)" for r in diffs))
    else:
        log("   sitio 100% conservado")

with open(OUT_LOG, "w") as fh:
    fh.write("\n".join(log_lines) + "\n")
print(f"\nLog: {OUT_LOG}")
