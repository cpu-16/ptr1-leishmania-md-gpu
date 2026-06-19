# -*- coding: utf-8 -*-
"""
Repara los chain IDs del PDB del tetrámero exportado por mdtraj.

PROBLEMA: mdtraj colapsa las 4 cadenas del tetrámero PTR1 en una sola ("A"),
porque el prmtop de AMBER numera los residuos de forma corrida (0..1151) sin
chain. Eso rompe cualquier selección "chain A/B/C/D" en PyMOL y deja las figuras
del sitio activo confusas (un solo "chain A" = todo el tetrámero).

ARREGLO (mínimo, no cambia coordenadas ni resSeq): reescribe la columna 22
(chainID) del PDB por rango de residuo. 288 residuos de proteína por protómero:
  A: resSeq    0–287   B: 288–575   C: 576–863   D: 864–1151
Ligandos (una copia por sitio, en orden de resSeq):
  NPH 1152→A 1153→B 1154→C 1155→D     HBI 1156→A 1157→B 1158→C 1159→D
Aguas/iones -> chain 'W' (irrelevantes para la figura del sitio).

Uso:  python src/reasignar_chains.py [pdb_in] [pdb_out]
Por defecto repara view/md_final_holo.pdb -> view/md_final_holo_chains.pdb
"""
import sys
import os

BASE = "."
PROT = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
RESPP = 288   # residuos de proteína por protómero


def chain_for(resname, resseq):
    if resname in PROT:
        idx = resseq // RESPP
        return chr(ord("A") + idx) if 0 <= idx <= 3 else "P"
    if resname == "NPH":
        return chr(ord("A") + (resseq - 1152)) if 1152 <= resseq <= 1155 else "L"
    if resname == "HBI":
        return chr(ord("A") + (resseq - 1156)) if 1156 <= resseq <= 1159 else "L"
    return "W"   # agua / iones


def main():
    pin = sys.argv[1] if len(sys.argv) > 1 else os.path.join(BASE, "results/system/prod/view/md_final_holo.pdb")
    pout = sys.argv[2] if len(sys.argv) > 2 else os.path.join(BASE, "results/system/prod/view/md_final_holo_chains.pdb")
    counts = {}
    out = []
    prot_resseqs = []
    for ln in open(pin):
        if ln.startswith(("ATOM", "HETATM")):
            resname = ln[17:20].strip()
            resseq = int(ln[22:26])
            ch = chain_for(resname, resseq)
            ln = ln[:21] + ch + ln[22:]
            if resname in PROT:
                prot_resseqs.append(resseq)
            if resname in PROT or resname in ("NPH", "HBI"):
                counts[ch] = counts.get(ch, 0) + 1
        out.append(ln)
    # GUARD (revisión Codex): la aritmética resSeq//288 asume numeración 0-based de
    # mdtraj/AMBER. Si el PDB viniera renumerado 1-based, todos los límites de cadena
    # se desplazarían y la asignación sería errónea de forma silenciosa -> abortar.
    lo, hi = min(prot_resseqs), max(prot_resseqs)
    if lo != 0 or hi != 4 * RESPP - 1:
        raise SystemExit(f"[ERROR] proteína con resSeq {lo}..{hi}; se esperaba 0..{4*RESPP-1} "
                         f"(4×{RESPP}, 0-based). El PDB no calza con la regla //{RESPP}; "
                         f"revisa la numeración antes de reasignar cadenas.")
    with open(pout, "w") as f:
        f.writelines(out)
    print(f"[OK] {pout}")
    print("átomos (proteína+ligandos) por chain:",
          {k: counts[k] for k in sorted(counts)})


if __name__ == "__main__":
    main()
