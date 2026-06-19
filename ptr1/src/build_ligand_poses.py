# -*- coding: utf-8 -*-
"""
Genera los 4 NADPH (residuo NPH) y 4 HBI en sus poses cristalográficas (una por
cadena A/B/C/D) a partir del tetrámero holo, para el ensamblaje del sistema.

- NADPH: toma el NAP cristalográfico de cada cadena y lo renombra a la nomenclatura
  NPH de Ryde (mismo MAP verificado por conectividad en build_nadph.py). Solo átomos
  pesados; tleap reconstruye los H desde la lib GAFF2 (NPH.lib) en cada pose.
- HBI: toma el HBI cristalográfico de cada cadena (ya usa nombres del CCD, idénticos
  a la lib HBI.gaff2). Solo átomos pesados; tleap reconstruye los H.

Salida: results/system/nph_chain{A..D}.pdb y hbi_chain{A..D}.pdb

Uso (env jic-folding):  python src/build_ligand_poses.py
"""
import os

BASE = "."
TET = os.path.join(BASE, "results", "tetramer", "ptr1_Lpanamensis_tetramer_holo.pdb")
OUT = os.path.join(BASE, "results", "system")
CHAINS = ["A", "B", "C", "D"]

# Mapeo NAP(cristal) -> NPH(Ryde), verificado por conectividad (ver build_nadph.py)
MAP_NAP = {
    "N1N": "N1N", "C2N": "C2N", "C3N": "C3N", "C4N": "C4N", "C5N": "C5N",
    "C6N": "C6N", "C7N": "C7N", "O7N": "O7N", "N7N": "N7N",
    "N1A": "N1A", "C2A": "C2A", "N3A": "N3A", "C4A": "C4A", "C5A": "C5A",
    "C6A": "C6A", "N6A": "N6A", "N7A": "N7A", "C8A": "C8A", "N9A": "N9A",
    "C1D": "C'N1", "C2D": "C'N2", "O2D": "O'N2", "C3D": "C'N3", "O3D": "O'N3",
    "C4D": "C'N4", "O4D": "O'N4", "C5D": "C'N5", "O5D": "O'N5",
    "C1B": "C'A1", "C2B": "C'A2", "O2B": "O'A2", "C3B": "C'A3", "O3B": "O'A3",
    "C4B": "C'A4", "O4B": "O'A4", "C5B": "C'A5", "O5B": "O'A5",
    "PN": "PN", "O1N": "OPN1", "O2N": "OPN2", "O3": "O3P",
    "PA": "PA", "O1A": "OPA1", "O2A": "OPA2",
    "P2B": "P'A2", "O1X": "OA22", "O2X": "OA23", "O3X": "OA24",
}

def fmt_name(n):
    return f"{n:<4s}" if len(n) >= 4 else f" {n:<3s}"

def read_lig(chain, resn):
    """Devuelve líneas (name, element, x,y,z) del ligando resn en la cadena dada."""
    out = []
    for ln in open(TET):
        if ln[:6].strip() not in ("ATOM", "HETATM"):
            continue
        if ln[21] != chain or ln[17:20].strip() != resn:
            continue
        name = ln[12:16].strip()
        el = ln[76:78].strip() or name[0]
        x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
        out.append((name, el, x, y, z))
    return out

def write_pdb(path, atoms, resn, chain, rename=None):
    with open(path, "w") as fh:
        for i, (name, el, x, y, z) in enumerate(atoms, 1):
            nm = rename[name] if rename else name
            fh.write(f"HETATM{i:5d} {fmt_name(nm)} {resn} {chain}   1    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}\n")
        fh.write("END\n")

for ch in CHAINS:
    nap = read_lig(ch, "NAP")
    hbi = read_lig(ch, "HBI")
    assert len(nap) == 48, f"NAP cadena {ch}: {len(nap)} átomos (esperaba 48)"
    assert len(hbi) == 17, f"HBI cadena {ch}: {len(hbi)} átomos (esperaba 17)"
    write_pdb(os.path.join(OUT, f"nph_chain{ch}.pdb"), nap, "NPH", ch, rename=MAP_NAP)
    write_pdb(os.path.join(OUT, f"hbi_chain{ch}.pdb"), hbi, "HBI", ch)
    print(f"cadena {ch}: NPH ({len(nap)} pesados) + HBI ({len(hbi)} pesados) escritos")

print("\nLigandos en pose listos para el ensamblaje en tleap.")
