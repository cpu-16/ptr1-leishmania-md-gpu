# -*- coding: utf-8 -*-
"""
Construye el ligando NADPH (residuo NPH de Ryde) en la pose cristalográfica del
NADP+ (NAP) de 1E92, renombrando los átomos pesados a la nomenclatura de Ryde.

Por qué: los parámetros publicados de Ryde (cargas RESP, results/ligands/ff_params/
nadph.prep) usan el residuo NPH con SUS nombres de átomo. El NAP cristalográfico usa
la nomenclatura del PDB (C1B/C1D/P2B...). Este script:
  1. Mapea cristal -> Ryde por correspondencia química (anillos coinciden; ribosas y
     fosfatos se renombran).
  2. VERIFICA el mapeo por conectividad (distancias en el cristal) para no confundir la
     ribosa de nicotinamida con la de adenina.
  3. Escribe un PDB con residuo NPH y nombres Ryde, en las coordenadas del cristal.

Luego tleap (loadAmberPrep nadph.prep) reconstruye los hidrógenos — incluido el
hidruro extra H4 que reduce NADP+ a NADPH — sobre estas coordenadas pesadas.

NOTA: el residuo NPH tiene carga total -4 (ambos fosfatos desprotonados). El 2'-fosfato
tiene pKa 6-7; -4 es la forma estándar a pH fisiológico. (Variante -3 = residuo NP3.)

Uso (env jic-folding):  python src/build_nadph.py
"""
import os
import numpy as np

BASE = "."
NAP_PDB = os.path.join(BASE, "results", "ligands", "nap_NADPplus_chainA.pdb")
OUT_PDB = os.path.join(BASE, "results", "ligands", "nadph_NPH_chainA.pdb")

# --- Mapeo químico cristal(NAP) -> Ryde(NPH) -----------------------------------
MAP = {
    # Anillo nicotinamida (nombres idénticos)
    "N1N": "N1N", "C2N": "C2N", "C3N": "C3N", "C4N": "C4N", "C5N": "C5N",
    "C6N": "C6N", "C7N": "C7N", "O7N": "O7N", "N7N": "N7N",
    # Anillo adenina (nombres idénticos)
    "N1A": "N1A", "C2A": "C2A", "N3A": "N3A", "C4A": "C4A", "C5A": "C5A",
    "C6A": "C6A", "N6A": "N6A", "N7A": "N7A", "C8A": "C8A", "N9A": "N9A",
    # Ribosa de la nicotinamida (sufijo D -> 'N) — verificado por conectividad: C1D-N1N
    "C1D": "C'N1", "C2D": "C'N2", "O2D": "O'N2", "C3D": "C'N3", "O3D": "O'N3",
    "C4D": "C'N4", "O4D": "O'N4", "C5D": "C'N5", "O5D": "O'N5",
    # Ribosa de la adenina (sufijo B -> 'A) — verificado por conectividad: C1B-N9A
    "C1B": "C'A1", "C2B": "C'A2", "O2B": "O'A2", "C3B": "C'A3", "O3B": "O'A3",
    "C4B": "C'A4", "O4B": "O'A4", "C5B": "C'A5", "O5B": "O'A5",
    # Fosfato lado nicotinamida
    "PN": "PN", "O1N": "OPN1", "O2N": "OPN2",
    # Oxígeno puente del pirofosfato
    "O3": "O3P",
    # Fosfato lado adenina
    "PA": "PA", "O1A": "OPA1", "O2A": "OPA2",
    # 2'-fosfato (en la ribosa de adenina)
    "P2B": "P'A2", "O1X": "OA22", "O2X": "OA23", "O3X": "OA24",
}

# Enlaces clave que DEBEN cumplirse en el cristal para validar la asignación de ribosas
# (nombres cristalográficos). Si alguno falla, el mapeo de sufijos B/D estaría invertido.
KEY_BONDS = [
    ("C1D", "N1N"),   # anomérico ribosa-nicotinamida (sufijo D) unido a la nicotinamida
    ("C1B", "N9A"),   # anomérico ribosa-adenina (sufijo B) unido a la adenina
    ("O5D", "PN"),    # 5'-O de ribosa nicotinamida esterifica el fosfato PN
    ("O5B", "PA"),    # 5'-O de ribosa adenina esterifica el fosfato PA
    ("O2B", "P2B"),   # 2'-O de ribosa adenina esterifica el 2'-fosfato
    ("O3", "PA"),     # oxígeno puente del pirofosfato
    ("O3", "PN"),     # oxígeno puente del pirofosfato
    ("O4D", "C1D"),   # anillo furanosa nicotinamida
    ("O4B", "C1B"),   # anillo furanosa adenina
]

def load_atoms(path):
    at = {}  # name -> (element, xyz, raw_line)
    for ln in open(path):
        if ln[:6].strip() not in ("ATOM", "HETATM"):
            continue
        name = ln[12:16].strip()
        el = ln[76:78].strip() or name[0]
        xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        at[name] = (el, xyz, ln)
    return at

atoms = load_atoms(NAP_PDB)
print(f"Átomos leídos del NAP cristalográfico: {len(atoms)}")

# --- 1) cobertura: todos los nombres del cristal están en el mapeo ---
faltan = [n for n in atoms if n not in MAP]
sobran = [k for k in MAP if k not in atoms]
assert not faltan, f"Átomos del cristal sin mapear: {faltan}"
assert not sobran, f"Mapeo con átomos inexistentes en el cristal: {sobran}"
assert len(set(MAP.values())) == len(MAP), "Colisión: dos átomos mapean al mismo nombre Ryde"
print(f"Cobertura OK: {len(MAP)} átomos mapeados 1:1, sin colisiones.")

# --- 2) consistencia de elementos (primera letra del nombre Ryde = elemento) ---
for crys, ryde in MAP.items():
    el_crys = atoms[crys][0]
    el_ryde = ryde.lstrip("'")[0]  # C'N1 -> C
    assert el_crys == el_ryde, f"Elemento no coincide: {crys}({el_crys}) -> {ryde}({el_ryde})"
print("Elementos consistentes en los 48 pares.")

# --- 3) verificación por conectividad (distancias en el cristal) ---
def dist(a, b):
    return float(np.linalg.norm(atoms[a][1] - atoms[b][1]))

print("\nVerificación de enlaces clave (cristal):")
ok = True
for a, b in KEY_BONDS:
    d = dist(a, b)
    bonded = d < 1.95
    flag = "OK " if bonded else "FALLA"
    if not bonded:
        ok = False
    print(f"  {flag}  {a:4s}-{b:4s} = {d:.2f} Å")
assert ok, "Algún enlace clave NO se cumple: revisar el mapeo de ribosas/fosfatos."
print("Conectividad consistente: asignación D(nicotinamida)/B(adenina) correcta.")

# --- 4) escribir PDB con residuo NPH y nombres Ryde ---
def fmt_name(n):
    # PDB columnas 13-16: nombres de <=3 chars empiezan en col 14; 4 chars en col 13.
    return f"{n:<4s}" if len(n) >= 4 else f" {n:<3s}"

with open(OUT_PDB, "w") as fh:
    serial = 1
    for crys, (el, xyz, raw) in atoms.items():
        ryde = MAP[crys]
        # Formato PDB estricto: nombre cols 13-16, altLoc col 17 (blanco),
        # resName cols 18-20 ("NPH"), chain col 22, resSeq cols 23-26, x desde col 31.
        line = (
            f"HETATM{serial:5d} {fmt_name(ryde)} NPH A   1    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {el:>2s}\n"
        )
        fh.write(line)
        serial += 1
    fh.write("END\n")

print(f"\nEscrito: {OUT_PDB}  ({len(atoms)} átomos pesados, residuo NPH)")
print("Siguiente: tleap reconstruye los H (incluido el hidruro de reducción) desde nadph.prep.")
