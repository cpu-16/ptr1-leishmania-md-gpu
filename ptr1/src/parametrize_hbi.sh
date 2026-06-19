#!/usr/bin/env bash
# ============================================================================
#  Parametrización del sustrato HBI (7,8-dihidrobiopterina) para MD clásica del
#  tetrámero PTR1 de L. panamensis. Reproducible, env jic-folding (AmberTools 24.8).
#
#  Estrategia: GAFF2 + cargas AM1-BCC. A diferencia del NADPH, aquí AM1-BCC SÍ es
#  fiable: molécula pequeña (30 átomos), neutra, sin fosfatos. Se parte de la
#  definición ideal del PDB Chemical Component Dictionary (CCD), que trae los H y
#  los órdenes de enlace correctos y usa los MISMOS nombres de átomo pesado que el
#  cristal 1E92 (N1, C2, N2, ...), lo que simplifica el ensamblaje posterior.
#
#  Requisito previo: HBI.cif descargado del CCD
#    (https://files.rcsb.org/ligands/download/HBI.cif -> results/ligands/ff_params/)
#
#  NOTA: la carga total del prmtop queda en ~0.002 e (redondeo de AM1-BCC en el
#  mol2 a 6 decimales). Es despreciable; se absorbe en la neutralización del
#  sistema completo con PME. Documentado en ESTADO_RAMA_B.md.
#
#  Uso:  bash src/parametrize_hbi.sh
# ============================================================================
set -euo pipefail

BASE="$(cd "$(dirname "$0")/.." && pwd)"
FF="$BASE/results/ligands/ff_params"
cd "$FF"

# --- Paso A: CIF del CCD -> PDB ideal con H + nombres CCD --------------------
# (el parser ccif de antechamber 24 aborta con HBI; parseamos el loop a mano.)
python3 - <<'PY'
out=[]; serial=1
for ln in open("HBI.cif"):
    p=ln.split()
    if len(p)>=18 and p[0]=="HBI" and p[17].isdigit():
        name=p[1]; el=p[3]
        x,y,z=float(p[12]),float(p[13]),float(p[14])   # coords ideales del CCD
        nm=f"{name:<4s}" if len(name)>=4 else f" {name:<3s}"
        out.append(f"HETATM{serial:5d} {nm} HBI A   1    "
                   f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}")
        serial+=1
open("HBI_ideal.pdb","w").write("\n".join(out)+"\nEND\n")
print(f"[A] PDB ideal: {len(out)} átomos")
PY

# --- Paso B: antechamber AM1-BCC (carga neta 0) -----------------------------
antechamber -i HBI_ideal.pdb -fi pdb -o HBI.gaff2.mol2 -fo mol2 \
  -at gaff2 -c bcc -nc 0 -rn HBI -dr no > antechamber_hbi.log 2>&1
echo "[B] HBI parametrizado con GAFF2 + AM1-BCC."

# --- Paso C: completar parámetros por analogía GAFF2 ------------------------
parmchk2 -i HBI.gaff2.mol2 -f mol2 -o HBI.gaff2.frcmod -s gaff2
if grep -q "ATTN" HBI.gaff2.frcmod; then
    echo "[C] AVISO: parmchk2 marcó ATTN — revisar HBI.gaff2.frcmod" >&2
else
    echo "[C] frcmod GAFF2 sin ATTN (todos los parámetros cubiertos)."
fi

# --- Paso D: topología del sustrato aislado (verificación) ------------------
cat > leap_HBI_gaff2.in <<'EOF'
source leaprc.gaff2
loadamberparams HBI.gaff2.frcmod
HBI = loadmol2 HBI.gaff2.mol2
check HBI
charge HBI
saveamberparm HBI HBI.gaff2.prmtop HBI.gaff2.inpcrd
savepdb HBI HBI.gaff2.pdb
quit
EOF
tleap -f leap_HBI_gaff2.in > tleap_HBI_gaff2.log 2>&1
echo "[D] tleap topología del sustrato: $(grep -oE 'Errors = [0-9]+' tleap_HBI_gaff2.log | tail -1)"

echo ""
echo "Artefactos para el ensamblaje del sistema:"
echo "  $FF/HBI.gaff2.mol2     (HBI GAFF2 + AM1-BCC, residuo HBI, nombres CCD)"
echo "  $FF/HBI.gaff2.frcmod   (parámetros GAFF2 del sustrato)"
