#!/usr/bin/env bash
# ============================================================================
#  Parametrización del cofactor NADPH (carga -4) para MD clásica del tetrámero
#  PTR1 de L. panamensis. Reproducible, env jic-folding (AmberTools 24.8).
#
#  Estrategia (decidida con verificación cruzada de Codex):
#    GAFF2 (parámetros completos por analogía, vía antechamber/parmchk2)
#    + cargas RESP publicadas de Ryde (Lund) para NADPH.
#  Evita el AM1-BCC frágil en una molécula grande y cargada (-4), y no toca los
#  tipos estándar de Amber (no contamina otros residuos del sistema).
#
#  Requisito previo: haber corrido  python src/build_nadph.py  (mapea el NAP
#  cristalográfico de 1E92 a la nomenclatura NPH de Ryde, verificado por
#  conectividad) -> results/ligands/nadph_NPH_chainA.pdb
#
#  Uso:  bash src/parametrize_nadph.sh
# ============================================================================
set -euo pipefail

BASE="$(cd "$(dirname "$0")/.." && pwd)"
LIG="$BASE/results/ligands"
FF="$LIG/ff_params"
TLE="$LIG/tleap"
mkdir -p "$TLE" "$FF"

# --- Paso A: tleap con el prep de Ryde -> mol2 (cargas RESP) + PDB con H ------
# Reconstruye los 26 H (incluido el hidruro de reducción en C4N) sobre las
# coordenadas pesadas del cristal y fija la carga -4 desde el CHARGE block.
cat > "$TLE/build_nadph.in" <<EOF
source leaprc.protein.ff14SB
source leaprc.RNA.OL3
loadAmberPrep $FF/nadph.prep
nph = loadPdb $LIG/nadph_NPH_chainA.pdb
saveMol2 nph $TLE/nadph.mol2 1
savePdb  nph $TLE/nadph_H.pdb
quit
EOF
tleap -f "$TLE/build_nadph.in" > "$TLE/build_nadph.log" 2>&1
echo "[A] mol2 con cargas RESP + PDB con H generados."

# --- Paso B: extraer cargas RESP del mol2 ------------------------------------
awk '/@<TRIPOS>ATOM/{f=1;next}/@<TRIPOS>BOND/{f=0}f{print $9}' \
    "$TLE/nadph.mol2" > "$FF/ryde.chg"

# --- Paso C: re-tipar a GAFF2 conservando las cargas RESP --------------------
cd "$FF"
antechamber -i "$TLE/nadph.mol2" -fi mol2 -o NPH.gaff2.mol2 -fo mol2 \
    -at gaff2 -c rc -cf ryde.chg -nc -4 -rn NPH -dr no > antechamber.log 2>&1
echo "[C] NADPH re-tipado a GAFF2 (cargas RESP conservadas)."

# --- Paso D: completar parámetros faltantes por analogía GAFF2 ---------------
parmchk2 -i NPH.gaff2.mol2 -f mol2 -o NPH.gaff2.frcmod -s gaff2
if grep -q "ATTN" NPH.gaff2.frcmod; then
    echo "[D] AVISO: parmchk2 marcó parámetros con ATTN — revisar NPH.gaff2.frcmod" >&2
else
    echo "[D] frcmod GAFF2 generado sin ATTN (todos los parámetros cubiertos)."
fi

# --- Paso E: topología del cofactor aislado (verificación) -------------------
cat > leap_NPH_gaff2.in <<EOF
source leaprc.protein.ff14SB
source leaprc.RNA.OL3
source leaprc.gaff2
loadamberparams NPH.gaff2.frcmod
NPH = loadmol2 NPH.gaff2.mol2
check NPH
charge NPH
saveamberparm NPH NPH.gaff2.prmtop NPH.gaff2.inpcrd
savepdb NPH NPH.gaff2.pdb
quit
EOF
tleap -f leap_NPH_gaff2.in > tleap_NPH_gaff2.log 2>&1
ERR=$(grep -oE "Errors = [0-9]+" tleap_NPH_gaff2.log | tail -1)
echo "[E] tleap topología del cofactor: $ERR (la advertencia de carga -4 es esperada)."

echo ""
echo "Artefactos para el ensamblaje del sistema completo:"
echo "  $FF/NPH.gaff2.mol2     (NADPH retipado GAFF2 + cargas RESP, residuo NPH)"
echo "  $FF/NPH.gaff2.frcmod   (parámetros GAFF2 del cofactor)"
echo "  $FF/NPH.gaff2.prmtop   (topología del cofactor aislado, para verificación)"
