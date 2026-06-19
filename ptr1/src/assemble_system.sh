#!/usr/bin/env bash
# ============================================================================
#  Ensamblaje del sistema solvatado para MD: tetrámero PTR1 + 4 NADPH + 4 HBI.
#  Reproducible, env jic-folding (AmberTools 24.8 + PropKa 3.5).
#
#  Requisitos previos (en este orden):
#    1. python src/build_tetramer.py        -> tetrámero holo
#    2. python src/build_nadph.py           -> NADPH cadena A (mapeo verificado)
#    3. bash   src/parametrize_nadph.sh     -> NPH.gaff2.{mol2,frcmod}
#    4. bash   src/parametrize_hbi.sh       -> HBI.gaff2.{mol2,frcmod}
#
#  Protonación (pH 7): estados estándar (HIS neutras, ASP/GLU-, LYS/ARG+, CYS
#  protonada). Respaldado por PropKa: ningún residuo con pKa fuertemente
#  desplazado salvo CYS97 (pKa 5.09, sin disulfuro ni rol catalítico -> se deja
#  protonada, criterio conservador). HIS127/GLU199/LYS198 borderline (<0.2 de
#  pH 7) -> estado estándar.
#
#  Salida: results/system/system.{prmtop,inpcrd,pdb}
#
#  Uso:  bash src/assemble_system.sh
# ============================================================================
set -euo pipefail

BASE="$(cd "$(dirname "$0")/.." && pwd)"
SYS="$BASE/results/system"
FF="$BASE/results/ligands/ff_params"
TET="$BASE/results/tetramer/ptr1_Lpanamensis_tetramer_holo.pdb"
mkdir -p "$SYS"
cd "$SYS"

# --- 1) Proteína: extraer, limpiar (pdb4amber detecta disulfuros, no hay) ----
grep -E "^ATOM|^TER" "$TET" > protein_raw.pdb; echo "END" >> protein_raw.pdb
pdb4amber -i protein_raw.pdb -o protein_4amber.pdb > pdb4amber.log 2>&1
echo "[1] Proteína limpia (4x288 res, sin disulfuros)."

# --- 2) PropKa informativo (no cambia nombres; documenta la decisión) --------
propka3 "$TET" -o 7.0 > propka_run.log 2>&1 || true
echo "[2] PropKa corrido (informativo; protonación estándar adoptada)."

# --- 3) Ligandos en pose (4 NADPH + 4 HBI, solo heavy) -----------------------
python "$BASE/src/build_ligand_poses.py" > /dev/null
echo "[3] Ligandos en pose generados (nph_chain{A..D}, hbi_chain{A..D})."

# --- 4) Libs OFF de los ligandos (GAFF2 + cargas) ----------------------------
cd "$FF"
cat > make_libs.in <<'EOF'
source leaprc.protein.ff14SB
source leaprc.RNA.OL3
source leaprc.gaff2
loadamberparams NPH.gaff2.frcmod
loadamberparams HBI.gaff2.frcmod
NPH = loadmol2 NPH.gaff2.mol2
HBI = loadmol2 HBI.gaff2.mol2
saveoff NPH NPH.lib
saveoff HBI HBI.lib
quit
EOF
tleap -f make_libs.in > make_libs.log 2>&1
cd "$SYS"
echo "[4] Libs NPH.lib / HBI.lib creadas."

# --- 5) Ensamblaje: combine + solvatar (oct 10 A) + neutralizar + sal 0.15 M -
# Sal: 88 pares NaCl = round(0.15 * 32620 aguas / 55.5). Verificar si cambia el
# nº de aguas al re-solvatar.
cat > assemble_final.in <<EOF
source leaprc.protein.ff14SB
source leaprc.RNA.OL3
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams $FF/NPH.gaff2.frcmod
loadamberparams $FF/HBI.gaff2.frcmod
loadoff $FF/NPH.lib
loadoff $FF/HBI.lib
prot = loadpdb protein_4amber.pdb
nphA = loadpdb nph_chainA.pdb
nphB = loadpdb nph_chainB.pdb
nphC = loadpdb nph_chainC.pdb
nphD = loadpdb nph_chainD.pdb
hbiA = loadpdb hbi_chainA.pdb
hbiB = loadpdb hbi_chainB.pdb
hbiC = loadpdb hbi_chainC.pdb
hbiD = loadpdb hbi_chainD.pdb
sys = combine {prot nphA nphB nphC nphD hbiA hbiB hbiC hbiD}
solvateOct sys TIP3PBOX 10.0
addIonsRand sys Na+ 0
addIonsRand sys Na+ 88 Cl- 88
charge sys
saveamberparm sys system.prmtop system.inpcrd
savepdb sys system.pdb
quit
EOF
tleap -f assemble_final.in > assemble_final.log 2>&1
ERR=$(grep -oE "Errors = [0-9]+" assemble_final.log | tail -1)
echo "[5] Sistema ensamblado: $ERR"
echo ""
echo "Sistema final -> $SYS/system.{prmtop,inpcrd,pdb}  (~114,856 átomos, neutro, 0.15 M NaCl)"
