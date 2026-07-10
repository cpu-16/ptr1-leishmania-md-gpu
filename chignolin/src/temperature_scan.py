#!/usr/bin/env python
"""Escaneo de temperatura para encontrar dónde Chignolin se pliega Y se despliega.

A 340 K (NPT) la variante CLN025 quedó 99.5% plegada -> casi no se desplegó ->
no se puede medir la cinética.

⚠️ NO concluir de aquí que el campo de fuerzas sobre-estabiliza el plegado. Esa
afirmación se retractó el 9 jul 2026: trayectorias cortas que arrancan plegadas se
quedan plegadas con casi cualquier campo de fuerzas.
Este script prueba varias temperaturas (en NVT, volumen fijo, para evitar que el
agua "hierva" a T alta), corre una trayectoria corta en cada una y mide la
fracción plegada. Sirve para elegir la T del estudio cinético principal:
la mejor es la más cercana a ~50% plegada (allí hay muchas transiciones).

Uso:
    python src/temperature_scan.py

Es METODOLOGÍA legítima del proyecto (estimar la temperatura de fusión aparente).
"""
import os

import mdtraj as md
import numpy as np
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer

from common import load_config, make_simulation

TEMPS = [360, 380, 400]     # K (NVT, así que no hay riesgo de ebullición)
EQ_NS = 0.2
PROD_NS = 10.0
OUTDIR = "results_scan"


def build_system_nvt(cfg):
    """Construye el sistema SIN barostato (NVT)."""
    fixer = PDBFixer(filename=cfg["input_pdb"])
    fixer.findMissingResidues(); fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues(); fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms(); fixer.addMissingAtoms()
    fixer.addMissingHydrogens(cfg.get("ph", 7.0))
    ff = app.ForceField(*cfg["forcefield"])
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.addSolvent(ff, model=cfg["water_model"],
                        padding=cfg["solvent_padding_nm"] * unit.nanometer,
                        ionicStrength=cfg["ionic_strength_M"] * unit.molar,
                        neutralize=True)
    system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME,
                             nonbondedCutoff=cfg["nonbonded_cutoff_nm"] * unit.nanometer,
                             constraints=app.HBonds,
                             hydrogenMass=cfg["hydrogen_mass_amu"] * unit.amu)
    return modeller, system   # sin MonteCarloBarostat -> NVT


def run_temp(cfg, temp_K, native_ca):
    modeller, system = build_system_nvt(cfg)
    integ = mm.LangevinMiddleIntegrator(temp_K * unit.kelvin,
                                        cfg["friction_per_ps"] / unit.picosecond,
                                        cfg["timestep_ps"] * unit.picoseconds)
    sim = make_simulation(modeller.topology, system, integ, cfg)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy()
    sim.context.setVelocitiesToTemperature(temp_K * unit.kelvin)
    dt = cfg["timestep_ps"]
    sim.step(int(EQ_NS * 1000 / dt))                     # equilibración corta

    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, f"T{temp_K}")
    with open(out + "_top.pdb", "w") as f:
        app.PDBFile.writeFile(
            sim.topology, sim.context.getState(getPositions=True).getPositions(), f)
    sim.reporters.append(app.DCDReporter(out + ".dcd", int(cfg["save_interval_ps"] / dt)))
    sim.step(int(PROD_NS * 1000 / dt))                   # producción

    t = md.load(out + ".dcd", top=out + "_top.pdb")
    t_ca = t.atom_slice(t.topology.select("name CA and protein"))
    rmsd = md.rmsd(t_ca, native_ca)
    folded = float(np.mean(rmsd < cfg["folded_rmsd_nm"]))
    return folded, float(rmsd.min()), float(np.median(rmsd))


def main():
    cfg = load_config()
    native = md.load(cfg["native_pdb"])
    native_ca = native.atom_slice(native.topology.select("name CA and protein"))
    print(f"{'T(K)':>6} {'%plegada':>10} {'RMSDmin(nm)':>12} {'RMSDmed(nm)':>12}")
    res = {}
    for T in TEMPS:
        folded, rmin, rmed = run_temp(cfg, T, native_ca)
        res[T] = folded
        print(f"{T:>6} {folded*100:>9.1f}% {rmin:>12.3f} {rmed:>12.3f}", flush=True)
    best = min(res, key=lambda T: abs(res[T] - 0.5))
    print(f"\nRECOMENDADA_PARA_CINETICA {best} K  (mas cercana a 50% plegada)")


main()
