#!/usr/bin/env python
"""PASO 1 — Preparar y equilibrar el sistema.

Toma la estructura experimental (PDB), la limpia, la mete en una caja de agua
con iones, minimiza la energía y la equilibra. Deja todo listo para correr
trayectorias. Solo se ejecuta UNA vez.

Uso:
    python src/prepare_system.py                 # usa config.yaml
    python src/prepare_system.py --pdb data/otra.pdb

Salidas (en results/prep/):
    system.xml        -> el System de OpenMM (campo de fuerzas + caja)
    system.pdb        -> topología + posiciones (referencia para análisis)
    equilibrated.xml  -> estado equilibrado (punto de partida del muestreo)
"""
import argparse
import os

import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer

from common import load_config, make_integrator, make_simulation, steps_for_ns


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--pdb", default=None, help="estructura de entrada (PDB)")
    args = ap.parse_args()

    cfg = load_config(args.config)
    pdb_in = args.pdb or cfg["input_pdb"]
    outdir = cfg["prep_dir"]
    os.makedirs(outdir, exist_ok=True)

    # --- 1. Limpiar/reparar la estructura con PDBFixer ---
    print(f"[1/6] Limpiando estructura {pdb_in} ...")
    fixer = PDBFixer(filename=pdb_in)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)   # quita ligandos/aguas cristalográficas
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(cfg.get("ph", 7.0))

    # --- 2. Campo de fuerzas + solvatación (agua + iones) ---
    print("[2/6] Añadiendo agua e iones ...")
    ff = app.ForceField(*cfg["forcefield"])
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.addSolvent(
        ff,
        model=cfg["water_model"],
        padding=cfg["solvent_padding_nm"] * unit.nanometer,
        ionicStrength=cfg["ionic_strength_M"] * unit.molar,
        neutralize=True,
    )
    n_atoms = modeller.topology.getNumAtoms()
    print(f"      Sistema solvatado: {n_atoms} átomos")

    # --- 3. Crear el System (PME + restricciones de H + HMR) ---
    print("[3/6] Construyendo el sistema ...")
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=cfg["nonbonded_cutoff_nm"] * unit.nanometer,
        constraints=app.HBonds,
        hydrogenMass=cfg["hydrogen_mass_amu"] * unit.amu,
    )
    # barostato -> ensemble NPT (presión y temperatura constantes)
    system.addForce(mm.MonteCarloBarostat(
        cfg["pressure_bar"] * unit.bar, cfg["temperature_K"] * unit.kelvin))

    integrator = make_integrator(cfg)
    sim = make_simulation(modeller.topology, system, integrator, cfg)
    sim.context.setPositions(modeller.positions)

    # --- 4. Minimización de energía ---
    print("[4/6] Minimizando energía ...")
    sim.minimizeEnergy()

    # --- 5. Equilibración ---
    eq_steps = steps_for_ns(cfg["equilibration_ns"], cfg)
    print(f"[5/6] Equilibrando {cfg['equilibration_ns']} ns ({eq_steps} pasos) ...")
    sim.context.setVelocitiesToTemperature(cfg["temperature_K"] * unit.kelvin)
    sim.reporters.append(app.StateDataReporter(
        os.path.join(outdir, "equilibration.log"), 1000,
        step=True, time=True, temperature=True, potentialEnergy=True,
        progress=True, totalSteps=eq_steps))
    sim.step(eq_steps)

    # --- 6. Guardar todo ---
    print("[6/6] Guardando system.xml, system.pdb y equilibrated.xml ...")
    with open(os.path.join(outdir, "system.xml"), "w") as f:
        f.write(mm.XmlSerializer.serialize(system))
    state = sim.context.getState(getPositions=True)
    with open(os.path.join(outdir, "system.pdb"), "w") as f:
        app.PDBFile.writeFile(sim.topology, state.getPositions(), f)
    sim.saveState(os.path.join(outdir, "equilibrated.xml"))

    print(f"\nListo. Sistema preparado en '{outdir}/'. "
          "Ahora corre: python src/adaptive_sampling.py")


if __name__ == "__main__":
    main()
