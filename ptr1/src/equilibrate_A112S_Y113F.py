# -*- coding: utf-8 -*-
"""
Minimización + equilibración escalonada del sistema PTR1 (tetrámero + 4 NADPH +
4 HBI) en OpenMM/CUDA. Protocolo conservador para un MODELO ENSAMBLADO, no por homología
(monómero de AlphaFold DB montado sobre el cristal 1E92; 75% id.)
con contactos de interfaz ajustados y loops AF no validados por cristal
(recomendación de la verificación cruzada con Codex).

Etapas:
  1. Minimización con soluto fuertemente restringido (relaja solvente/iones).
  2. Minimización con backbone+ligandos moderadamente restringidos.
  3. NVT: calentamiento 0 -> 300 K (300 ps), restraints en backbone+ligandos.
  4. NPT: 1 ns, bajando restraints de backbone+ligandos por etapas.
  5. NPT: 2 ns, restraints solo en Cα.
  6. NPT libre: 5 ns (piloto). Go/no-go: RMSD por cadena, interfaces, ligandos.

Entrada:  results/system/system.{prmtop,inpcrd}
Salida:   results/system/equil/  (estados .xml, trayectoria .dcd, log .csv)

Uso (env jic-folding):  python src/equilibrate_system.py
"""
import os
import sys
from openmm import app
import openmm as mm
from openmm import unit

BASE = "."
SYS = os.path.join(BASE, "results", "system")
OUT = os.path.join(SYS, "equil_A112S_Y113F")
os.makedirs(OUT, exist_ok=True)

TEMP = 300.0 * unit.kelvin
PRESSURE = 1.0 * unit.bar
TIMESTEP = 2.0 * unit.femtoseconds          # conservador (sin HMR) para equilibrar
FRICTION = 1.0 / unit.picosecond
CUTOFF = 1.0 * unit.nanometer

def kcalA2(k):
    """k en kcal/mol/Å² -> kJ/mol/nm²."""
    return k * 418.4

def make_restraint(system, positions, sel_indices, k_kcal):
    """Añade un restraint armónico posicional a los átomos sel_indices.
    Devuelve el índice de la fuerza y su Global parameter 'k' para modularlo."""
    force = mm.CustomExternalForce("0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
    force.addGlobalParameter("k", kcalA2(k_kcal) * unit.kilojoule_per_mole / unit.nanometer**2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    for i in sel_indices:
        x, y, z = positions[i].value_in_unit(unit.nanometer)
        force.addParticle(i, [x, y, z])
    fidx = system.addForce(force)
    return fidx, force

def select(topology, kind):
    """Índices de átomos por categoría."""
    solute_res = {"WAT", "HOH", "Na+", "Cl-", "NA", "CL"}
    idx = []
    for atom in topology.atoms():
        rn = atom.residue.name
        is_solute = rn not in solute_res
        if kind == "solute_heavy" and is_solute and atom.element != app.element.hydrogen:
            idx.append(atom.index)
        elif kind == "backbone_lig" and is_solute and (
                atom.name in ("CA", "C", "N", "O") or atom.residue.name in ("NPH", "HBI")):
            if atom.element != app.element.hydrogen:
                idx.append(atom.index)
        elif kind == "calpha" and atom.name == "CA":
            idx.append(atom.index)
    return idx

def build_system():
    prm = app.AmberPrmtopFile(os.path.join(SYS, "mutant_A112S_Y113F.prmtop"))
    inp = app.AmberInpcrdFile(os.path.join(SYS, "mutant_A112S_Y113F.inpcrd"))
    system = prm.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=CUTOFF,
                              constraints=app.HBonds, rigidWater=True)
    return prm, inp, system

def make_sim(prm, system, positions, boxvectors, platform_name="CUDA"):
    integrator = mm.LangevinMiddleIntegrator(TEMP, FRICTION, TIMESTEP)
    platform = mm.Platform.getPlatformByName(platform_name)
    props = {"Precision": "mixed"} if platform_name in ("CUDA", "OpenCL") else {}
    sim = app.Simulation(prm.topology, system, integrator, platform, props)
    sim.context.setPositions(positions)
    if boxvectors is not None:
        sim.context.setPeriodicBoxVectors(*boxvectors)
    return sim

def report_state(sim, tag):
    st = sim.context.getState(getEnergy=True)
    pe = st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    print(f"  [{tag}] E_pot = {pe:.1f} kcal/mol", flush=True)

def main():
    print("== Construyendo sistema ==", flush=True)
    prm, inp, system = build_system()
    pos0 = inp.getPositions()
    box = inp.boxVectors

    # restraint sobre backbone (Cα,C,N,O) + ligandos pesados, k modulable.
    # Las cadenas laterales quedan libres desde el inicio -> relajan los clashes
    # del modelo ensamblado; el backbone protege el plegamiento.
    fidx, force = make_restraint(system, pos0, select(prm.topology, "backbone_lig"), 100.0)

    # plataforma con fallback
    plat = "CUDA"
    try:
        sim = make_sim(prm, system, pos0, box, plat)
    except Exception as e:
        print(f"  CUDA no disponible ({e}); usando CPU", flush=True)
        plat = "CPU"; sim = make_sim(prm, system, pos0, box, plat)
    print(f"== Plataforma: {plat} ==", flush=True)

    # --- Etapa 1: min con backbone+lig fuertemente restringido k=100 ---
    print("== Etapa 1: minimización (backbone+lig k=100; libera solvente+sidechains) ==", flush=True)
    report_state(sim, "pre-min1")
    sim.minimizeEnergy(maxIterations=5000)
    report_state(sim, "post-min1")

    # --- Etapa 2: min con restraint suavizado a k=10 ---
    print("== Etapa 2: minimización (backbone+lig k=10) ==", flush=True)
    sim.context.setParameter("k", kcalA2(10.0))
    sim.minimizeEnergy(maxIterations=5000)
    report_state(sim, "post-min2")

    # --- Etapa 3: NVT calentamiento 0->300 K, restraint k=10 ---
    print("== Etapa 3: NVT calentamiento 0->300 K (300 ps, k=10) ==", flush=True)
    sim.context.setVelocitiesToTemperature(5 * unit.kelvin)
    n_steps = int(300 * unit.picoseconds / TIMESTEP)
    n_chunks = 30
    for c in range(n_chunks):
        T = (5 + (300 - 5) * (c + 1) / n_chunks) * unit.kelvin
        sim.integrator.setTemperature(T)
        sim.step(n_steps // n_chunks)
    report_state(sim, "post-NVT")

    # añadir barostato para NPT
    system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMP))
    sim.context.reinitialize(preserveState=True)

    # --- Etapa 4: NPT 1 ns, bajando restraint 10 -> 2 ---
    print("== Etapa 4: NPT 1 ns (restraint 10 -> 2) ==", flush=True)
    sim.reporters.append(app.StateDataReporter(
        os.path.join(OUT, "equil.csv"), 5000, step=True, time=True,
        potentialEnergy=True, temperature=True, density=True, volume=True, speed=True))
    sim.reporters.append(app.DCDReporter(os.path.join(OUT, "equil.dcd"), 10000))
    for k in (8.0, 6.0, 4.0, 2.0):
        sim.context.setParameter("k", kcalA2(k))
        sim.step(int(250 * unit.picoseconds / TIMESTEP))
        report_state(sim, f"NPT k={k}")

    # --- Etapa 5: NPT 2 ns, restraint solo simbólico (k=1) ---
    print("== Etapa 5: NPT 2 ns (k=1) ==", flush=True)
    sim.context.setParameter("k", kcalA2(1.0))
    sim.step(int(2000 * unit.picoseconds / TIMESTEP))
    report_state(sim, "post-NPT-restr")
    sim.saveState(os.path.join(OUT, "equil_restrained.xml"))

    # --- Etapa 6: NPT libre 5 ns (piloto) ---
    print("== Etapa 6: NPT libre 5 ns (piloto) ==", flush=True)
    sim.context.setParameter("k", 0.0)
    sim.step(int(5000 * unit.picoseconds / TIMESTEP))
    report_state(sim, "post-piloto")
    sim.saveState(os.path.join(OUT, "equil_free.xml"))
    with open(os.path.join(OUT, "system_final.pdb"), "w") as fh:
        st = sim.context.getState(getPositions=True)
        app.PDBFile.writeFile(prm.topology, st.getPositions(), fh)
    print("== Equilibración + piloto COMPLETOS ==", flush=True)
    print(f"   Estados: {OUT}/equil_free.xml ; trayectoria: equil.dcd ; log: equil.csv", flush=True)

if __name__ == "__main__":
    main()
