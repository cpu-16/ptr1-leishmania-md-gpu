# -*- coding: utf-8 -*-
"""
Producción de MD del sistema PTR1 (tetrámero + 4 NADPH + 4 HBI) en OpenMM/CUDA.

Parte del estado equilibrado (results/system/equil/equil_free.xml) y corre NPT a
300 K con HMR (4 amu / paso de 4 fs) para acelerar en la RTX 4060. Resumible: usa
un checkpoint .chk; si existe y se pasa --resume, continúa desde ahí (espíritu de
trabajo independiente/reanudable del proyecto, como run_md.py de Chignolin).

Uso:
  python src/produce_md.py --ns 100                 # 100 ns desde el equilibrado
  python src/produce_md.py --ns 100 --resume        # continuar desde el checkpoint
  python src/produce_md.py --ns 300 --out results/system/prod

Salida: <out>/prod.dcd, prod.csv (log), prod.chk (checkpoint), prod_final.xml
"""
import os
import sys
import argparse
from openmm import app
import openmm as mm
from openmm import unit

BASE = "."
SYS = os.path.join(BASE, "results", "system")

TEMP = 300.0 * unit.kelvin
PRESSURE = 1.0 * unit.bar
FRICTION = 1.0 / unit.picosecond
CUTOFF = 1.0 * unit.nanometer
TIMESTEP = 4.0 * unit.femtoseconds          # HMR permite 4 fs
HMASS = 4.0 * unit.amu
SAVE_PS = 50.0                               # frame cada 50 ps
CHK_PS = 1000.0                              # checkpoint cada 1 ns


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", type=float, default=100.0, help="nanosegundos de producción")
    ap.add_argument("--out", default=os.path.join(SYS, "prod_Y113F"))
    ap.add_argument("--resume", action="store_true", help="continuar desde el checkpoint")
    ap.add_argument("--platform", default="CUDA")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    prm = app.AmberPrmtopFile(os.path.join(SYS, "mutant_Y113F.prmtop"))
    # HMR: redistribuye masa a los H -> integración estable a 4 fs
    system = prm.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=CUTOFF,
                              constraints=app.HBonds, rigidWater=True,
                              hydrogenMass=HMASS)
    system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMP))

    integrator = mm.LangevinMiddleIntegrator(TEMP, FRICTION, TIMESTEP)
    platform = mm.Platform.getPlatformByName(args.platform)
    props = {"Precision": "mixed"} if args.platform in ("CUDA", "OpenCL") else {}
    sim = app.Simulation(prm.topology, system, integrator, platform, props)

    chk = os.path.join(args.out, "prod.chk")
    append = False
    if args.resume and os.path.exists(chk):
        # Robustez: si prod.chk quedó truncado/corrupto (p.ej. escritura cortada por
        # un apagado), intentar el respaldo prod.chk.bak antes de rendirse.
        try:
            sim.loadCheckpoint(chk)
        except Exception as e:
            bak = chk + ".bak"
            print(f"!! Falló cargar {chk} ({e}); intentando {bak}...", flush=True)
            sim.loadCheckpoint(bak)
        done_ps = sim.context.getState().getTime().value_in_unit(unit.picoseconds)
        # Validar que el tiempo cargado es coherente con el objetivo: detecta un
        # checkpoint basura (t<=0 o t mayor que la producción pedida) antes de correr.
        assert 0 < done_ps <= args.ns * 1000, (
            f"Tiempo del checkpoint fuera de rango: {done_ps:.0f} ps "
            f"(objetivo {args.ns * 1000:.0f} ps)")
        print(f"== Reanudando desde checkpoint (t = {done_ps:.0f} ps) ==", flush=True)
        append = True
    else:
        # arrancar desde el estado equilibrado
        equil = os.path.join(SYS, "equil_Y113F", "equil_free.xml")
        with open(equil) as fh:
            state = mm.XmlSerializer.deserialize(fh.read())
        # Cargar SOLO posiciones y caja, NO el estado completo: el estado equilibrado
        # guarda el parámetro global 'k' del restraint de equilibración, que el System
        # de producción (sin restraints) no tiene -> setState() falla con ese 'k'.
        sim.context.setPositions(state.getPositions())
        sim.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
        sim.context.setVelocitiesToTemperature(TEMP)   # re-sampleo de velocidades para HMR
        # Reiniciar el reloj a 0: el estado equilibrado trae el tiempo acumulado de
        # la equilibración (~13 ns). Sin esto, el conteo de pasos y --resume contarían
        # mal y la producción podría completar <args.ns reales (bug que marcó Codex).
        sim.context.setTime(0.0 * unit.picoseconds)
        print(f"== Producción nueva desde el equilibrado ({args.ns} ns) ==", flush=True)

    save_steps = int(SAVE_PS * unit.picoseconds / TIMESTEP)
    chk_steps = int(CHK_PS * unit.picoseconds / TIMESTEP)
    total_steps = int(args.ns * 1000 * unit.picoseconds / TIMESTEP)

    sim.reporters.append(app.DCDReporter(os.path.join(args.out, "prod.dcd"), save_steps, append=append))
    sim.reporters.append(app.StateDataReporter(
        os.path.join(args.out, "prod.csv"), save_steps, step=True, time=True,
        potentialEnergy=True, kineticEnergy=True, temperature=True, density=True,
        volume=True, speed=True, append=append))
    sim.reporters.append(app.CheckpointReporter(chk, chk_steps))
    # progreso a consola
    sim.reporters.append(app.StateDataReporter(
        sys.stdout, chk_steps, step=True, time=True, temperature=True, speed=True,
        remainingTime=True, totalSteps=total_steps))

    # cuántos pasos faltan
    cur = int(sim.context.getState().getTime() / TIMESTEP) if append else 0
    remaining = max(0, total_steps - cur)
    print(f"== Corriendo {remaining} pasos ({remaining*TIMESTEP.value_in_unit(unit.picoseconds)/1000:.1f} ns) en {args.platform} ==", flush=True)
    sim.step(remaining)

    sim.saveState(os.path.join(args.out, "prod_final.xml"))
    sim.saveCheckpoint(chk)
    # sentinel de "producción completa" (lo usa la cola para no relanzar)
    done_ps = sim.context.getState().getTime().value_in_unit(unit.picoseconds)
    with open(os.path.join(args.out, "prod_DONE.txt"), "w") as f:
        f.write(f"Produccion completa: {done_ps/1000:.1f} ns (objetivo {args.ns} ns)\n")
    print(f"== Producción COMPLETA ({done_ps/1000:.1f} ns) -> {args.out}/prod.dcd ==", flush=True)


if __name__ == "__main__":
    main()
