#!/usr/bin/env python
"""PASO 2 — Una "unidad de trabajo": corre UNA trayectoria corta.

Esto es exactamente lo que correría CADA computadora del equipo (igual que un
"work unit" de Folding@Home). Arranca desde un estado dado, con una semilla
aleatoria propia, y graba una trayectoria. Es independiente de las demás.

Uso (normalmente lo llama adaptive_sampling.py, pero puedes correrlo a mano):
    python src/run_md.py --start results/prep/equilibrated.xml \
                         --seed 1 --out results/trajectories/test

Salidas:
    <out>.dcd         -> la trayectoria
    <out>.log         -> energías/temperatura
    <out>_final.xml   -> estado final (sirve como semilla de otra ronda)
"""
import argparse
import os

from openmm import app, unit

from common import (load_config, load_system, make_integrator,
                    make_simulation, steps_for_ns, frames_interval)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--start", required=True, help="estado inicial (.xml)")
    ap.add_argument("--seed", type=int, required=True, help="semilla aleatoria / id")
    ap.add_argument("--out", required=True, help="prefijo de salida")
    ap.add_argument("--ns", type=float, default=None, help="longitud en ns (def: config)")
    args = ap.parse_args()

    cfg = load_config(args.config)
    ns = args.ns if args.ns is not None else cfg["traj_length_ns"]

    # cargar el sistema y la topología preparados en el paso 1
    system = load_system(os.path.join(cfg["prep_dir"], "system.xml"))
    pdb = app.PDBFile(os.path.join(cfg["prep_dir"], "system.pdb"))

    integrator = make_integrator(cfg)
    integrator.setRandomNumberSeed(args.seed)
    sim = make_simulation(pdb.topology, system, integrator, cfg)

    # arrancar desde el estado indicado y darle velocidades aleatorias propias
    sim.loadState(args.start)
    sim.context.setVelocitiesToTemperature(
        cfg["temperature_K"] * unit.kelvin, args.seed)

    steps = steps_for_ns(ns, cfg)
    every = frames_interval(cfg)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    sim.reporters.append(app.DCDReporter(args.out + ".dcd", every))
    sim.reporters.append(app.StateDataReporter(
        args.out + ".log", every, step=True, time=True,
        temperature=True, potentialEnergy=True, speed=True))

    print(f"Corriendo {ns} ns ({steps} pasos), semilla {args.seed} -> {args.out}.dcd")
    sim.step(steps)
    sim.saveState(args.out + "_final.xml")
    print("Trayectoria terminada:", args.out)


if __name__ == "__main__":
    main()
