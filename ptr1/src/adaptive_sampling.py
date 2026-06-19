#!/usr/bin/env python
"""PASO 3 — El "mini Folding@Home": reparte muchas trayectorias cortas.

Este es el coordinador. En cada RONDA lanza varias trayectorias independientes
(llamando a run_md.py) y luego elige desde dónde arrancar la siguiente ronda.
Esa es la idea del muestreo adaptativo: en vez de una simulación larga, muchas
cortas y baratas que en conjunto exploran mejor el plegamiento.

  - En UNA sola máquina con 1 GPU: deja parallel_workers=1 (las corre en fila).
  - En VARIAS máquinas (el "cluster de estudiantes"): cada quien corre run_md.py
    sobre las semillas que le toquen y guarda en una carpeta compartida
    (Syncthing/NFS). Ver el README, sección "Repartir entre computadoras".

Uso:
    python src/adaptive_sampling.py
"""
import argparse
import glob
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

from common import load_config


def run_one(job):
    """Lanza run_md.py como subproceso (= una unidad de trabajo)."""
    config, start, seed, out, py = job
    here = os.path.dirname(os.path.abspath(__file__))
    cmd = [py, os.path.join(here, "run_md.py"),
           "--config", config, "--start", start,
           "--seed", str(seed), "--out", out]
    subprocess.run(cmd, check=True)
    return out


def select_seeds(cfg, round_idx, n_needed):
    """Elige los estados iniciales de la siguiente ronda.

    Ronda 0: todas arrancan del estado equilibrado.
    Rondas siguientes: re-siembra desde los estados finales de la ronda previa.

    >>> PUNTO DE MEJORA (adaptive sampling "de verdad") <<<
    Aquí es donde, más adelante, conviene clusterizar lo ya muestreado y
    priorizar los estados POCO visitados (microestados de baja población o de la
    frontera plegado/desplegado). Eso acelera mucho ver eventos de plegamiento.
    De momento usamos una estrategia simple (round-robin sobre estados finales).
    """
    eq = os.path.join(cfg["prep_dir"], "equilibrated.xml")
    if round_idx == 0:
        return [eq] * n_needed
    finals = sorted(glob.glob(
        os.path.join(cfg["traj_dir"], f"round{round_idx-1}_*_final.xml")))
    if not finals:
        return [eq] * n_needed
    reps = (finals * (n_needed // len(finals) + 1))[:n_needed]
    return reps


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.yaml")
    args = ap.parse_args()
    cfg = load_config(args.config)
    os.makedirs(cfg["traj_dir"], exist_ok=True)

    py = sys.executable
    n_traj = cfg["trajectories_per_round"]
    rounds = cfg["rounds"]
    workers = cfg["parallel_workers"]

    total_ns = n_traj * rounds * cfg["traj_length_ns"]
    print(f"Plan: {rounds} rondas x {n_traj} trayectorias x "
          f"{cfg['traj_length_ns']} ns = {total_ns} ns en total\n")

    seed = 0
    for r in range(rounds):
        seeds_states = select_seeds(cfg, r, n_traj)
        jobs = []
        for i in range(n_traj):
            out = os.path.join(cfg["traj_dir"], f"round{r}_{i:03d}")
            if os.path.exists(out + "_final.xml"):
                print(f"  (ya existe {out}, se omite)")
                seed += 1
                continue
            jobs.append((args.config, seeds_states[i], seed, out, py))
            seed += 1
        print(f"=== Ronda {r}: {len(jobs)} trayectorias, {workers} en paralelo ===")
        if jobs:
            with ThreadPoolExecutor(max_workers=workers) as ex:
                list(ex.map(run_one, jobs))

    print(f"\nMuestreo terminado. Trayectorias en '{cfg['traj_dir']}/'.")
    print("Ahora corre:  python src/analysis_rmsd.py   (chequeo rápido)")
    print("        y/o:  python src/build_msm.py       (cinética / MSM)")


if __name__ == "__main__":
    main()
