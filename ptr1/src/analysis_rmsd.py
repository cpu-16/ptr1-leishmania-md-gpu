#!/usr/bin/env python
"""Chequeo rápido — ¿la proteína llega a su forma plegada?

Calcula el RMSD de cada frame respecto a la estructura experimental nativa.
Si ves valores que bajan cerca de 0 (y a veces suben), ¡estás viendo
plegamiento y desplegamiento! Es la primera señal de que el muestreo funciona.

Uso:
    python src/analysis_rmsd.py
Salidas (en results/msm/):
    rmsd_hist.png    -> histograma de RMSD (debería tener un pico "plegado" bajo)
    rmsd_traj.png    -> RMSD vs tiempo de una trayectoria de ejemplo
"""
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

from common import load_config


def main():
    cfg = load_config()
    top = os.path.join(cfg["prep_dir"], "system.pdb")
    native = md.load(cfg["native_pdb"])
    os.makedirs(cfg["msm_dir"], exist_ok=True)

    dcds = sorted(glob.glob(os.path.join(cfg["traj_dir"], "*.dcd")))
    if not dcds:
        raise SystemExit("No hay trayectorias. Corre primero adaptive_sampling.py")

    # Comparar SOLO los Cα. La nativa no tiene hidrógenos y la simulación sí,
    # pero ambas tienen 10 Cα; recortamos a Cα en las dos para que los índices
    # coincidan. (md.rmsd ya hace el mejor alineamiento rígido.)
    nat_ca = native.atom_slice(native.topology.select("name CA and protein"))
    all_rmsd = []
    first_traj_rmsd = None
    for d in dcds:
        t = md.load(d, top=top)
        t_ca = t.atom_slice(t.topology.select("name CA and protein"))
        rmsd = md.rmsd(t_ca, nat_ca)
        all_rmsd.append(rmsd)
        if first_traj_rmsd is None:
            first_traj_rmsd = rmsd

    rmsd = np.concatenate(all_rmsd) * 10.0   # nm -> Å para mostrar
    folded = np.mean(rmsd < cfg["folded_rmsd_nm"] * 10.0)
    print(f"Frames totales: {rmsd.size}")
    print(f"RMSD mínimo: {rmsd.min():.2f} Å | mediana: {np.median(rmsd):.2f} Å")
    print(f"Fracción 'plegada' (RMSD < {cfg['folded_rmsd_nm']*10:.1f} Å): {folded:.1%}")

    plt.figure()
    plt.hist(rmsd, bins=60)
    plt.xlabel("RMSD a la estructura nativa (Å)")
    plt.ylabel("frecuencia")
    plt.title("Distribución de RMSD (¿hay pico plegado a la izquierda?)")
    plt.savefig(os.path.join(cfg["msm_dir"], "rmsd_hist.png"), dpi=130)

    plt.figure()
    t_ps = np.arange(first_traj_rmsd.size) * cfg["save_interval_ps"]
    plt.plot(t_ps / 1000.0, first_traj_rmsd * 10.0)
    plt.xlabel("tiempo (ns)")
    plt.ylabel("RMSD (Å)")
    plt.title("RMSD vs tiempo (1ª trayectoria)")
    plt.savefig(os.path.join(cfg["msm_dir"], "rmsd_traj.png"), dpi=130)
    print(f"Gráficas guardadas en {cfg['msm_dir']}/")


if __name__ == "__main__":
    main()
