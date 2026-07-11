#!/usr/bin/env python
"""PASO 4: construir el Markov State Model (MSM) a partir de las trayectorias.

El MSM estima los estados metaestables y sus poblaciones (termodinámica). Cuando
hay transiciones suficientes entre plegado y desplegado, también da los tiempos
medios de transición (cinética); si el muestreo no las cubre, el MSM solo separa
subestados del plegado y esos tiempos NO representan la cinética real (ver la
advertencia automática más abajo).

Flujo: features (distancias Cα) -> TICA -> clustering -> MSM -> PCCA+ (macroestados)
       -> poblaciones + tiempo de plegamiento (MFPT).

NOTA: la API de deeptime puede variar entre versiones. Si algo falla, revisa
https://deeptime-ml.github.io (sección Markov state models). Está escrito para
deeptime >= 0.4.

Uso:
    python src/build_msm.py
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

from common import load_config

from deeptime.decomposition import TICA
from deeptime.clustering import KMeans
from deeptime.markov import TransitionCountEstimator
from deeptime.markov.msm import MaximumLikelihoodMSM


def ca_distance_features(traj):
    """Distancias entre todos los pares de Cα: descriptor simple y robusto."""
    ca = traj.topology.select("name CA")
    pairs = [(ca[i], ca[j]) for i in range(len(ca)) for j in range(i + 1, len(ca))]
    return md.compute_distances(traj, pairs)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.yaml")
    args = ap.parse_args()

    cfg = load_config(args.config)
    top = os.path.join(cfg["prep_dir"], "system.pdb")
    os.makedirs(cfg["msm_dir"], exist_ok=True)

    dcds = sorted(glob.glob(os.path.join(cfg["traj_dir"], "*.dcd")))
    if not dcds:
        raise SystemExit("No hay trayectorias. Corre primero adaptive_sampling.py")

    # --- cargar trayectorias (solo proteína) y calcular features + RMSD nativo ---
    print(f"[1/5] Cargando {len(dcds)} trayectorias y calculando features ...")
    native = md.load(cfg["native_pdb"])
    nat_ca = native.atom_slice(native.topology.select("name CA and protein"))
    prot_sel = md.load(top).topology.select("protein")
    feats, rmsds = [], []
    for d in dcds:
        t = md.load(d, top=top).atom_slice(prot_sel)
        feats.append(ca_distance_features(t))
        t_ca = t.atom_slice(t.topology.select("name CA"))   # solo Cα, alinea con la nativa
        rmsds.append(md.rmsd(t_ca, nat_ca))

    # --- TICA: reduce dimensiones quedándose con los modos lentos ---
    print("[2/5] TICA ...")
    tica = TICA(lagtime=cfg["tica_lag"], dim=cfg["tica_dim"]).fit_fetch(feats)
    ys = [tica.transform(f) for f in feats]

    # --- clustering en microestados ---
    print(f"[3/5] Clustering en {cfg['n_clusters']} microestados ...")
    km = KMeans(n_clusters=cfg["n_clusters"], max_iter=200,
                fixed_seed=1).fit_fetch(np.concatenate(ys))
    dtrajs = [km.transform(y).astype(np.int64) for y in ys]

    # --- MSM al lag elegido ---
    print(f"[4/5] Construyendo MSM (lag = {cfg['msm_lag']} frames) ...")
    counts = TransitionCountEstimator(
        lagtime=cfg["msm_lag"], count_mode="sliding").fit_fetch(dtrajs)
    msm = MaximumLikelihoodMSM().fit_fetch(counts.submodel_largest())

    ps_per_frame = cfg["save_interval_ps"]
    its_ns = msm.timescales()[:5] * ps_per_frame / 1000.0
    print("      Timescales implícitas (ns):",
          np.array2string(its_ns, precision=2))

    # --- PCCA+: agrupar en 2 macroestados (plegado / desplegado) y sacar cinética ---
    print("[5/5] PCCA+ (2 macroestados) + cinética ...")
    pcca = msm.pcca(n_metastable_sets=2)
    micro_assign = pcca.assignments               # macroestado de cada microestado
    pops = msm.stationary_distribution            # población por microestado

    # RMSD medio de cada microestado -> el de menor RMSD es "plegado"
    all_dtrajs = np.concatenate(dtrajs)
    all_rmsd = np.concatenate(rmsds)
    micro_rmsd = np.array([
        all_rmsd[all_dtrajs == k].mean() if np.any(all_dtrajs == k) else np.nan
        for k in range(cfg["n_clusters"])
    ])
    # active[i] es el SÍMBOLO original del microestado en el índice interno i del
    # submodelo conexo del MSM. micro_assign (PCCA+) y pops están indexados por ese
    # índice interno i, no por el símbolo original: hay que traducir símbolo -> i.
    active = list(msm.count_model.state_symbols)
    internal_rmsd = np.array([micro_rmsd[s] for s in active])   # RMSD por índice interno
    # macroestado plegado = el que contiene el microestado (superviviente) de menor RMSD
    folded_internal = int(np.nanargmin(internal_rmsd))
    folded_macro = micro_assign[folded_internal]

    folded_set = [i for i in range(len(active)) if micro_assign[i] == folded_macro]
    unfolded_set = [i for i in range(len(active)) if micro_assign[i] != folded_macro]

    p_folded = pops[folded_set].sum()
    print(f"      Población PLEGADA estimada: {p_folded:.1%}")

    # tiempos medios de primer paso (cinética): desplegado->plegado y vuelta
    mfpt_fold = msm.mfpt(unfolded_set, folded_set) * ps_per_frame / 1000.0
    mfpt_unfold = msm.mfpt(folded_set, unfolded_set) * ps_per_frame / 1000.0
    print(f"      Tiempo medio de PLEGAMIENTO (MFPT): {mfpt_fold:.0f} ns")
    print(f"      Tiempo medio de DESPLEGAMIENTO    : {mfpt_unfold:.0f} ns")
    print("      (Chignolin/CLN025 experimental: plegado ~ cientos de ns)")

    # Aviso de validez: si casi todo está plegado por RMSD, el PCCA+ separó
    # SUB-estados del plegado y el MFPT de arriba NO es la cinética real.
    rmsd_folded_frac = float(np.mean(all_rmsd < cfg["folded_rmsd_nm"]))
    if rmsd_folded_frac > 0.90:
        print(f"\n  *** ADVERTENCIA: {rmsd_folded_frac:.0%} de los frames está plegado "
              "(por RMSD).")
        print("      El MSM está dominado por el estado plegado: probablemente separó")
        print("      SUB-estados del plegado, no plegado/desplegado. Por tanto el MFPT")
        print("      de arriba NO es el tiempo de plegamiento real. Hace falta más")
        print("      muestreo (o mayor T / partir de desplegado) para medir la cinética.")

    # --- guardar resumen y una gráfica de timescales ---
    # Si >90% está plegado, el MFPT no es cinética real (ver advertencia): NO lo
    # persistimos como un número válido, para que nadie reutilice cifras retractadas.
    kinetics_valid = rmsd_folded_frac <= 0.90
    with open(os.path.join(cfg["msm_dir"], "resumen.txt"), "w") as f:
        f.write(f"Trayectorias: {len(dcds)}\n")
        f.write(f"Timescales (ns): {its_ns}\n")
        f.write(f"Poblacion plegada: {p_folded:.4f}\n")
        f.write(f"Fraccion de frames plegados (RMSD): {rmsd_folded_frac:.4f}\n")
        if kinetics_valid:
            f.write(f"MFPT plegamiento (ns): {mfpt_fold:.1f}\n")
            f.write(f"MFPT desplegamiento (ns): {mfpt_unfold:.1f}\n")
        else:
            f.write("MFPT plegamiento (ns): NA (muestreo insuficiente: >90% de los "
                    "frames plegado; el MSM separa subestados del plegado, no la cinetica real)\n")
            f.write("MFPT desplegamiento (ns): NA (idem)\n")

    plt.figure()
    plt.bar(range(1, len(its_ns) + 1), its_ns)
    plt.xlabel("modo")
    plt.ylabel("timescale implícita (ns)")
    plt.title("Procesos lentos del MSM")
    plt.savefig(os.path.join(cfg["msm_dir"], "timescales.png"), dpi=130)
    print(f"\nResultados en {cfg['msm_dir']}/ (resumen.txt, timescales.png)")


if __name__ == "__main__":
    main()
