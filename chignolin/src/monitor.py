#!/usr/bin/env python
"""Dashboard de terminal — monitorea el muestreo en vivo.

Junta en una sola pantalla (refrescada cada pocos segundos):
  - trayectorias hechas / corriendo / planificadas
  - nanosegundos acumulados y barra de progreso
  - velocidad actual (ns/día) leída de los .log de OpenMM
  - tiempo estimado restante (ETA)
  - estado de la GPU (si nvidia-smi está disponible)

Uso (en OTRA terminal, mientras corre adaptive_sampling.py):
    python src/monitor.py
    python src/monitor.py --once              # imprime una vez y sale
    python src/monitor.py --config config_test.yaml
"""
import argparse
import glob
import os
import shutil
import subprocess
import time

from common import load_config


def gpu_status():
    """Devuelve la línea de estado de la GPU vía nvidia-smi, o None."""
    if shutil.which("nvidia-smi") is None:
        return None
    try:
        out = subprocess.check_output(
            ["nvidia-smi",
             "--query-gpu=utilization.gpu,memory.used,memory.total,temperature.gpu,power.draw",
             "--format=csv,noheader,nounits"], text=True).strip()
        u, mu, mt, temp, pw = [x.strip() for x in out.split(",")]
        return f"uso {u}% | mem {mu}/{mt} MB | {temp}°C | {pw} W"
    except Exception:
        return None


def last_speed(logfile):
    """Última velocidad (ns/día) registrada en un .log de StateDataReporter."""
    try:
        with open(logfile) as f:
            lines = f.readlines()
        if len(lines) < 2:
            return None
        header = lines[0].lstrip("#").strip().split(",")
        idx = next((i for i, h in enumerate(header) if "Speed" in h), None)
        if idx is None:
            return None
        return float(lines[-1].strip().split(",")[idx])
    except Exception:
        return None


def partial_ns(logfile, cfg):
    """ns ya simulados en una trayectoria EN CURSO (cuenta líneas del .log)."""
    try:
        with open(logfile) as f:
            n = max(0, sum(1 for _ in f) - 1)   # menos la cabecera
        return min(n * cfg["save_interval_ps"] / 1000.0, cfg["traj_length_ns"])
    except Exception:
        return 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--once", action="store_true")
    ap.add_argument("--interval", type=float, default=5.0)
    args = ap.parse_args()

    cfg = load_config(args.config)
    tdir = cfg["traj_dir"]
    planned = cfg["trajectories_per_round"] * cfg["rounds"]
    ns_each = cfg["traj_length_ns"]
    ns_plan = planned * ns_each

    while True:
        done = glob.glob(os.path.join(tdir, "*_final.xml"))
        logs = glob.glob(os.path.join(tdir, "*.log"))
        running = [l for l in logs if not os.path.exists(l[:-4] + "_final.xml")]
        n_done = len(done)
        partial = [partial_ns(l, cfg) for l in running]   # progreso de las en curso
        ns_done = n_done * ns_each + sum(partial)
        speeds = [s for s in (last_speed(l) for l in running) if s]
        speed = sum(speeds)

        frac = (ns_done / ns_plan) if ns_plan else 0.0
        bar = "#" * int(frac * 30) + "-" * (30 - int(frac * 30))

        os.system("clear")
        print("=" * 54)
        print("  MONITOR — muestreo (mini Folding@Home)")
        print("=" * 54)
        print(f"Config: {args.config}")
        cur = f" | en curso al {max(partial)/ns_each*100:.0f}%" if (running and ns_each) else ""
        print(f"Trayectorias: {n_done}/{planned} completas{cur}")
        print(f"[{bar}] {frac*100:5.1f}%")
        print(f"ns acumulados: {ns_done:.0f} / {ns_plan:.0f} ns")
        if speed:
            print(f"velocidad actual: {speed:.0f} ns/día ({len(speeds)} traj. activas)")
            if ns_done < ns_plan:
                eta_h = (ns_plan - ns_done) / speed * 24.0
                print(f"ETA aprox.: {eta_h:.1f} h a esta velocidad")
        else:
            print("velocidad: (aún sin datos en los .log)")
        g = gpu_status()
        print(f"GPU: {g}" if g else "GPU: nvidia-smi no disponible")
        print(f"\n(refresca cada {args.interval:.0f}s — Ctrl+C para salir)")

        if args.once or n_done >= planned:
            break
        time.sleep(args.interval)


if __name__ == "__main__":
    main()
