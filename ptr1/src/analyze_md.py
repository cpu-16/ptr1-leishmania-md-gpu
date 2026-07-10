# -*- coding: utf-8 -*-
"""
Análisis de la MD de PTR1 (tetrámero + 4 NADPH + 4 HBI) con cpptraj + matplotlib.

Calcula y grafica:
  - RMSD Cα: global y por cadena (estabilidad del modelo ensamblado).
  - RMSF por residuo (flexibilidad; señala loops AF y extremos).
  - RMSD de cada ligando (NADPH y HBI) tras alinear la proteína → ¿se salen del bolsillo?
  - Distancia catalítica de transferencia de hidruro C4N(NADPH)↔C6(HBI) por sitio.
  - Radio de giro de la proteína (compactación global).

cpptraj hace el cálculo pesado (lee prmtop + dcd); matplotlib grafica los .dat.

Uso:
  python src/analyze_md.py --traj results/system/prod/prod.dcd
  python src/analyze_md.py --traj results/system/equil/equil.dcd --dt 0.5 --out results/system/equil/analysis

--dt = ps por frame guardado (prod=50, equil=20 según DCDReporter). Se usa para el eje en ns.
"""
import os
import sys
import argparse
import subprocess
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "."
SYS = os.path.join(BASE, "results", "system")
PRMTOP = os.path.join(SYS, "system.prmtop")

# Rangos de residuo por cadena (4×288). Ajustar si cambia el ensamblaje.
CHAIN_RANGES = {"A": (1, 288), "B": (289, 576), "C": (577, 864), "D": (865, 1152)}


def residue_labels(prmtop):
    """Lista de (resnum_1based, label) desde el prmtop."""
    txt = open(prmtop).read()
    m = re.search(r"%FLAG RESIDUE_LABEL\b.*?%FORMAT\([^)]*\)\s*(.*?)%FLAG", txt, re.S)
    labels = m.group(1).split()
    return [(i + 1, lab) for i, lab in enumerate(labels)]


def find_ligands(prmtop):
    """Devuelve listas de resnums (1-based) de NPH y HBI, en orden de cadena."""
    labs = residue_labels(prmtop)
    nph = [i for i, l in labs if l == "NPH"]
    hbi = [i for i, l in labs if l == "HBI"]
    return nph, hbi


def count_frames(top, traj):
    """Cuenta frames de la trayectoria vía cpptraj (no hay mdtraj en el env)."""
    inp = f"parm {top}\ntrajin {traj}\nrun\nquit\n"
    r = subprocess.run(["cpptraj"], input=inp, capture_output=True, text=True)
    m = re.search(r"occur on (\d+) frames", r.stdout)
    return int(m.group(1)) if m else None


def load_dat(path, ycol=1):
    """Lee un .dat de cpptraj (col0=frame). Devuelve (frames, valores)."""
    if not os.path.exists(path):
        return None, None
    d = np.loadtxt(path, comments="#")
    if d.ndim == 1:
        d = d.reshape(-1, 1) if d.size else d.reshape(0, 1)
    return d[:, 0], d[:, ycol]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--traj", required=True, help="trayectoria (.dcd/.nc)")
    ap.add_argument("--top", default=PRMTOP)
    ap.add_argument("--dt", type=float, default=50.0, help="ps por frame guardado")
    ap.add_argument("--start-frame", type=int, default=0,
                    help="frame inicial (1-based). Negativo = relativo al final "
                         "(p.ej. -250 = últimos 250 frames, para evaluar solo el piloto libre)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    out = args.out or os.path.join(os.path.dirname(args.traj), "analysis")
    os.makedirs(out, exist_ok=True)

    nph, hbi = find_ligands(args.top)
    print(f"NADPH en residuos {nph}; HBI en residuos {hbi}")
    # emparejamiento de sitio por orden (k-ésimo NPH con k-ésimo HBI). VERIFICAR con
    # las distancias resultantes: si una sale grande, el sitio real es cruzado.
    pairs = list(zip(nph, hbi))

    # --- rango de frames (para evaluar solo el piloto libre si se pide) ---
    trajin = f"trajin {args.traj}"
    if args.start_frame:
        sf = args.start_frame
        if sf < 0:
            total = count_frames(args.top, args.traj)
            sf = max(1, (total + sf + 1)) if total else 1
            print(f"Trayectoria con {total} frames; analizando desde el frame {sf}")
        trajin = f"trajin {args.traj} {sf} last"

    # --- generar script cpptraj ---
    lines = [f"parm {args.top}", trajin, "autoimage"]
    # referencia = primer frame
    lines.append("reference %s 1" % args.traj)
    # RMSD Cα global + por cadena (fit a proteína)
    lines.append("rms rms_all reference :1-1152@CA out %s/rmsd_all.dat mass" % out)
    for ch, (a, b) in CHAIN_RANGES.items():
        lines.append(f"rms rms_{ch} reference :{a}-{b}@CA out {out}/rmsd_{ch}.dat nofit")
    # RMSF por residuo (tras fit a proteína)
    lines.append("rms reference :1-1152@CA")
    lines.append("atomicfluct out %s/rmsf.dat :1-1152@CA byres" % out)
    # radio de giro de la proteína
    lines.append("radgyr radg :1-1152@CA out %s/radgyr.dat mass" % out)
    # RMSD de cada ligando tras alinear la proteína (sin re-fit -> mide salida del bolsillo)
    lines.append("rms reference :1-1152@CA")
    for k, r in enumerate(nph):
        lines.append(f"rms lignph_{k} reference :{r}&!@H= out {out}/rmsd_nph_{k}.dat nofit")
    for k, r in enumerate(hbi):
        lines.append(f"rms lighbi_{k} reference :{r}&!@H= out {out}/rmsd_hbi_{k}.dat nofit")
    # distancia catalítica C4N(NADPH)↔C6(HBI) por sitio
    for k, (rn, rh) in enumerate(pairs):
        lines.append(f"distance dcat_{k} :{rn}@C4N :{rh}@C6 out {out}/dist_cat_{k}.dat")
    lines.append("run")
    lines.append("quit")
    cpptraj_in = os.path.join(out, "analyze.cpptraj")
    open(cpptraj_in, "w").write("\n".join(lines) + "\n")

    # --- correr cpptraj ---
    print("== Corriendo cpptraj ==")
    r = subprocess.run(["cpptraj", "-i", cpptraj_in], capture_output=True, text=True)
    if r.returncode != 0:
        print("cpptraj FALLÓ:\n", r.stderr[-2000:]); sys.exit(1)

    def to_ns(frames):
        return frames * args.dt / 1000.0

    # --- 1) RMSD por cadena ---
    plt.figure(figsize=(8, 5))
    for ch in CHAIN_RANGES:
        f, v = load_dat(f"{out}/rmsd_{ch}.dat")
        if f is not None:
            plt.plot(to_ns(f), v, label=f"cadena {ch}", lw=1)
    fg, vg = load_dat(f"{out}/rmsd_all.dat")
    if fg is not None:
        plt.plot(to_ns(fg), vg, "k--", label="global", lw=1.5)
    plt.xlabel("tiempo (ns)"); plt.ylabel("RMSD Cα (Å)")
    plt.title("Estabilidad estructural por cadena (PTR1 L. panamensis)")
    plt.legend(); plt.tight_layout(); plt.savefig(f"{out}/rmsd_por_cadena.png", dpi=130); plt.close()

    # --- 2) RMSF por residuo (cadena A como representativa) ---
    rf = load_dat(f"{out}/rmsf.dat")
    if rf[0] is not None:
        resid, fluct = rf
        plt.figure(figsize=(9, 4))
        a, b = CHAIN_RANGES["A"]
        mask = (resid >= a) & (resid <= b)
        plt.plot(resid[mask] - a + 1, fluct[mask], lw=0.9)
        plt.xlabel("residuo (cadena A)"); plt.ylabel("RMSF (Å)")
        plt.title("Flexibilidad por residuo — picos = loops / extremos")
        plt.tight_layout(); plt.savefig(f"{out}/rmsf.png", dpi=130); plt.close()

    # --- 3) RMSD de ligandos ---
    plt.figure(figsize=(8, 5))
    for k in range(len(nph)):
        f, v = load_dat(f"{out}/rmsd_nph_{k}.dat")
        if f is not None:
            plt.plot(to_ns(f), v, label=f"NADPH {chr(65+k)}", lw=1)
    for k in range(len(hbi)):
        f, v = load_dat(f"{out}/rmsd_hbi_{k}.dat")
        if f is not None:
            plt.plot(to_ns(f), v, "--", label=f"HBI {chr(65+k)}", lw=1)
    plt.xlabel("tiempo (ns)"); plt.ylabel("RMSD ligando (Å, proteína alineada)")
    plt.title("Permanencia de cofactor/sustrato en el bolsillo")
    plt.legend(fontsize=8, ncol=2); plt.tight_layout()
    plt.savefig(f"{out}/rmsd_ligandos.png", dpi=130); plt.close()

    # --- 4) Distancia catalítica ---
    plt.figure(figsize=(8, 5))
    cat_final = {}
    for k in range(len(pairs)):
        f, v = load_dat(f"{out}/dist_cat_{k}.dat")
        if f is not None:
            plt.plot(to_ns(f), v, label=f"sitio {chr(65+k)}", lw=1)
            cat_final[chr(65 + k)] = (float(np.mean(v)), float(np.std(v)))
    plt.axhspan(3.0, 4.0, color="green", alpha=0.1, label="rango competente (~3-4 Å)")
    plt.xlabel("tiempo (ns)"); plt.ylabel("d(C4N_NADPH – C6_HBI) (Å)")
    plt.title("Distancia de transferencia de hidruro por sitio activo")
    plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(f"{out}/dist_catalitica.png", dpi=130); plt.close()

    # --- resumen ---
    with open(f"{out}/resumen.txt", "w") as fh:
        fh.write("# Análisis MD — PTR1 L. panamensis (tetrámero holo)\n")
        fh.write(f"Trayectoria: {args.traj}  (dt={args.dt} ps/frame)\n\n")
        if fg is not None:
            fh.write(f"RMSD Cα global: media {np.mean(vg):.2f} Å, final {vg[-1]:.2f} Å\n")
        for ch in CHAIN_RANGES:
            f, v = load_dat(f"{out}/rmsd_{ch}.dat")
            if f is not None:
                fh.write(f"  cadena {ch}: media {np.mean(v):.2f} Å, máx {np.max(v):.2f} Å\n")
        fh.write("\nDistancia catalítica C4N–C6 (media ± sd) por sitio:\n")
        for ch, (m, s) in cat_final.items():
            flag = "  competente" if m < 4.5 else "  ALEJADO (revisar emparejamiento/sitio)"
            fh.write(f"  sitio {ch}: {m:.2f} ± {s:.2f} Å{flag}\n")
        fh.write("\nNADPH/HBI RMSD final (proteína alineada):\n")
        for k in range(len(nph)):
            f, v = load_dat(f"{out}/rmsd_nph_{k}.dat")
            if f is not None:
                fh.write(f"  NADPH {chr(65+k)}: {v[-1]:.2f} Å\n")
        for k in range(len(hbi)):
            f, v = load_dat(f"{out}/rmsd_hbi_{k}.dat")
            if f is not None:
                fh.write(f"  HBI {chr(65+k)}: {v[-1]:.2f} Å\n")

    print(f"== Análisis listo -> {out}/  (rmsd_por_cadena.png, rmsf.png, rmsd_ligandos.png, dist_catalitica.png, resumen.txt) ==")


if __name__ == "__main__":
    main()
