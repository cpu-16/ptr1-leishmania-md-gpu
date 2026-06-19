# -*- coding: utf-8 -*-
"""
Análisis comparativo WT vs mutantes (Y113F, A112S+Y113F) del anclaje del sustrato
por la Tyr114 en PTR1 de L. panamensis. Implementa las exigencias de rigor de Codex:

 - WT limitado a los PRIMEROS 30 ns (comparación justa con los mutantes de 30 ns).
 - Métricas POR CADENA A-D (4 observaciones semi-independientes), no solo el promedio.
 - H-bond ancla·O10(sustrato) con distancia <3.5 Å Y ángulo >120° (ambas direcciones).
 - Block averaging (bloques 1/2/5/10 ns) para estimar incertidumbre con 1 réplica.
 - ¿Se desestabiliza el sustrato sin el anclaje? -> distancia catalítica C4N–C6 y
   distancia ancla·O10, por sitio.
 - Correlación temporal: geometría catalítica cuando el H-bond está activo vs inactivo.
 - Competencia por agua: agua dentro de 3.5 Å del O10 (¿reemplaza el contacto?).

Todas las distancias/ángulos con mdtraj periodic=True sobre coordenadas CRUDAS
(la caja es octaedro truncado; cpptraj y/o superponer-antes dan artefactos PBC).
El residuo ancla se identifica por proximidad del anillo aromático (CZ) al O10 en el
primer frame -> robusto a la numeración y a la mutación Tyr->Phe (en el mutante el
ancla será PHE, sin -OH, y el H-bond es 0 por construcción química).

Uso:  python src/analyze_mutant_compare.py
Corre sobre las condiciones cuya trayectoria ya exista; salta las que falten.
Salida: results/system/mutant_compare/  (tabla .txt, .json, figuras .png)
"""
import os
import json
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "."
SYS = os.path.join(BASE, "results", "system")
OUT = os.path.join(SYS, "mutant_compare")
os.makedirs(OUT, exist_ok=True)

DT_PS = 50.0        # ps por frame (SAVE_PS en produce_md.py)
HB_DIST = 0.35      # nm  (3.5 Å)  criterio de distancia donor-aceptor
HB_ANGLE = 120.0    # grados        criterio de ángulo D-H...A
WAT_DIST = 0.35     # nm  agua O ... O10
# Posición de la sustitución F114Y por cadena A-D (numeración mdtraj = prmtop-1).
# Se identifica el ancla por ESTA posición fija, NO por proximidad (hay otras
# aromáticas cerca del O10 que confunden). Debe ser TYR (WT) o PHE (mutantes).
ANCHOR_RESSEQ = [113, 401, 689, 977]

# (label, prmtop, dcd, max_ns)  -> WT recortado a 30 ns para comparar parejo
CONDITIONS = [
    ("WT",          os.path.join(SYS, "system.prmtop"),             os.path.join(SYS, "prod", "prod.dcd"),             30.0),
    ("Y113F",       os.path.join(SYS, "mutant_Y113F.prmtop"),       os.path.join(SYS, "prod_Y113F", "prod.dcd"),       None),
    ("A112S_Y113F", os.path.join(SYS, "mutant_A112S_Y113F.prmtop"), os.path.join(SYS, "prod_A112S_Y113F", "prod.dcd"), None),
]
CHAINS = "ABCD"


def bonded_h(top, o_idx):
    """Índice del H enlazado a un oxígeno (para el criterio de ángulo)."""
    if o_idx is None:
        return None
    for a, b in top.bonds:
        if a.index == o_idx and b.element is not None and b.element.symbol == "H":
            return b.index
        if b.index == o_idx and a.element is not None and a.element.symbol == "H":
            return a.index
    return None


def atom_in(res, name):
    for a in res.atoms:
        if a.name == name:
            return a.index
    return None


def block_average(boolarr, fpns):
    """Media±sd de la ocupancia (%) en bloques de 1/2/5/10 ns."""
    n = len(boolarr)
    out = {}
    for blk in (1, 2, 5, 10):
        bs = int(blk * fpns)
        nb = n // bs if bs >= 1 else 0
        if bs < 1 or nb < 2:
            continue
        bm = [100.0 * boolarr[i * bs:(i + 1) * bs].mean() for i in range(nb)]
        out[f"{blk}ns"] = [round(float(np.mean(bm)), 1), round(float(np.std(bm)), 1)]
    return out


def analyze(label, prmtop, dcd, max_ns):
    t = md.load_dcd(dcd, top=prmtop)
    if max_ns:
        t = t[:int(max_ns * 1000 / DT_PS)]
    top = t.topology
    n = t.n_frames
    fpns = 1000.0 / DT_PS
    hbi = sorted((r for r in top.residues if r.name == "HBI"), key=lambda r: r.resSeq)
    nph = sorted((r for r in top.residues if r.name == "NPH"), key=lambda r: r.resSeq)
    wat_O = top.select("water and name O")
    arom = [(r, atom_in(r, "CZ")) for r in top.residues if r.name in ("TYR", "PHE")]
    arom = [(r, i) for r, i in arom if i is not None]

    sites = {}
    for k in range(len(hbi)):
        site = CHAINS[k]
        o10 = atom_in(hbi[k], "O10")
        c6 = atom_in(hbi[k], "C6")
        c4n = atom_in(nph[k], "C4N")
        # ancla = la posición de la sustitución F114Y en esta cadena (numeración fija),
        # NO por proximidad. Ventana ±1 por si difiere la numeración mdtraj/prmtop.
        anchor_res = None
        for rs in (ANCHOR_RESSEQ[k], ANCHOR_RESSEQ[k] + 1):
            for r in top.residues:
                if r.resSeq == rs and r.name in ("TYR", "PHE"):
                    anchor_res = r
                    break
            if anchor_res:
                break
        if anchor_res is None:
            sites[site] = dict(anchor="NA", occ=None, blocks={}, dcat_mean=None,
                               dcat_sd=None, anchor_o10_mean=None, wat_occ=None,
                               dcat_on=None, dcat_off=None)
            continue
        anchor_cz = atom_in(anchor_res, "CZ")
        oh = atom_in(anchor_res, "OH")  # None si es PHE

        # H-bond ancla·O10 (solo posible si el ancla es TYR con -OH)
        if anchor_res.name == "TYR" and oh is not None and o10 is not None:
            dist = md.compute_distances(t, [[oh, o10]], periodic=True)[:, 0]
            hh, ho10 = bonded_h(top, oh), bonded_h(top, o10)
            ang_ok = np.zeros(n, bool)
            if hh is not None:   # Tyr dona: O(Tyr)-H...O10
                ang_ok |= np.degrees(md.compute_angles(t, [[oh, hh, o10]], periodic=True)[:, 0]) > HB_ANGLE
            if ho10 is not None:  # O10 dona: O10-H...O(Tyr)
                ang_ok |= np.degrees(md.compute_angles(t, [[o10, ho10, oh]], periodic=True)[:, 0]) > HB_ANGLE
            if hh is None and ho10 is None:
                ang_ok = np.ones(n, bool)
            hbond = (dist < HB_DIST) & ang_ok
            anchor_o10 = md.compute_distances(t, [[oh, o10]], periodic=True)[:, 0] * 10
        else:
            hbond = np.zeros(n, bool)  # PHE: químicamente imposible
            anchor_o10 = md.compute_distances(t, [[anchor_cz, o10]], periodic=True)[:, 0] * 10  # CZ·O10 (ring proxy)

        dcat = md.compute_distances(t, [[c4n, c6]], periodic=True)[:, 0] * 10
        neigh = md.compute_neighbors(t, WAT_DIST, [o10], haystack_indices=wat_O)
        wat_occ = 100.0 * np.mean([len(x) > 0 for x in neigh])
        on, off = hbond, ~hbond
        sites[site] = dict(
            anchor=f"{anchor_res.name}{anchor_res.resSeq}",
            occ=round(100.0 * hbond.mean(), 1),
            blocks=block_average(hbond, fpns),
            dcat_mean=round(float(dcat.mean()), 2), dcat_sd=round(float(dcat.std()), 2),
            anchor_o10_mean=round(float(anchor_o10.mean()), 2),
            wat_occ=round(wat_occ, 1),
            dcat_on=round(float(dcat[on].mean()), 2) if on.any() else None,
            dcat_off=round(float(dcat[off].mean()), 2) if off.any() else None,
        )
    return dict(label=label, ns=round(n / fpns, 1), nframes=n, sites=sites)


def main():
    res = []
    for label, prm, dcd, mx in CONDITIONS:
        if not os.path.exists(dcd):
            print(f"[salta] {label}: aún no existe {dcd}")
            continue
        print(f"[analiza] {label} ({dcd}) ...")
        res.append(analyze(label, prm, dcd, mx))

    if not res:
        print("No hay trayectorias todavía."); return
    json.dump(res, open(os.path.join(OUT, "mutant_compare.json"), "w"), indent=2)

    # tabla de texto
    lines = ["COMPARACIÓN WT vs MUTANTES — anclaje Tyr114·sustrato (PTR1 L. panamensis)\n"]
    for r in res:
        lines.append(f"\n=== {r['label']}  ({r['ns']} ns, {r['nframes']} frames) ===")
        lines.append(f"{'sitio':6} {'ancla':9} {'H-bond%':8} {'±(bloq5ns)':11} {'dcat(Å)':9} "
                     f"{'dcat_on':8} {'dcat_off':9} {'ancla·O10':10} {'agua%':6}")
        for s in CHAINS:
            d = r["sites"].get(s)
            if not d:
                continue
            b5 = d["blocks"].get("5ns", [None, None])[1]
            lines.append(f"{s:6} {d['anchor']:9} {d['occ']:<8} {str(b5):<11} "
                         f"{d['dcat_mean']:<9} {str(d['dcat_on']):<8} {str(d['dcat_off']):<9} "
                         f"{d['anchor_o10_mean']:<10} {d['wat_occ']:<6}")
    txt = "\n".join(lines)
    open(os.path.join(OUT, "mutant_compare.txt"), "w").write(txt + "\n")
    print(txt)

    # figura: ocupancia del H-bond por sitio y condición
    labels = [r["label"] for r in res]
    x = np.arange(4); w = 0.8 / max(1, len(res))
    plt.figure(figsize=(8, 4.5))
    for j, r in enumerate(res):
        occ = [r["sites"].get(s, {}).get("occ", 0) for s in CHAINS]
        err = [r["sites"].get(s, {}).get("blocks", {}).get("5ns", [0, 0])[1] for s in CHAINS]
        plt.bar(x + j * w, occ, w, yerr=err, capsize=3, label=r["label"])
    plt.xticks(x + w * (len(res) - 1) / 2, [f"sitio {s}" for s in CHAINS])
    plt.ylabel("Ocupancia H-bond ancla·O10 (%)"); plt.ylim(0, 100)
    plt.title("Anclaje del sustrato por la posición 114: WT (Tyr) vs mutantes (Phe)")
    plt.legend(); plt.tight_layout()
    plt.savefig(os.path.join(OUT, "hbond_ocupancia_por_sitio.png"), dpi=140); plt.close()
    print(f"\n[OK] -> {OUT}/  (mutant_compare.txt, .json, hbond_ocupancia_por_sitio.png)")


if __name__ == "__main__":
    main()
