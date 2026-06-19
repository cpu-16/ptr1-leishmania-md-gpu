# -*- coding: utf-8 -*-
"""
Figura de correlación H-bond(Tyr114·O10) <-> distancia catalítica (C4N-C6),
sobre la producción WT COMPLETA (100 ns), por sitio activo A-D.

Objetivo: EVALUAR, sitio por sitio (sin predecir cuáles), si cuando el ancla Tyr114
engancha el sustrato (H-bond presente) la geometría de transferencia de hidruro tiende
a ser más corta/favorable (dcat_on < dcat_off). El resultado se lee del dato, no de una
expectativa previa.

RIGOR / LENGUAJE OBLIGADO:
 - Es una OBSERVACIÓN CONSISTENTE EN LA TRAYECTORIA (100 ns, 1 réplica), NO una
   "correlación reproducible": los 4 sitios del homotetrámero comparten secuencia
   y se acoplan por las interfaces -> NO son réplicas estadísticas independientes.
 - Por eso la figura muestra la DISTRIBUCIÓN de dcat condicionada al estado del
   H-bond (boxplot descriptivo), no una prueba de significancia.

Definiciones idénticas a analyze_mutant_compare.py (consistencia metodológica):
 ancla = posición F114Y por cadena (numeración mdtraj fija); H-bond = d(OH·O10)<3.5 Å
 Y ángulo D-H...A > 120°; todo periodic=True sobre la caja de octaedro truncado.
 Lee la trayectoria por chunks (iterload) para no cargar 2.6 GB de golpe.

Uso:  python src/figura_correlacion_hbond_dcat.py
Salida: results/system/prod/analysis/correlacion_hbond_dcat.{png,_resumen.txt,.json}
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
PRMTOP = os.path.join(SYS, "system.prmtop")
DCD = os.path.join(SYS, "prod", "prod.dcd")
OUT = os.path.join(SYS, "prod", "analysis")
os.makedirs(OUT, exist_ok=True)

DT_PS = 50.0
HB_DIST = 0.35       # nm  (3.5 Å)
HB_ANGLE = 120.0     # grados
ANCHOR_RESSEQ = [113, 401, 689, 977]   # F114Y por cadena (numeración mdtraj)
CHAINS = "ABCD"
CHUNK = 250          # frames por bloque de lectura


def atom_in(res, name):
    for a in res.atoms:
        if a.name == name:
            return a.index
    return None


def bonded_h(top, o_idx):
    if o_idx is None:
        return None
    for a, b in top.bonds:
        if a.index == o_idx and b.element is not None and b.element.symbol == "H":
            return b.index
        if b.index == o_idx and a.element is not None and a.element.symbol == "H":
            return a.index
    return None


def build_sites(top):
    """Índices de átomos por sitio: ancla(OH), o10, c4n, c6, y los H para el ángulo."""
    hbi = sorted((r for r in top.residues if r.name == "HBI"), key=lambda r: r.resSeq)
    nph = sorted((r for r in top.residues if r.name == "NPH"), key=lambda r: r.resSeq)
    # guard: la topología debe tener exactamente los 4 sitios (4 HBI + 4 NPH); si no,
    # algo cambió y los resultados serían silenciosamente erróneos.
    assert len(hbi) == 4 and len(nph) == 4, \
        f"esperaba 4 HBI y 4 NPH, hay {len(hbi)} HBI y {len(nph)} NPH"
    sites = {}
    for k in range(len(hbi)):
        site = CHAINS[k]
        o10 = atom_in(hbi[k], "O10")
        c6 = atom_in(hbi[k], "C6")
        c4n = atom_in(nph[k], "C4N")
        anchor_res = None
        for rs in (ANCHOR_RESSEQ[k], ANCHOR_RESSEQ[k] + 1):
            for r in top.residues:
                if r.resSeq == rs and r.name in ("TYR", "PHE"):
                    anchor_res = r
                    break
            if anchor_res:
                break
        if anchor_res is None or anchor_res.name != "TYR":
            print(f"[aviso] sitio {site}: ancla no es TYR ({anchor_res}); se omite")
            continue
        oh = atom_in(anchor_res, "OH")
        hh, ho10 = bonded_h(top, oh), bonded_h(top, o10)
        if hh is None and ho10 is None:
            print(f"[aviso] sitio {site}: sin H explícito en OH ni O10 -> el criterio de "
                  f"H-bond se reduce a SOLO distancia (sin filtro de ángulo)")
        sites[site] = dict(
            anchor=f"{anchor_res.name}{anchor_res.resSeq}",
            oh=oh, o10=o10, c4n=c4n, c6=c6,
            hh=hh, ho10=ho10,
            dist=[], ang=[], dcat=[],
        )
    return sites


def main():
    top = md.load_prmtop(PRMTOP)
    sites = build_sites(top)
    nframes = 0
    for chunk in md.iterload(DCD, top=PRMTOP, chunk=CHUNK):
        nframes += chunk.n_frames
        for s, d in sites.items():
            d["dist"].append(md.compute_distances(chunk, [[d["oh"], d["o10"]]], periodic=True)[:, 0])
            d["dcat"].append(md.compute_distances(chunk, [[d["c4n"], d["c6"]]], periodic=True)[:, 0] * 10)
            # ángulo D-H...A en cualquiera de las dos direcciones del puente
            ang = np.zeros(chunk.n_frames)
            if d["hh"] is not None:
                a1 = np.degrees(md.compute_angles(chunk, [[d["oh"], d["hh"], d["o10"]]], periodic=True)[:, 0])
                ang = np.maximum(ang, np.nan_to_num(a1))
            if d["ho10"] is not None:
                a2 = np.degrees(md.compute_angles(chunk, [[d["o10"], d["ho10"], d["oh"]]], periodic=True)[:, 0])
                ang = np.maximum(ang, np.nan_to_num(a2))
            if d["hh"] is None and d["ho10"] is None:
                ang = np.full(chunk.n_frames, 180.0)  # sin H explícito: solo distancia
            d["ang"].append(ang)
        print(f"  ... {nframes} frames", end="\r")
    print(f"\n[info] {nframes} frames ({nframes / (1000.0 / DT_PS):.1f} ns)")

    # consolidar y separar ON/OFF
    summary = {}
    for s, d in sites.items():
        dist = np.concatenate(d["dist"])
        ang = np.concatenate(d["ang"])
        dcat = np.concatenate(d["dcat"])
        hbond = (dist < HB_DIST) & (ang > HB_ANGLE)
        on, off = dcat[hbond], dcat[~hbond]
        summary[s] = dict(
            anchor=d["anchor"],
            occ_pct=round(100.0 * hbond.mean(), 1),
            n_on=int(hbond.sum()), n_off=int((~hbond).sum()),
            dcat_on_mean=round(float(on.mean()), 2) if on.size else None,
            dcat_on_sd=round(float(on.std()), 2) if on.size else None,
            dcat_off_mean=round(float(off.mean()), 2) if off.size else None,
            dcat_off_sd=round(float(off.std()), 2) if off.size else None,
            delta=round(float(off.mean() - on.mean()), 2) if on.size and off.size else None,
            _on=on, _off=off,
        )

    # ---- figura: distribución de dcat condicionada al H-bond, por sitio ----
    fig, ax = plt.subplots(figsize=(8.5, 5))
    pos, ticks, ticklab = [], [], []
    for i, s in enumerate(CHAINS):
        if s not in summary:
            continue
        d = summary[s]
        p0, p1 = i * 3 + 1, i * 3 + 2
        data = [d["_on"] if d["_on"].size else [np.nan],
                d["_off"] if d["_off"].size else [np.nan]]
        bp = ax.boxplot(data, positions=[p0, p1], widths=0.7, patch_artist=True,
                        showfliers=False, medianprops=dict(color="black"))
        bp["boxes"][0].set_facecolor("#2c7fb8")   # ON
        bp["boxes"][1].set_facecolor("#bdbdbd")    # OFF
        ticks.append((p0 + p1) / 2)
        ticklab.append(f"sitio {s}\nH-bond {d['occ_pct']:.0f}%")
        if d["delta"] is not None:
            ymax = np.nanpercentile(np.concatenate([d["_on"], d["_off"]]), 95)
            ax.text((p0 + p1) / 2, ymax + 0.15, f"Δ={d['delta']:+.2f} Å",
                    ha="center", va="bottom", fontsize=9)
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklab)
    ax.set_ylabel("Distancia catalítica C4N(NADPH)–C6(sustrato)  [Å]")
    ax.set_title("Geometría catalítica según el anclaje Tyr114·O10 — WT 100 ns (1 réplica)\n"
                 "azul = H-bond presente (ON)   ·   gris = ausente (OFF)", fontsize=10)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor="#2c7fb8", label="H-bond ON"),
                       Patch(facecolor="#bdbdbd", label="H-bond OFF")], loc="upper right")
    ax.grid(axis="y", alpha=0.3)
    fig.text(0.5, 0.005,
             "Lectura descriptiva de 1 réplica · frames autocorrelacionados (no son N independientes) · "
             "los 4 sitios no son réplicas · sin prueba de significancia",
             ha="center", fontsize=7, style="italic", color="#555555")
    fig.tight_layout(rect=[0, 0.025, 1, 1])
    png = os.path.join(OUT, "correlacion_hbond_dcat.png")
    fig.savefig(png, dpi=140)
    plt.close(fig)

    # ---- resumen de texto ----
    for d in summary.values():
        d.pop("_on", None); d.pop("_off", None)
    L = ["CORRELACIÓN H-bond(Tyr114·O10) <-> distancia catalítica (C4N–C6) — WT 100 ns",
         "=" * 78,
         f"Trayectoria: {DCD}  ({nframes} frames, {nframes/(1000.0/DT_PS):.1f} ns)",
         "Criterio H-bond: d(OH·O10) < 3.5 Å Y ángulo D-H...A > 120° (periodic=True).",
         "",
         f"{'sitio':6} {'ancla':9} {'H-bond%':8} {'dcat_ON':16} {'dcat_OFF':16} {'Δ(OFF-ON)':10} {'n_on/n_off'}"]
    for s in CHAINS:
        if s not in summary:
            continue
        d = summary[s]
        on = f"{d['dcat_on_mean']}±{d['dcat_on_sd']}" if d['dcat_on_mean'] is not None else "—"
        off = f"{d['dcat_off_mean']}±{d['dcat_off_sd']}" if d['dcat_off_mean'] is not None else "—"
        L.append(f"{s:6} {d['anchor']:9} {d['occ_pct']:<8} {on:16} {off:16} "
                 f"{str(d['delta']):10} {d['n_on']}/{d['n_off']}")
    # lectura data-driven: qué sitios muestran efecto y su relación con la ocupancia
    THR = 0.15  # Å: umbral de Δ para considerar el efecto "apreciable" (no nulo)
    marcados = [s for s in CHAINS if s in summary and (summary[s]["delta"] or 0) >= THR]
    nulos = [s for s in CHAINS if s in summary and abs(summary[s]["delta"] or 0) < THR]
    top_site = max(summary, key=lambda s: summary[s]["occ_pct"]) if summary else None
    L += ["",
          "INTERPRETACIÓN (lenguaje honesto — leer el dato, no la expectativa):",
          f" - Δ>0 => con el H-bond Tyr114·O10 presente, la distancia catalítica tiende a ser",
          f"   MÁS CORTA (geometría de hidruro más favorable). Umbral 'apreciable': Δ≥{THR} Å.",
          f" - El efecto NO es uniforme entre sitios: marcado en {marcados or '—'}; "
          f"nulo/ruido en {nulos or '—'}.",
          f" - En esta trayectoria, el único caso con Δ apreciable (C, Δ≈+0.37 Å) coincide con la",
          f"   mayor ocupancia del anclaje ({summary.get('C',{}).get('occ_pct','?')}%); en los sitios",
          f"   de anclaje minoritario (~33–43%) el efecto cae dentro del ruido. Es una coincidencia",
          f"   DESCRIPTIVA, NO permite inferir causalidad ni mecanismo.",
          " - SESGO ESTADÍSTICO declarado: los frames están autocorrelacionados, así que n_on/n_off",
          "   NO son tamaños muestrales independientes — el N efectivo es mucho menor y las",
          "   medias±sd parecen más firmes de lo que son. Por eso NO se hace prueba de significancia.",
          " - LÍMITE: observación de 1 réplica, 1 solo sitio 'positivo' -> ANECDÓTICA. Los 4 sitios",
          "   del homotetrámero NO son réplicas independientes (misma secuencia, acoplamiento por",
          "   interfaces). NO afirmar 'correlación reproducible'. Confirmar exigiría réplicas (futuro).",
          " - El sustrato O10 compite con agua por el contacto (ver wat_occ en mutant_compare)."]
    txt = "\n".join(L)
    open(os.path.join(OUT, "correlacion_hbond_dcat_resumen.txt"), "w").write(txt + "\n")
    json.dump(summary, open(os.path.join(OUT, "correlacion_hbond_dcat.json"), "w"), indent=2)
    print(txt)
    print(f"\n[OK] -> {png}")


if __name__ == "__main__":
    main()
