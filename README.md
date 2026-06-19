# De AlphaFold a la dinámica molecular accesible en una GPU de consumo

**Calibración con Chignolin · Aplicación a la PTR1 de _Leishmania panamensis_**

![License](https://img.shields.io/badge/c%C3%B3digo-MIT-blue)
![Figures](https://img.shields.io/badge/figuras-CC--BY--4.0-lightgrey)
![Python](https://img.shields.io/badge/python-3.10+-3776AB)
![MD](https://img.shields.io/badge/MD-OpenMM-E37222)
![MSM](https://img.shields.io/badge/MSM-deeptime-5A4FCF)
![GPU](https://img.shields.io/badge/GPU-NVIDIA%20RTX%204060-76B900)

> **Una GPU de consumo basta para recuperar la _estructura_ de una proteína;
> medir su _dinámica_ completa, todavía no.**

Proyecto de la **Jornada de Iniciación Científica (JIC) 2026**. Estudia, con dinámica
molecular (MD) corrida en una sola **GPU de consumo (NVIDIA RTX 4060)** y software libre,
hasta dónde llega la biología estructural computacional con recursos accesibles:

- **Calibración — Chignolin (CLN025):** una mini-proteína de plegamiento conocido se usa
  para medir el alcance del método contra estructura experimental y contra AlphaFold2.
- **Aplicación — PTR1 de _Leishmania panamensis_:** la pteridina reductasa 1 (PTR1) es un
  blanco para fármacos contra el parásito, sin equivalente directo en humano, cuya
  estructura no se ha determinado experimentalmente. Se modela por homología y se estudia
  la estabilidad y la flexibilidad de su sitio activo.

---

## 🔑 Resultado central

Con una sola GPU de consumo y herramientas de código abierto se puede **recuperar la
estructura** de una proteína con calidad comparable a la experimental, pero **medir su
cinética/dinámica completa** exige un muestreo (microsegundos, múltiples réplicas) que
excede lo que rinde una GPU de consumo. AlphaFold y la MD son **complementarias**: una da
la «foto», la otra intenta la «película».

---

## 📁 Estructura del repositorio

```
.
├── chignolin/              # Calibración del método (pipeline reproducible, config-driven)
│   ├── config.yaml         #   parámetros de la corrida real (450 ns, 340 K)
│   ├── config_test.yaml    #   smoke test (~2 ns) para validar el pipeline en minutos
│   ├── scripts/            #   download_structure.sh (baja 5AWL del PDB)
│   └── src/                #   preparación → muestreo → análisis → MSM (OpenMM + deeptime)
├── ptr1/                   # Aplicación: PTR1 de L. panamensis (tetrámero holo con NADPH)
│   └── src/                #   construcción del modelo, parametrización, MD, mutantes, análisis
├── figures/                # Figuras científicas (PNG, 300 dpi)
│   ├── chignolin/          #   superposición exp/AF/MD, histograma y trayectoria de RMSD, escaneo de T
│   └── ptr1/               #   tetrámero, sitio activo, superposición vs 1E92, red catalítica, convergencia
├── environment.yml         # Entorno conda (OpenMM + CUDA + mdtraj + deeptime)
├── LICENSE                 # MIT (código)
└── LICENSE-figures         # CC-BY-4.0 (figuras)
```

---

## ⚙️ Requisitos

```bash
conda env create -f environment.yml && conda activate jic-folding
```

El entorno cubre **OpenMM (CUDA), mdtraj y deeptime**. El módulo de PTR1 necesita además
**AmberTools** (`tleap`, `antechamber`) para parametrizar los cofactores NADPH/HBI, y
**PyMOL** para los _renders_. Plataforma probada: NVIDIA RTX 4060, CUDA por conda-forge,
precisión `mixed`.

---

## ▶️ Reproducir la calibración (Chignolin)

Cada paso lee `chignolin/config.yaml` (o `--config config_test.yaml` para el smoke test):

```bash
cd chignolin
bash scripts/download_structure.sh     # baja 5awl.pdb desde el PDB
python src/prepare_system.py           # 1) limpia, solvata y equilibra
python src/adaptive_sampling.py        # 2) muchas trayectorias cortas (estilo Folding@home)
python src/analysis_rmsd.py            # 3) chequeo de RMSD
python src/build_msm.py                # 4) MSM + PCCA+ + MFPT
python src/monitor.py                  # (opcional, otra terminal) progreso/ETA/GPU
```

El diseño imita a Folding@home: **muchas trayectorias cortas e independientes**, resumibles,
en vez de una sola larga. Las unidades de trabajo se detectan por su `*_final.xml`.

## ▶️ Reproducir la aplicación (PTR1)

Los scripts de `ptr1/src/` documentan el método completo: construcción del tetrámero holo
por homología (`build_tetramer.py`), parametrización de cofactores
(`parametrize_nadph.sh`, `parametrize_hbi.sh`, `build_nadph.py`), equilibración y
producción (`equilibrate_system.py`, `produce_md.py`), controles _in silico_ de mutantes
(`equilibrate_Y113F.py`, `produce_Y113F.py`, `*_A112S_Y113F.py`) y análisis
(`analyze_md.py`, `analyze_mutant_compare.py`, `convergencia_bloques.py`,
`figura_correlacion_hbond_dcat.py`). Por defecto usan `BASE = "."` (ejecutar desde `ptr1/`
o ajustar `BASE`).

> ⚠️ Reproducir PTR1 de punta a punta requiere las estructuras de partida (monómero de
> AlphaFold DB + cristal **1E92**) y AmberTools; ver «Qué NO incluye este repositorio».

---

## 📊 Resultados (honestos)

**Chignolin — calibración**
- 450 ns de MD a 340 K. **RMSD mínimo 0.15 Å** vs experimental (mediana ~0.91 Å);
  **AlphaFold2 (ColabFold) 0.90 Å**, pLDDT ~94. La estructura se recupera bien.
- A 340 K la proteína quedó **~99.5 % plegada**: no se observó suficiente
  desplegado/replegado para medir la cinética. El MFPT del MSM **no** representa la
  cinética real de plegamiento (que exige escalas de microsegundos), sino una separación
  de subestados del estado plegado.

**PTR1 de _L. panamensis_ — aplicación**
- Modelo por homología: monómero de AlphaFold DB → **tetrámero holo** sobre el cristal
  **1E92** (_L. major_, 73.9 % de identidad) + NADPH + sustrato. **114 856 átomos**.
- **100 ns** estables (RMSD Cα ≈ 2.0 Å); ligandos retenidos; **tríada catalítica
  Ser111–Tyr193–Lys197 conservada 83–100 %** en los cuatro sitios; distancia de hidruro
  **≈ 3.8–4.0 Å** («geometría compatible con la catálisis», _no_ actividad).
- Validación cruzada: MD vs 1E92 **RMSD 1.05 Å**; AF3 vs MD 0.73 Å y AF3 vs 1E92 0.28 Å.
  ⚠️ AF3 usó cuatro cristales de PTR1 como plantilla → es un **chequeo de consistencia /
  triangulación**, no una validación ciega.
- **Diferencia entre especies (preliminar):** la Tyr114 de _L. panamensis_ (sustitución
  F114Y respecto a _L. major_) parece formar un **contacto polar accesorio** con el
  sustrato, ausente en _L. major_. Controles _in silico_ (30 ns) lo señalan como
  **accesorio, no esencial**; la hipótesis de epistasis (“−OH migrado”) **no se sostiene**.
  Es resultado de **una sola réplica**.

---

## ⚖️ Limitaciones (explícitas)

- **Una réplica por condición** — sin estadística entre réplicas; varias observaciones son
  preliminares.
- **MD clásica** — no modela la reacción química: «geometría compatible» ≠ actividad.
- **Modelo por homología** — PTR1 es una aproximación estructural, no una estructura
  experimental nueva.
- **No hubo cribado virtual.** La metodología es base para futuros estudios; el cribado
  quedó **fuera del alcance** de este trabajo.
- Chignolin queda **sobre-estabilizada** a 340 K con este campo de fuerza.

---

## 📦 Qué NO incluye este repositorio

Para mantenerlo liviano y reproducible, **no** se versionan los datos pesados ni los
intermedios: trayectorias (`*.dcd`), checkpoints (`*.chk`), topologías/estados grandes,
ni las carpetas `results/`. El repositorio contiene el **código del método y las figuras**;
reproducir las corridas implica generar/descargar las estructuras de partida y los datos
crudos localmente.

---

## 🧰 Stack científico

OpenMM · AMBER ff14SB / GAFF2 · TIP3P · AlphaFold2 (ColabFold) · AlphaFold3 · deeptime (MSM)
· mdtraj · PyMOL · AmberTools · GPU NVIDIA RTX 4060.

## 📄 Licencia

- **Código:** [MIT](LICENSE).
- **Figuras (`figures/`):** [CC-BY-4.0](LICENSE-figures).

## 🔖 Cómo citar

Si usas este código o estas figuras, cita el artículo asociado del proyecto
(Jornada de Iniciación Científica 2026, en evaluación).
