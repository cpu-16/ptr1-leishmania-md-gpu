# De AlphaFold a la dinámica molecular accesible en una GPU de consumo

**Calibración con Chignolin · Aplicación a la PTR1 de _Leishmania panamensis_**

![License](https://img.shields.io/badge/c%C3%B3digo-MIT-blue)
![Figures](https://img.shields.io/badge/figuras%20y%20media-CC--BY--4.0-lightgrey)
![Python](https://img.shields.io/badge/python-3.10+-3776AB)
![MD](https://img.shields.io/badge/MD-OpenMM-E37222)
![AlphaFold](https://img.shields.io/badge/IA-AlphaFold%202%20%2F%203-0A7E8C)
![GPU](https://img.shields.io/badge/GPU-NVIDIA%20RTX%204060-76B900)

> **Una GPU de consumo basta para recuperar la _estructura_ de una proteína;
> medir su _dinámica_ completa, todavía no.**

Proyecto de la **Jornada de Iniciación Científica (JIC) 2026**. Estudia, con dinámica
molecular (MD) corrida en una sola **GPU de consumo (NVIDIA RTX 4060)** y software libre,
hasta dónde llega la biología estructural computacional con recursos accesibles, y la
contrasta con la predicción estructural por IA (**AlphaFold 2 y 3**):

- **Calibración — Chignolin (CLN025):** una mini-proteína de plegamiento conocido se usa
  para medir el alcance del método contra estructura experimental y contra AlphaFold2.
- **Aplicación — PTR1 de _Leishmania panamensis_:** la pteridina reductasa 1 (PTR1) es una
  enzima clave del parásito, sin equivalente directo en humano, cuya estructura no se ha
  determinado experimentalmente. Se modela por homología y se estudia la estabilidad y la
  flexibilidad de su sitio activo.

---

## 🔑 Resultado central

Con una sola GPU de consumo y herramientas de código abierto se puede **recuperar la
estructura** de una proteína con calidad comparable a la experimental, pero **medir su
cinética/dinámica completa** exige un muestreo (microsegundos, múltiples réplicas) que
excede lo que rinde una GPU de consumo. AlphaFold y la MD son **complementarias**: una da
la «foto», la otra intenta la «película».

---

## 🎬 Animaciones de las simulaciones (`media/`)

| Chignolin — plegada (MD) | PTR1 — tetrámero holo, 100 ns (MD) |
|:---:|:---:|
| ![Chignolin MD](media/chignolin_md.gif) | ![PTR1 MD 100 ns](media/ptr1_md_100ns.gif) |

_También disponibles en `.mp4` (mayor calidad)._

---

## 🧬 Modelos de AlphaFold (`alphafold/`)

El proyecto contrasta la MD con dos generaciones de predicción estructural por IA:

**AlphaFold2 — Chignolin** · `alphafold/chignolin_af2/` (vía [ColabFold](https://github.com/sokrypton/ColabFold))
- 5 modelos `.pdb` + sus `scores`. pLDDT ≈ 94; **RMSD 0.90 Å** frente a la estructura experimental.

**AlphaFold3 — PTR1, tetrámero holo con NADPH** · `alphafold/ptr1_af3/` (vía [AlphaFold Server](https://alphafoldserver.com/))
- 5 modelos `.cif` + `summary_confidences`. **pTM / ipTM = 0.95.**
- AF3 vs nuestro modelo de MD: **0.73 Å** · AF3 vs cristal **1E92**: **0.28 Å**.
- ⚠️ **Caveat honesto:** AF3 usó **cuatro cristales de PTR1 como plantilla** (1E92, 2XOX,
  1P33, 1E7W — incluidos en `ptr1_af3/templates/`). Por eso es un **chequeo de consistencia
  / triangulación**, no una validación a ciegas. Análisis completo en
  [`ptr1_af3/RESULTADO_AF3.md`](alphafold/ptr1_af3/RESULTADO_AF3.md).

Referencias: AlphaFold2 — Jumper et al., _Nature_ 2021 · AlphaFold3 — Abramson et al.,
_Nature_ 2024 · ColabFold — Mirdita et al., _Nat. Methods_ 2022. Los datos de AF3 se
distribuyen bajo los términos de [`ptr1_af3/terms_of_use.md`](alphafold/ptr1_af3/terms_of_use.md).

---

## 📁 Estructura del repositorio

```
.
├── chignolin/          # Calibración del método (pipeline reproducible, config-driven)
│   ├── config.yaml     #   parámetros de la corrida real (450 ns, 340 K)
│   ├── config_test.yaml#   smoke test (~2 ns) para validar el pipeline en minutos
│   ├── scripts/        #   download_structure.sh (baja 5AWL del PDB)
│   └── src/            #   preparación → muestreo → análisis → MSM (OpenMM + deeptime)
├── ptr1/               # Aplicación: PTR1 de L. panamensis (tetrámero holo con NADPH)
│   └── src/            #   construcción del modelo, parametrización, MD, mutantes, análisis
├── alphafold/          # Modelos de IA estructural
│   ├── chignolin_af2/  #   AlphaFold2 (ColabFold): 5 modelos + scores
│   └── ptr1_af3/       #   AlphaFold3: modelos .cif, confianzas, plantillas, figuras, análisis
├── media/              # Animaciones de las trayectorias (mp4 + gif)
├── paper/              # Artículo divulgativo JIC (versión ciega, PDF)
├── figures/            # Figuras científicas (PNG, 300 dpi)
│   ├── chignolin/      #   superposición exp/AF/MD, RMSD, escaneo de temperatura
│   └── ptr1/           #   tetrámero, sitio activo, superposición vs 1E92, red catalítica
├── environment.yml     # Entorno conda (OpenMM + CUDA + mdtraj + deeptime)
├── LICENSE             # MIT (código)
└── LICENSE-figures     # CC-BY-4.0 (figuras y media)
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

El diseño imita a Folding@home: **muchas trayectorias cortas e independientes**, resumibles.

## ▶️ Reproducir la aplicación (PTR1)

Los scripts de `ptr1/src/` documentan el método completo: construcción del tetrámero holo
por homología (`build_tetramer.py`), parametrización de cofactores
(`parametrize_nadph.sh`, `parametrize_hbi.sh`, `build_nadph.py`), equilibración y
producción (`equilibrate_system.py`, `produce_md.py`), controles _in silico_ de mutantes
(`*_Y113F.py`, `*_A112S_Y113F.py`) y análisis (`analyze_md.py`,
`analyze_mutant_compare.py`, `convergencia_bloques.py`, `figura_correlacion_hbond_dcat.py`).
Por defecto usan `BASE = "."` (ejecutar desde `ptr1/` o ajustar `BASE`).

> ⚠️ Reproducir PTR1 de punta a punta requiere las estructuras de partida (monómero de
> AlphaFold DB + cristal **1E92**) y AmberTools; ver «Qué incluye y qué no».

---

## 📊 Resultados (honestos)

**Chignolin — calibración**
- 450 ns de MD a 340 K. **RMSD mínimo 0.15 Å** vs experimental (mediana ~0.91 Å);
  **AlphaFold2 0.90 Å**, pLDDT ~94. La estructura se recupera bien.
- A 340 K la proteína quedó **~99.5 % plegada**: no se observó suficiente
  desplegado/replegado para medir la cinética. El MFPT del MSM **no** representa la
  cinética real de plegamiento (que exige escalas de microsegundos).

**PTR1 de _L. panamensis_ — aplicación**
- Modelo por homología: monómero de AlphaFold DB → **tetrámero holo** sobre el cristal
  **1E92** (_L. major_, 73.9 % de identidad) + NADPH + sustrato. **114 856 átomos**.
- **100 ns** estables (RMSD Cα ≈ 2.0 Å); ligandos retenidos; **tríada catalítica
  Ser111–Tyr193–Lys197 conservada 83–100 %** en los cuatro sitios; distancia de hidruro
  **≈ 3.8–4.0 Å** («geometría compatible con la catálisis», _no_ actividad).
- Validación cruzada: MD vs 1E92 **RMSD 1.05 Å**; AF3 vs MD 0.73 Å y AF3 vs 1E92 0.28 Å
  (chequeo de consistencia; ver sección AlphaFold).
- **Diferencia entre especies (preliminar):** la Tyr114 de _L. panamensis_ (sustitución
  F114Y respecto a _L. major_) parece formar un **contacto polar accesorio** con el
  sustrato. Controles _in silico_ (30 ns) lo señalan como **accesorio, no esencial**; la
  hipótesis de epistasis no se sostiene. Es resultado de **una sola réplica**.

## ⚖️ Limitaciones (explícitas)

- **Una réplica por condición** — sin estadística entre réplicas; varias observaciones son preliminares.
- **MD clásica** — no modela la reacción química: «geometría compatible» ≠ actividad.
- **Modelo por homología** — PTR1 es una aproximación estructural, no una estructura nueva.
- **No hubo cribado virtual.** La metodología es base para futuros estudios; el cribado quedó **fuera del alcance**.
- Chignolin queda **sobre-estabilizada** a 340 K con este campo de fuerza.

## 📦 Qué incluye y qué no

**Incluye:** código del método (Chignolin y PTR1), figuras (300 dpi), animaciones,
modelos de AlphaFold 2/3 con sus métricas de confianza, y el artículo divulgativo (ciego).

**No incluye** (por tamaño): trayectorias de MD (`*.dcd`), checkpoints, topologías/estados
grandes, las MSAs y los `*_full_data_*.json` de AF3, ni la carpeta `results/`. El
repositorio contiene el método y los resultados; reproducir las corridas implica generar
los datos crudos localmente.

---

## 🧰 Stack científico

OpenMM · AMBER ff14SB / GAFF2 · TIP3P · AlphaFold2 (ColabFold) · AlphaFold3 (AlphaFold Server)
· deeptime (MSM) · mdtraj · PyMOL · AmberTools · GPU NVIDIA RTX 4060.

## 📄 Licencia

- **Código:** [MIT](LICENSE).
- **Figuras y animaciones (`figures/`, `media/`):** [CC-BY-4.0](LICENSE-figures).

## 🔖 Cómo citar

Si usas este código, las figuras o los modelos, cita el artículo asociado del proyecto
(Jornada de Iniciación Científica 2026, en evaluación); ver `paper/`.
