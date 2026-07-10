# Reproducibilidad

## Entorno

```bash
conda env create -f environment.yml && conda activate jic-folding
```

Cubre **OpenMM (CUDA), mdtraj y deeptime**. El módulo de PTR1 necesita además **AmberTools**
(`tleap`, `antechamber`) para parametrizar NADPH/HBI, y **PyMOL** para los renders.

## Chignolin (calibración)

```bash
cd chignolin
bash scripts/download_structure.sh     # baja 5awl.pdb del PDB
python src/prepare_system.py           # limpia, solvata, equilibra
python src/adaptive_sampling.py        # muchas trayectorias cortas, resumibles
python src/analysis_rmsd.py            # chequeo de RMSD
python src/build_msm.py                # MSM + PCCA+ + MFPT
```

Cada paso lee `config.yaml`; usa `--config config_test.yaml` para un smoke test (~2 ns) en
minutos. El muestreo imita a Folding@home: trabajos cortos e independientes (se detectan por
su `*_final.xml` y se omiten al reanudar).

## PTR1 de _L. panamensis_ (aplicación)

Requiere las estructuras de partida (monómero de AlphaFold DB + cristal **1E92**) y AmberTools.
Los scripts de `ptr1/src/` usan `BASE = "."` (ejecutar desde `ptr1/` o ajustar `BASE`):

```
build_tetramer.py            # tetrámero holo ensamblado por superposición sobre 1E92
parametrize_nadph.sh / parametrize_hbi.sh / build_nadph.py   # cofactores
equilibrate_system.py → produce_md.py                        # 100 ns
analyze_md.py · convergencia_bloques.py · figura_correlacion_hbond_dcat.py
# Controles in silico:
equilibrate_Y113F.py / produce_Y113F.py · *_A112S_Y113F.py · analyze_mutant_compare.py
```

## Datos NO incluidos (por tamaño)

Trayectorias `*.dcd`, checkpoints `*.chk`, topologías/estados grandes, las MSAs y los
`*_full_data_*.json` de AF3, y la carpeta `results/`. El repositorio trae el **método y los
resultados**; los datos crudos se generan localmente.

## Hardware de referencia

1× **NVIDIA RTX 4060**, CUDA por conda-forge, precisión `mixed`. HMR a 4 amu → paso de 4 fs.
