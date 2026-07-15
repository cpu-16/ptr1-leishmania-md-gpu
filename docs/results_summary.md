# Resumen de resultados

> Métricas principales con sus caveats. El relato completo está en el README.

## Chignolin: calibración

| Métrica | Valor | Nota |
|---|---|---|
| RMSD mediana (MD vs experimental) | **0.91 Å** | estructura representativa del ensemble: 0.85 Å |
| AlphaFold2 (ColabFold) vs experimental | **0.90 Å** | pLDDT ~94 |
| MD total | 450 ns @ 340 K | OpenMM · AMBER ff14SB · TIP3P |
| Estado plegado a 340 K | ~99.5 % | trayectorias de 15 ns que parten del plegado; no prueba sobre-estabilización |

⚠️ A 340 K no se observó suficiente desplegado/replegado: el **MFPT del MSM no representa la
cinética real** de plegamiento (eso exige escalas de microsegundos). Es el límite que el
proyecto cuantifica, no un resultado fallido.

## PTR1 de _Leishmania panamensis_: aplicación

| Métrica | Valor | Nota |
|---|---|---|
| Sistema | tetrámero holo, **114 856 átomos** | + NADPH + sustrato |
| Identidad de secuencia vs _L. major_ | 73.9 % | plantilla cristalográfica 1E92 |
| MD de producción | **4 réplicas × 100 ns** | RMSD Cα monomérico 1.70–2.12 Å; esqueleto estable y reproducible |
| Red del sitio activo | **75–100 %** ocupancia | contactos de sustrato, cofactor y catalíticas Tyr194/Lys198; el contacto con el sustrato **varía entre réplicas** |
| Contactos con el sustrato | variables entre réplicas | réplica 3: un sitio se aleja a 16.7 Å; geometría compatible ≠ actividad |
| MD vs cristal 1E92 | 1.05 Å | convergencia, **no** validación ciega |
| AF3 vs MD / AF3 vs 1E92 | 0.73 / 0.28 Å | AF3 usó cristales como plantilla |
| Mutantes _in silico_ (Y113F, A112S+Y113F) | 30 ns c/u | Tyr114 = contacto accesorio, preliminar |

⚠️ La variante natural tiene 4 réplicas independientes; los controles mutantes, una réplica. «Geometría compatible» ≠ actividad catalítica. La lectura de
la Tyr114 es preliminar (una réplica).

⚠️ **Corrección (jul 2026).** Un barrido previo a 400 K corría a volumen constante, con una
sobrepresión de unos 200 bar. A presión constante y 400 K durante 100 ns, Chignolin **sí se
despliega** (27.5 % de fotogramas nativos, RMSD hasta 8.5 Å). La fracción nativa a 400 K no es
estable con una sola trayectoria: dos corridas NVT que solo difieren en la semilla y en la
equilibración dan 98 % y 47 % sobre los mismos 10 ns. El límite que este trabajo cuantifica es
el tiempo de simulación acumulable, no el campo de fuerza.
