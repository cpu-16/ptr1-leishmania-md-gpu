# Resumen de resultados

> Métricas principales con sus caveats. El relato completo está en el README.

## Chignolin — calibración

| Métrica | Valor | Nota |
|---|---|---|
| RMSD mínimo (MD vs experimental) | **0.15 Å** | mediana ~0.91 Å |
| AlphaFold2 (ColabFold) vs experimental | **0.90 Å** | pLDDT ~94 |
| MD total | 450 ns @ 340 K | OpenMM · AMBER ff14SB · TIP3P |
| Estado plegado a 340 K | ~99.5 % | sobre-estabilizada |

⚠️ A 340 K no se observó suficiente desplegado/replegado: el **MFPT del MSM no representa la
cinética real** de plegamiento (eso exige escalas de microsegundos). Es el límite que el
proyecto cuantifica, no un resultado fallido.

## PTR1 de _Leishmania panamensis_ — aplicación

| Métrica | Valor | Nota |
|---|---|---|
| Sistema | tetrámero holo, **114 856 átomos** | + NADPH + sustrato |
| Identidad de secuencia vs _L. major_ | 73.9 % | plantilla cristalográfica 1E92 |
| MD de producción | **100 ns** | RMSD Cα ≈ 2.0 Å (estable) |
| Tríada catalítica conservada | **83–100 %** | en los 4 protómeros |
| Distancia de hidruro | 3.8–4.0 Å | geometría compatible ≠ actividad |
| MD vs cristal 1E92 | 1.05 Å | convergencia, **no** validación ciega |
| AF3 vs MD / AF3 vs 1E92 | 0.73 / 0.28 Å | AF3 usó cristales como plantilla |
| Mutantes _in silico_ (Y113F, A112S+Y113F) | 30 ns c/u | Tyr114 = contacto accesorio, preliminar |

⚠️ Una réplica por condición. «Geometría compatible» ≠ actividad catalítica. No se realizó
cribado virtual (trabajo futuro). La lectura de la Tyr114 es preliminar (una réplica).
