# Validación AlphaFold 3 — PTR1 tetrámero holo de *L. panamensis*

**Job:** `PTR1_Lpanamensis_tetramero_NADPH` (alphafoldserver.com, 31 may 2026)
**Entrada:** 4×288 residuos (UniProt A0A088SA10) + 4× NADPH (CCD_NDP). Modelo de mayor confianza: `model_0`.

## Métricas de confianza
- **pTM = 0.95, ipTM = 0.95, ranking_score = 0.97** · `fraction_disordered` 0.03 · `has_clash` 0.0.
- `chain_pair_iptm` proteína–su NADPH = **0.97** (cofactor muy confiado); interfaces proteína–proteína 0.89–0.93 (tetrámero sólido).

## Superposición (monómero, una subunidad)
| Comparación | RMSD pliegue | Centroide NADPH |
|---|---|---|
| AF3 vs nuestro modelo MD (panamensis) | **0.73 Å** | **0.9 Å** |
| AF3 vs cristal 1E92 (major, experimental) | **0.28 Å** | **0.3 Å** |

(El alineamiento del tetrámero completo falla por la simetría de las 4 cadenas idénticas; el monómero es la comparación válida.)

## Interpretación HONESTA (caveat clave)
AF3 **usó como plantillas 4 estructuras cristalográficas de PTR1**: **1E92, 2XOX, 1P33, 1E7W** (confirmado en `extracted/templates/`). Por tanto la predicción **NO es a ciegas** — AF3 conocía el pliegue y el sitio del cofactor de PTR1, incluido el propio 1E92.

**Lo que SÍ vale (consistencia / triangulación):**
- Nuestro modelo homólogo + 100 ns de MD **se mantuvo consistente** con lo que AF3 predice usando las mejores plantillas experimentales (RMSD 0.73 Å). → la MD no derivó a algo artificial.
- Nuestro **NADPH trasplantado** (0.9 Å de AF3, 0.3 Å del cristal) está **bien colocado**. → valida el montaje del cofactor.
- El **tetrámero** se ensambla con confianza (ipTM 0.95).
- Cuatro vías convergen a <1 Å: cristal 1E92 + modelo homólogo + MD + AF3. Argumento de robustez para la JIC.

**Lo que NO vale:**
- NO es validación *independiente* (AF3 usó 1E92 como plantilla; la coincidencia es esperada).
- NO prueba nada nuevo sobre rasgos específicos de *panamensis* (Tyr114, etc.) más allá de lo que plantilla+secuencia ya implican.

**Conclusión:** chequeo de consistencia *aprobado* — el modelo es sólido y de calidad publicable, pero presentarlo como "convergencia con solapamiento de método", no como confirmación a ciegas.

Figura: `af3_vs_cristal_monomero.png`.
