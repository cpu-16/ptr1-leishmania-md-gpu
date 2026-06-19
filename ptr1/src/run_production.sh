#!/usr/bin/env bash
# ============================================================================
#  Runner de producción robusto (el go/no-go del piloto ya pasó). Corre la
#  producción de 100 ns con reintentos vía --resume y logging persistente.
#  Lanzar desacoplado:
#    setsid bash src/run_production.sh >/dev/null 2>&1 < /dev/null &
# ============================================================================
set -uo pipefail
BASE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$BASE"
PROD="results/system/prod"
mkdir -p "$PROD"
exec > >(tee -a "$PROD/production.log") 2>&1
trap 'echo "[$(date "+%F %T")] [prod] FIN (exit $?)"' EXIT
stamp() { date "+%F %T"; }

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate jic-folding

NS=${NS:-100}
echo "[$(stamp)] [prod] arrancando producción de ${NS} ns (PID $$)..."

for intento in 1 2 3 4 5; do
    if [ -f "$PROD/prod_DONE.txt" ]; then
        echo "[$(stamp)] [prod] ya completa (sentinel)."; break
    fi
    echo "[$(stamp)] [prod] intento $intento/5..."
    # Si ya existe checkpoint, SIEMPRE continuar (--resume) — así es seguro pausar y
    # relanzar este mismo script sin perder lo avanzado. Solo arranca de cero si no hay .chk.
    RESUME_FLAG=""
    [ -f "$PROD/prod.chk" ] && RESUME_FLAG="--resume"
    [ "$intento" -gt 1 ] && sleep 30
    python src/produce_md.py --ns "$NS" $RESUME_FLAG >> "$PROD/run.log" 2>&1
    if [ -f "$PROD/prod_DONE.txt" ]; then
        echo "[$(stamp)] [prod] producción COMPLETA."; break
    fi
    echo "[$(stamp)] [prod] intento $intento no completó; reintentará con --resume."
done

if [ ! -f "$PROD/prod_DONE.txt" ]; then
    echo "[$(stamp)] [prod] NO completó tras 5 intentos. Revisar run.log."
    touch "$PROD/PRODUCCION_FALLIDA.txt"; exit 1
fi
echo "[$(stamp)] [prod] LISTO -> $PROD/prod.dcd"
