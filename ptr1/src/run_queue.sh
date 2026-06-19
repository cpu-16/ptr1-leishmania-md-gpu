#!/usr/bin/env bash
# ============================================================================
#  Cola automática robusta: espera el fin de la equilibración+piloto, valida
#  (go/no-go multi-criterio) y lanza la producción sola, con reintentos.
#  Diseñada para "fire-and-forget" toda la noche (revisada con Codex).
#
#  Lanzar desacoplada:
#    setsid bash src/run_queue.sh >/dev/null 2>&1 < /dev/null &
#
#  Todo el registro va a results/system/prod/queue.log (heartbeat cada minuto).
#  Si algo va mal, deja un archivo-bandera en results/system/prod/ explicando qué.
# ============================================================================
set -uo pipefail   # sin -e: manejamos los errores explícitamente (hay greps condicionales)

BASE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$BASE"
PROD="results/system/prod"
EQUIL="results/system/equil"
mkdir -p "$PROD"
LOG="$PROD/queue.log"

# logging persistente + marca de salida (detecta muerte silenciosa)
exec > >(tee -a "$LOG") 2>&1
trap 'echo "[$(date "+%F %T")] [cola] FIN (exit $?)"' EXIT
stamp() { date "+%F %T"; }

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate jic-folding

EQUIL_FREE="$EQUIL/equil_free.xml"
EQUIL_DCD="$EQUIL/equil.dcd"
EQUIL_CSV="$EQUIL/equil.csv"
ANALYSIS="$EQUIL/analysis"
WAIT_TIMEOUT=${WAIT_TIMEOUT:-86400}   # 24 h máximo de espera

bandera() { touch "$PROD/$1"; echo "[$(stamp)] [cola] bandera: $1"; }

echo "[$(stamp)] [cola] iniciada (PID $$). Esperando equil_free.xml..."

# --- 1) esperar fin de equilibración, con timeout y heartbeat ---
elapsed=0
while [ ! -f "$EQUIL_FREE" ]; do
    sleep 60; elapsed=$((elapsed+60))
    echo "[$(stamp)] [cola] esperando equilibración... (${elapsed}s)"
    if [ "$elapsed" -ge "$WAIT_TIMEOUT" ]; then
        echo "[$(stamp)] [cola] TIMEOUT esperando equilibración."
        bandera "EQUIL_TIMEOUT.txt"; exit 1
    fi
done
sleep 10
echo "[$(stamp)] [cola] equilibración + piloto TERMINARON."

# --- 2) chequeo de NaN/inf en el log de la simulación (explosión) ---
if grep -qiE "nan|inf" "$EQUIL_CSV" 2>/dev/null; then
    echo "[$(stamp)] [cola] NaN/inf en equil.csv -> el piloto explotó."
    bandera "PILOTO_NAN.txt"; exit 1
fi

# --- 3) integridad del DCD (frames > 0 via cpptraj) ---
NF=$(printf 'parm results/system/system.prmtop\ntrajin %s\nrun\nquit\n' "$EQUIL_DCD" \
     | cpptraj 2>/dev/null | grep -oE "occur on [0-9]+ frames" | grep -oE "[0-9]+")
if [ -z "${NF:-}" ] || [ "$NF" -lt 50 ]; then
    echo "[$(stamp)] [cola] DCD inválido o muy corto (frames=${NF:-0})."
    bandera "EQUIL_DCD_CORRUPTO.txt"; exit 1
fi
echo "[$(stamp)] [cola] DCD OK ($NF frames)."

# --- 4) analizar SOLO el tramo del piloto libre (últimos ~250 frames = 5 ns) ---
echo "[$(stamp)] [cola] analizando el piloto (últimos 250 frames)..."
if ! python src/analyze_md.py --traj "$EQUIL_DCD" --dt 20 --start-frame -250 --out "$ANALYSIS"; then
    echo "[$(stamp)] [cola] analyze_md.py falló."
    bandera "ANALISIS_FALLIDO.txt"; exit 1
fi

# --- 5) GO/NO-GO multi-criterio: RMSD (estructura) + T y densidad (termodinámica) ---
RESUMEN="$ANALYSIS/resumen.txt"
RMSD=$(grep -i "RMSD Cα global" "$RESUMEN" 2>/dev/null | grep -oE "media [0-9.]+" | grep -oE "[0-9.]+" | head -1)
# T (col 4) y densidad (col 6) promediadas sobre las últimas 100 líneas del piloto
TEMP=$(tail -n 100 "$EQUIL_CSV" 2>/dev/null | awk -F',' 'NF>=6{s+=$4;n++} END{if(n)printf "%.1f",s/n}')
DENS=$(tail -n 100 "$EQUIL_CSV" 2>/dev/null | awk -F',' 'NF>=6{s+=$6;n++} END{if(n)printf "%.4f",s/n}')
echo "[$(stamp)] [cola] piloto: RMSD=${RMSD:-?} Å  T=${TEMP:-?} K  densidad=${DENS:-?} g/mL"

GO=$(python3 -c "
rmsd='${RMSD:-99}'; temp='${TEMP:-999}'; dens='${DENS:-0}'
try:
    ok = float(rmsd)<5.0 and 280<=float(temp)<=320 and 0.95<=float(dens)<=1.10
except: ok=False
print('GO' if ok else 'NOGO')")

if [ "$GO" != "GO" ]; then
    echo "[$(stamp)] [cola] NO-GO: el piloto no pasó los criterios. Producción NO lanzada."
    bandera "PILOTO_REVISAR.txt"; exit 1
fi
echo "[$(stamp)] [cola] GO: el piloto pasó. Lanzando producción de 100 ns."

# --- 6) producción con reintentos (--resume), parando si ya completó ---
for intento in 1 2 3; do
    if [ -f "$PROD/prod_DONE.txt" ]; then
        echo "[$(stamp)] [cola] producción ya completa (sentinel presente)."; break
    fi
    echo "[$(stamp)] [cola] producción intento $intento/3..."
    if [ "$intento" -eq 1 ]; then
        python src/produce_md.py --ns 100 > "$PROD/run.log" 2>&1
    else
        sleep 30
        python src/produce_md.py --ns 100 --resume >> "$PROD/run.log" 2>&1
    fi
    if [ -f "$PROD/prod_DONE.txt" ]; then
        echo "[$(stamp)] [cola] producción COMPLETA."; break
    fi
    echo "[$(stamp)] [cola] intento $intento no completó; reintentando con --resume..."
done

if [ ! -f "$PROD/prod_DONE.txt" ]; then
    echo "[$(stamp)] [cola] producción NO completó tras 3 intentos."
    bandera "PRODUCCION_FALLIDA.txt"; exit 1
fi
echo "[$(stamp)] [cola] TODO LISTO. Producción en $PROD/prod.dcd"
