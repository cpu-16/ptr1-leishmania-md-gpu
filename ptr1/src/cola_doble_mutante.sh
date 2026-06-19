#!/usr/bin/env bash
# Encolador: espera a que termine la MD del mutante Y113F y entonces lanza el
# doble mutante A112S+Y113F (control de epistasis sugerido por Codex). Una sola
# GPU -> corridas secuenciales. Pensado para correr desacoplado (setsid).
cd "$(dirname "$0")/.."
PY=python
LOG=results/system/cola_doble.log
DONE_Y113F=results/system/prod_Y113F/prod_DONE.txt
DONE_DOBLE=results/system/prod_A112S_Y113F/prod_DONE.txt

echo "[$(date '+%F %T')] encolador iniciado; esperando que termine Y113F..." >> "$LOG"
fails=0
while [ ! -f "$DONE_Y113F" ]; do
  if ! pgrep -f "equilibrate_Y113F|produce_Y113F" >/dev/null 2>&1; then
    fails=$((fails+1))                       # 2 fallos seguidos (10 min) sin proceso ni DONE -> abortar
    if [ "$fails" -ge 2 ] && [ ! -f "$DONE_Y113F" ]; then
      echo "[$(date '+%F %T')] Y113F ya no corre y no hay prod_DONE.txt -> abortando cola (revisar run_Y113F.log)" >> "$LOG"
      exit 1
    fi
  else
    fails=0
  fi
  sleep 300
done

echo "[$(date '+%F %T')] Y113F COMPLETO. Lanzando doble mutante A112S+Y113F..." >> "$LOG"
$PY src/equilibrate_A112S_Y113F.py >> results/system/run_A112S_Y113F.log 2>&1 \
  && $PY src/produce_A112S_Y113F.py --ns 30 >> results/system/run_A112S_Y113F.log 2>&1
rc=$?
echo "[$(date '+%F %T')] doble mutante terminado (exit $rc)." >> "$LOG"
[ -f "$DONE_DOBLE" ] && echo "[$(date '+%F %T')] prod_DONE.txt del doble presente -> OK." >> "$LOG"
