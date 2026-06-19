#!/usr/bin/env bash
# Watcher de temperatura GPU para la producción MD de PTR1 (Rama B, JIC 2026).
#
# Registra temp/clock/potencia/throttle cada INTERVAL segundos y:
#   - Alerta (log + notificación de escritorio) si la GPU pasa de ALERT °C.
#   - Vigila el proceso de producción: avisa si MUERE o si COMPLETA (prod_DONE.txt).
#
# Pensado para correr desacoplado de la sesión:
#   setsid nohup bash src/watch_temp.sh >/dev/null 2>&1 < /dev/null &
#
# Ver en vivo:  tail -f results/system/prod/temp_watch.log

PROD="./results/system/prod"
LOG="$PROD/temp_watch.log"
DONE="$PROD/prod_DONE.txt"
ALERT=90            # °C: umbral de alerta temprana (shutdown del HW ronda ~95-97 °C)
INTERVAL=60         # segundos entre lecturas
NOTIFY_COOLDOWN=600 # s mínimo entre notificaciones de calor para no spamear
last_notify=0

notify() {  # notificación de escritorio; no falla si no hay sesión gráfica accesible
  DISPLAY=:0 DBUS_SESSION_BUS_ADDRESS="unix:path=/run/user/$(id -u)/bus" \
    notify-send -u critical "MD PTR1 — $1" "$2" 2>/dev/null || true
}

echo "[$(date '+%F %T')] watcher iniciado (alerta >= ${ALERT} C, cada ${INTERVAL}s)" >> "$LOG"

while true; do
  ts=$(date '+%F %T')
  read -r temp clock power < <(nvidia-smi \
    --query-gpu=temperature.gpu,clocks.gr,power.draw \
    --format=csv,noheader,nounits 2>/dev/null | head -1 | tr -d ' ' | tr ',' ' ')
  thr=$(nvidia-smi --query-gpu=clocks_event_reasons.sw_thermal_slowdown \
    --format=csv,noheader 2>/dev/null | head -1 | tr -d ' ')

  if [ -z "$temp" ]; then
    echo "[$ts] !! no se pudo leer la GPU (nvidia-smi sin respuesta)" >> "$LOG"
  else
    flag=""
    if [ "$temp" -ge "$ALERT" ] 2>/dev/null; then
      flag="  <== ALERTA >= ${ALERT}C"
      now=$(date '+%s')
      if [ $(( now - last_notify )) -ge "$NOTIFY_COOLDOWN" ]; then
        notify "GPU caliente: ${temp}C" "Throttle térmico. Revisa ventilación. (clock ${clock} MHz)"
        last_notify=$now
      fi
    fi
    echo "[$ts] GPU ${temp}C | clock ${clock} MHz | ${power} W | thermal_slowdown=${thr}${flag}" >> "$LOG"
  fi

  # ¿Sigue viva la producción?
  if ! pgrep -f "produce_md.py --ns 100 --resume" >/dev/null 2>&1; then
    if [ -f "$DONE" ]; then
      echo "[$ts] ✅ Producción COMPLETA: $(tr -d '\n' < "$DONE")" >> "$LOG"
      notify "Producción COMPLETA" "$(tr -d '\n' < "$DONE")"
    else
      echo "[$ts] ⚠️  El proceso de producción YA NO ESTÁ y no hay prod_DONE.txt (¿se cortó?)." >> "$LOG"
      notify "Producción detenida" "El proceso terminó sin prod_DONE.txt. Revisa run_resume.log."
    fi
    echo "[$ts] watcher finalizado." >> "$LOG"
    exit 0
  fi

  sleep "$INTERVAL"
done
