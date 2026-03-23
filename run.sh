#!/usr/bin/env bash
# ============================================================================
# CLIF-LungCx-Epi | Run all analysis scripts
#
# Usage:
#   chmod +x run.sh
#   ./run.sh
#
# Logs are saved to the output/run_<SITE>_<DATE>/ folder alongside results.
# No patient-level PHI is printed — only aggregate counts, file paths, and
# config info.
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Read site name from config.json
CONFIG_FILE="config/config.json"
if [ ! -f "$CONFIG_FILE" ]; then
  echo "ERROR: $CONFIG_FILE not found. Copy config_template.json and fill in your site details."
  exit 1
fi

SITE_NAME=$(python3 -c "import json; print(json.load(open('$CONFIG_FILE'))['site_name'])" 2>/dev/null \
  || grep -o '"site_name"[[:space:]]*:[[:space:]]*"[^"]*"' "$CONFIG_FILE" | sed 's/.*"site_name"[[:space:]]*:[[:space:]]*"//;s/"$//')

if [ -z "$SITE_NAME" ]; then
  echo "ERROR: Could not read site_name from $CONFIG_FILE"
  exit 1
fi

DATE_STAMP=$(date +"%Y%m%d")
OUT_DIR="output/run_${SITE_NAME}_${DATE_STAMP}"
mkdir -p "$OUT_DIR"

LOG_FILE="${OUT_DIR}/pipeline_${SITE_NAME}_${DATE_STAMP}.log"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

# Spinner for long-running steps
spin() {
  local pid=$1
  local msg=$2
  local chars='⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏'
  local i=0
  local start_time=$SECONDS
  while kill -0 "$pid" 2>/dev/null; do
    local elapsed=$(( SECONDS - start_time ))
    local mins=$(( elapsed / 60 ))
    local secs=$(( elapsed % 60 ))
    printf "\r  %s %s [%02d:%02d]" "${chars:i%${#chars}:1}" "$msg" "$mins" "$secs"
    i=$(( i + 1 ))
    sleep 0.2
  done
  wait "$pid"
  local exit_code=$?
  local elapsed=$(( SECONDS - start_time ))
  local mins=$(( elapsed / 60 ))
  local secs=$(( elapsed % 60 ))
  if [ $exit_code -eq 0 ]; then
    printf "\r  ✓ %s [%02d:%02d]\n" "$msg" "$mins" "$secs"
  else
    printf "\r  ✗ %s FAILED [%02d:%02d]\n" "$msg" "$mins" "$secs"
  fi
  return $exit_code
}

echo ""
echo "  CLIF-LungCx-Epi Pipeline"
echo "  Site: $SITE_NAME"
echo "  Log:  $LOG_FILE"
echo "  ─────────────────────────────────"
echo ""

log "=== CLIF-LungCx-Epi pipeline started ==="
log "Site: $SITE_NAME"
log "Working directory: $SCRIPT_DIR"
log "Output directory: $OUT_DIR"
log "Log file: $LOG_FILE"

# Step 0: Install packages
log "--- Step 0: Checking packages ---"
Rscript code/00_renv_restore.R >> "$LOG_FILE" 2>&1 &
spin $! "Checking packages"
log "--- Step 0: Complete ---"

# Run full pipeline (02 sources 01 internally so all objects stay in one R process)
log "--- Running full pipeline (01 → 02) ---"
Rscript code/02_federated_clusters.R >> "$LOG_FILE" 2>&1 &
spin $! "Running analysis pipeline"
log "--- Pipeline complete ---"

echo ""
log "=== Pipeline finished successfully ==="
echo "  Output: $OUT_DIR/"
echo "  Log:    $LOG_FILE"
echo ""
