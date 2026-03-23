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

log "=== CLIF-LungCx-Epi pipeline started ==="
log "Site: $SITE_NAME"
log "Working directory: $SCRIPT_DIR"
log "Output directory: $OUT_DIR"
log "Log file: $LOG_FILE"

# Step 0: Install packages (skip if already done)
log "--- Step 0: Checking packages ---"
Rscript code/00_renv_restore.R >> "$LOG_FILE" 2>&1
log "--- Step 0: Complete ---"

# Step 1: Cohort identification + exposome linkage
log "--- Step 1: Cohort identification (01_lungcx_cohort.R) ---"
Rscript code/01_lungcx_cohort.R >> "$LOG_FILE" 2>&1
log "--- Step 1: Complete ---"

# Step 2: Trajectory clustering + association models
log "--- Step 2: Federated clusters (02_federated_clusters.R) ---"
Rscript code/02_federated_clusters.R >> "$LOG_FILE" 2>&1
log "--- Step 2: Complete ---"

log "=== Pipeline finished successfully ==="
log "Review log: $LOG_FILE"
log "Review output: $OUT_DIR/"
