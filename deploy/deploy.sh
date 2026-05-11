#!/usr/bin/env bash
# Deploy this Shiny app to Posit Connect using a per-app YAML under deploy/apps/.
# Prerequisites: rsconnect R package, `rsconnect add` for your server (API key not stored in repo).
# Connect app settings: set EXPRS_MAIN_CONFIG to the same path as APP_YAML (repo-relative), e.g. deploy/apps/heart.yaml
#
# Usage:
#   export CONNECT_API_KEY=...   # or use rsconnect stored credentials
#   ./deploy/deploy.sh deploy/apps/heart.yaml
#   ./deploy/deploy.sh heart     # shorthand for deploy/apps/heart.yaml
#   ./deploy/deploy.sh heart --dry-run
#   ./deploy/deploy.sh heart --verbose    # or -v: bundle step detail + writeManifest(verbose=TRUE)
#   ./deploy/deploy.sh heart --debug      # implies --verbose; print every bundle path; bash set -x during upload

set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

DRY_RUN=0
VERBOSE=0
DEBUG=0
POS=()
for a in "$@"; do
  case "$a" in
    --dry-run) DRY_RUN=1 ;;
    --verbose|-v) VERBOSE=1 ;;
    --debug) DEBUG=1; VERBOSE=1 ;;
    *) POS+=("$a") ;;
  esac
done

if [[ ${#POS[@]} -lt 1 ]]; then
  echo "Usage: $0 <deploy/apps/foo.yaml|foo> [--dry-run] [--verbose|-v] [--debug]" >&2
  exit 1
fi

APP_ARG="${POS[0]}"

if [[ "$APP_ARG" != *.yaml && "$APP_ARG" != *.yml ]]; then
  APP_YAML="deploy/apps/${APP_ARG}.yaml"
else
  APP_YAML="$APP_ARG"
fi

if [[ ! -f "$APP_YAML" ]]; then
  echo "Deploy app file not found: $APP_YAML" >&2
  exit 1
fi

export EXPRS_DEPLOY_REPO="$REPO_ROOT"

R_EXTRA=()
[[ "$DRY_RUN" -eq 1 ]] && R_EXTRA+=(--dry-run)
[[ "$VERBOSE" -eq 1 ]] && R_EXTRA+=(--verbose)
[[ "$DEBUG" -eq 1 ]] && R_EXTRA+=(--debug)

if [[ "$DRY_RUN" -eq 1 ]]; then
  Rscript deploy/write_manifest.R "$APP_YAML" "${R_EXTRA[@]}"
  exit 0
fi

Rscript deploy/write_manifest.R "$APP_YAML" "${R_EXTRA[@]}"

ENV_FILE="$REPO_ROOT/deploy/.cache/last/rsconnect_deploy.env"
if [[ ! -f "$ENV_FILE" ]]; then
  echo "Missing $ENV_FILE after write_manifest.R" >&2
  exit 1
fi
# shellcheck disable=SC1090
set -a
source "$ENV_FILE"
set +a

if [[ -z "${RSCONNECT_ACCOUNT:-}" || -z "${RSCONNECT_TITLE:-}" || -z "${RSCONNECT_APP_ID:-}" ]]; then
  echo "deploy block in YAML must set deploy.account, deploy.title, deploy.app_id" >&2
  exit 1
fi

# Upload phase: keep -v so rsconnect shows progress (separate from manifest --verbose/--debug).
if [[ "$DEBUG" -eq 1 ]]; then
  set -x
fi

rsconnect deploy manifest "$REPO_ROOT/manifest.json" \
  --name "$RSCONNECT_ACCOUNT" \
  --title "$RSCONNECT_TITLE" \
  --app-id "$RSCONNECT_APP_ID" \
  -v
