#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash run.sh imggenetics --data-dir data --output-dir outputs
#   bash run.sh python scripts/fmri/compute_connectivity_msdl.py --help
#
# This is a thin convenience wrapper that:
#   1) creates/activates .venv
#   2) installs requirements + this repo (editable)
#   3) runs the command you pass after it

PYTHON_BIN="${PYTHON_BIN:-python3}"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

VENV_DIR=".venv"
if [ ! -d "$VENV_DIR" ]; then
  "$PYTHON_BIN" -m venv "$VENV_DIR"
fi

# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

python -m pip install --upgrade pip >/dev/null
pip install -r requirements.txt >/dev/null
pip install -e . >/dev/null

if [ "$#" -eq 0 ]; then
  echo "[ERROR] No command provided."
  echo "Example:"
  echo "  bash run.sh imggenetics --data-dir data --output-dir outputs"
  exit 2
fi

if [ "$1" = "imggenetics" ]; then
  shift
  python -m imggenetics "$@"
else
  "$@"
fi
