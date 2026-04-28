#!/usr/bin/env bash
set -euo pipefail

VENV_PATH=".venv"
KERNEL_NAME="bathy2strat"
KERNEL_DISPLAY_NAME="Python (.venv bathy2strat)"
USE_LOCK=0

if [[ "${1:-}" == "--lock" ]]; then
  USE_LOCK=1
fi

if [[ ! -d "$VENV_PATH" ]]; then
  echo "Creating virtual environment at $VENV_PATH"
  python -m venv "$VENV_PATH"
fi

PYTHON_EXE="$VENV_PATH/bin/python"
if [[ ! -x "$PYTHON_EXE" ]]; then
  echo "Python executable not found at $PYTHON_EXE" >&2
  exit 1
fi

REQUIREMENTS_FILE="requirements.txt"
if [[ "$USE_LOCK" -eq 1 ]]; then
  REQUIREMENTS_FILE="requirements-lock.txt"
fi

if [[ ! -f "$REQUIREMENTS_FILE" ]]; then
  echo "Requirements file not found: $REQUIREMENTS_FILE" >&2
  exit 1
fi

echo "Installing dependencies from $REQUIREMENTS_FILE"
"$PYTHON_EXE" -m pip install --upgrade pip
"$PYTHON_EXE" -m pip install -r "$REQUIREMENTS_FILE"

echo "Registering Jupyter kernel '$KERNEL_DISPLAY_NAME'"
"$PYTHON_EXE" -m ipykernel install --user --name "$KERNEL_NAME" --display-name "$KERNEL_DISPLAY_NAME"

echo "Verifying core imports"
"$PYTHON_EXE" -c "import numpy, pandas, scipy, matplotlib, netCDF4, pyproj, geopandas; print('Environment OK')"

echo
echo "Setup complete."
echo "In VS Code notebooks, select kernel: $KERNEL_DISPLAY_NAME"
