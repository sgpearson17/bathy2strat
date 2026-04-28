param(
    [string]$VenvPath = ".venv",
    [string]$KernelName = "bathy2strat",
    [string]$KernelDisplayName = "Python (.venv bathy2strat)",
    [switch]$UseLockFile
)

$ErrorActionPreference = "Stop"

if (-not (Test-Path $VenvPath)) {
    Write-Host "Creating virtual environment at $VenvPath"
    python -m venv $VenvPath
}

$pythonExe = Join-Path $VenvPath "Scripts/python.exe"
if (-not (Test-Path $pythonExe)) {
    throw "Python executable not found at $pythonExe"
}

$requirementsFile = if ($UseLockFile) { "requirements-lock.txt" } else { "requirements.txt" }
if (-not (Test-Path $requirementsFile)) {
    throw "Requirements file not found: $requirementsFile"
}

Write-Host "Installing dependencies from $requirementsFile"
& $pythonExe -m pip install --upgrade pip
& $pythonExe -m pip install -r $requirementsFile

Write-Host "Registering Jupyter kernel '$KernelDisplayName'"
& $pythonExe -m ipykernel install --user --name $KernelName --display-name $KernelDisplayName

Write-Host "Verifying core imports"
& $pythonExe -c "import numpy, pandas, scipy, matplotlib, netCDF4, pyproj, geopandas; print('Environment OK')"

Write-Host ""
Write-Host "Setup complete."
Write-Host "In VS Code notebooks, select kernel: $KernelDisplayName"
