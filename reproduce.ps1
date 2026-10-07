param(
    [string]$PythonExe = 'python',
    [switch]$SkipDownload
)
$ErrorActionPreference = 'Stop'
$projectRoot = $PSScriptRoot
Push-Location -LiteralPath $projectRoot
try {
    & $PythonExe -m venv .venv
    if ($LASTEXITCODE -ne 0) { throw 'Could not create the Python environment.' }
    $analysisPython = Join-Path $projectRoot '.venv/Scripts/python.exe'
    & $analysisPython -m pip install -r requirements-lock.txt
    if ($LASTEXITCODE -ne 0) { throw 'Dependency installation failed.' }
    $stages = @('02_prepare.py','03_model.py','04_validate.py','05_notebook.py')
    if (-not $SkipDownload) { $stages = @('01_download.py') + $stages }
    foreach ($stage in $stages) {
        & $analysisPython (Join-Path 'scripts' $stage)
        if ($LASTEXITCODE -ne 0) { throw "Analysis stopped at $stage" }
    }
    Write-Output 'Completed. Open Research_Report.html or Erlotinib_Response_Study.ipynb.'
}
finally { Pop-Location }
