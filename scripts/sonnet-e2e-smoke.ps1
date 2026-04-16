param(
    [string]$Fixture = ".\\tests\\fixtures\\interdigital_capacitor.sonx",
    [string]$WorkDir = ".\\output\\sonnet_e2e",
    [double]$FreqHz = 1.0e9,
    [switch]$VerboseSolver
)

$ErrorActionPreference = "Stop"

if (-not (Test-Path $Fixture)) {
    throw "Fixture not found: $Fixture"
}

New-Item -ItemType Directory -Force -Path $WorkDir | Out-Null

$stem = [System.IO.Path]::GetFileNameWithoutExtension($Fixture)
$configPath = Join-Path $WorkDir ("$stem.json")
$meshPath = Join-Path $WorkDir ("$stem.msh")
$stepPath = Join-Path $WorkDir ("${stem}_geometry.step")
$solveOut = Join-Path $WorkDir "solve"

Write-Host "[1/3] Convert Sonnet -> REM (single-frequency)"
cargo run -p rem-cli -- --project $Fixture --format sonnet19 --output $WorkDir --freq-min $FreqHz --freq-max $FreqHz --freq-step 1.0 --output-step | Out-Null
if ($LASTEXITCODE -ne 0) {
    throw "Conversion failed"
}

if (-not (Test-Path $configPath)) {
    throw "Missing config output: $configPath"
}
if (-not (Test-Path $meshPath)) {
    throw "Missing mesh output: $meshPath"
}
if (-not (Test-Path $stepPath)) {
    throw "Missing STEP output: $stepPath"
}

$stepItem = Get-Item $stepPath
if ($stepItem.Length -le 0) {
    throw "STEP output is empty: $stepPath"
}

Write-Host "[2/3] Run REM solver"
$solveArgs = @("run", "-p", "rem-cli", "--", $configPath, "-o", $solveOut)
if ($VerboseSolver) {
    $solveArgs += "-v"
}
cargo @solveArgs
if ($LASTEXITCODE -ne 0) {
    throw "Solver failed"
}

Write-Host "[3/3] Validate solve output path"
if (-not (Test-Path $solveOut)) {
    throw "Solver output directory missing: $solveOut"
}

$postproDir = Join-Path $solveOut "postpro"
if (-not (Test-Path $postproDir)) {
    throw "Solver postpro directory missing: $postproDir"
}

$portCsv = Join-Path $postproDir "port-S.csv"
if (-not (Test-Path $portCsv)) {
    throw "Expected S-parameter CSV missing: $portCsv"
}

$csvItem = Get-Item $portCsv
if ($csvItem.Length -le 0) {
    throw "S-parameter CSV is empty: $portCsv"
}

$touchstone = Get-ChildItem $postproDir -Filter "s_params.s*p" -ErrorAction SilentlyContinue | Select-Object -First 1
if (-not $touchstone) {
    throw "Expected Touchstone file missing under: $postproDir"
}
if ($touchstone.Length -le 0) {
    throw "Touchstone file is empty: $($touchstone.FullName)"
}

Write-Host "E2E smoke succeeded"
Write-Host " - Config: $configPath"
Write-Host " - Mesh:   $meshPath"
Write-Host " - STEP:   $stepPath"
Write-Host " - Solve:  $solveOut"
Write-Host " - CSV:    $portCsv"
Write-Host " - SNP:    $($touchstone.FullName)"
