param(
    [string]$SourceDir = ".\\testdata\\sonnet",
    [string]$OutputDir = ".\\examples\\sonnet",
    [switch]$OutputStep
)

$ErrorActionPreference = "Stop"

if (-not (Test-Path $SourceDir)) {
    throw "Source directory not found: $SourceDir"
}

New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null

$sonxFiles = Get-ChildItem -Recurse $SourceDir -Filter *.sonx
if (-not $sonxFiles) {
    throw "No .sonx files found under $SourceDir"
}

Write-Host "Converting $($sonxFiles.Count) Sonnet projects to $OutputDir ..."
foreach ($file in $sonxFiles) {
    $args = @(
        "run", "-p", "rem-cli", "--",
        "--project", $file.FullName,
        "--format", "sonnet19",
        "--output", $OutputDir
    )

    if ($OutputStep) {
        $args += "--output-step"
    }

    Write-Host " - $($file.Name)"
    cargo @args | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "Conversion failed for $($file.FullName)"
    }
}

Write-Host "Done. Converted files are in $OutputDir"
