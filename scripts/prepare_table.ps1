$ErrorActionPreference = 'Stop'
$releaseRoot = Split-Path -Parent $PSScriptRoot
$tablePath = Join-Path $releaseRoot 'large_assets\cr.bin'
$checksumPath = Join-Path $releaseRoot 'large_assets\cr.bin.sha256'
$expectedHash = ((Get-Content -LiteralPath $checksumPath -Raw).Trim() -split '\s+')[0]
$stream = [System.IO.File]::OpenRead($tablePath)
$sha256 = [System.Security.Cryptography.SHA256]::Create()
try {
    $actualHash = [System.BitConverter]::ToString($sha256.ComputeHash($stream)).Replace('-', '')
} finally {
    $stream.Dispose()
    $sha256.Dispose()
}
if ($actualHash -ine $expectedHash) { throw 'cr.bin SHA-256 verification failed.' }
foreach ($variant in @('v2.0-openmp', 'v2.0-mpi')) {
    Copy-Item -LiteralPath $tablePath -Destination (Join-Path $releaseRoot "src\$variant\cr.bin") -Force
}
Write-Host 'Verified and prepared cr.bin for both source trees.'
