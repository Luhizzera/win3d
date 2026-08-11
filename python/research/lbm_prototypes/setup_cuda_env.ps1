# Roda uma vez, depois de "pip install -r requirements.txt", para deixar o
# numba.cuda capaz de achar a NVVM/cudart instaladas via pip (numba espera o
# layout de um toolkit CUDA "de verdade" ou conda; os pacotes nvidia-*-cu12
# do pip ficam espalhados em pastas separadas dentro de site-packages).
#
# Cria .venv\cuda_home\{nvvm,bin} como junctions apontando pros pacotes
# certos. O lbm_cylinder_2d_cuda.py ja aponta CUDA_HOME pra essa pasta
# automaticamente, sem precisar setar variavel de ambiente manualmente.

$ErrorActionPreference = "Stop"
$venv = Join-Path $PSScriptRoot ".venv"
$root = Join-Path $venv "cuda_home"

New-Item -ItemType Directory -Force -Path $root | Out-Null

$nvvmTarget = Join-Path $venv "Lib\site-packages\nvidia\cuda_nvcc\nvvm"
$binTarget = Join-Path $venv "Lib\site-packages\nvidia\cuda_runtime\bin"

foreach ($pair in @(@("nvvm", $nvvmTarget), @("bin", $binTarget))) {
    $link = Join-Path $root $pair[0]
    $target = $pair[1]
    if (-not (Test-Path $target)) {
        throw "Pacote nao encontrado: $target (rode pip install -r requirements.txt primeiro)"
    }
    if (Test-Path $link) { Remove-Item $link -Force }
    New-Item -ItemType Junction -Path $link -Target $target | Out-Null
    Write-Host "OK: $link -> $target"
}
