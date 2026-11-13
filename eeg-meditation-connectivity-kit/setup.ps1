# Usage (PowerShell):  .\setup.ps1 -Source "C:\path\to\matlab\project" -Holder "Your Name"
param(
  [Parameter(Mandatory)] [string]$Source,
  [Parameter(Mandatory)] [string]$Holder
)
$ErrorActionPreference = 'Stop'
$Kit   = $PSScriptRoot
$Dest  = Join-Path ([Environment]::GetFolderPath('Desktop')) 'eeg-meditation-connectivity'
$Limit = 10MB
$Excluded = @('results','results_clean','results_state','data','FULL_BACKUP','.git')

if (-not (Get-Command git -ErrorAction SilentlyContinue)) { throw "git is not installed." }
if (-not (git config user.name) -or -not (git config user.email)) {
  throw 'Set your git identity first: git config --global user.name "Your Name"; git config --global user.email "you@example.com"'
}
if (-not (Test-Path $Source -PathType Container)) { throw "Source folder not found: $Source" }
$Source = (Resolve-Path $Source).Path.TrimEnd('\')
if ((Test-Path $Dest) -and (Get-ChildItem $Dest -Force)) { throw "$Dest exists and is not empty. Move it aside first." }
New-Item -ItemType Directory -Force (Join-Path $Dest 'src') | Out-Null

Write-Host "== Copying .m files from $Source =="
$copied = 0; $skipped = 0
Get-ChildItem $Source -Recurse -File -Filter *.m | ForEach-Object {
  $rel = $_.FullName.Substring($Source.Length + 1)
  $parts = $rel -split '[\\/]'
  if ($parts[0..($parts.Count-2)] | Where-Object { $Excluded -contains $_ }) { return }
  if ($_.Length -gt $Limit) { Write-Host "  skip (>10 MB): $rel"; $script:skipped++; return }
  $target = Join-Path (Join-Path $Dest 'src') $rel
  New-Item -ItemType Directory -Force (Split-Path $target) | Out-Null
  Copy-Item $_.FullName $target
  $script:copied++
}
Write-Host "  copied: $copied   skipped: $skipped"
if ($copied -eq 0) { throw "No .m files found. Check the source path." }

Copy-Item "$Kit\templates\README.md"     "$Dest\README.md"
Copy-Item "$Kit\templates\gitignore.txt" "$Dest\.gitignore"
(Get-Content "$Kit\templates\LICENSE.txt" -Raw).Replace('__YEAR__', (Get-Date).Year).Replace('__HOLDER__', $Holder) |
  Set-Content "$Dest\LICENSE" -NoNewline
New-Item -ItemType Directory -Force "$Dest\.githooks" | Out-Null
Copy-Item "$Kit\templates\pre-commit.txt" "$Dest\.githooks\pre-commit"

Push-Location $Dest
try {
  git init -q -b main
  git config core.hooksPath .githooks
  git add -A
  git update-index --chmod=+x .githooks/pre-commit

  Write-Host "`n== Files staged for the first commit =="
  git ls-files | ForEach-Object { '{0,10}  {1}' -f (Get-Item $_).Length, $_ }
  $bad = git ls-files | Where-Object { $_ -notmatch '\.m$|^README\.md$|^LICENSE$|^\.gitignore$|^\.githooks/pre-commit$' }
  if ($bad) { Write-Host "`n!! Unexpected files staged - NOT committing:"; $bad; exit 1 }

  git commit -q -m "Initial commit: EEG meditation connectivity MATLAB pipeline"
  Write-Host "`n== Committed $((git ls-files).Count) files =="
  git log --oneline -1
} finally { Pop-Location }
