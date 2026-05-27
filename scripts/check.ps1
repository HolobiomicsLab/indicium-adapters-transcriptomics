# Run the same checks as the GitHub Actions CI: ruff lint + pytest with
# coverage. Use this before every push so a lint or test regression is
# caught locally — pytest alone does not exercise ruff, and the CI runs
# both. Aborts on the first failure.
#
# Usage:  ./scripts/check.ps1

$ErrorActionPreference = "Stop"

Write-Host "[check] ruff check ." -ForegroundColor Cyan
py -3.13 -m ruff check .
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "[check] pytest" -ForegroundColor Cyan
py -3.13 -m pytest
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "[check] All CI checks passed locally." -ForegroundColor Green
