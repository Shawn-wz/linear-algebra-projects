#Requires -Version 7.0
<#
.SYNOPSIS
    Test harness for linear-equations-solver.cpp.

.DESCRIPTION
    Runs every case file through the solver and checks three independent things.

    1. Expected classification.  Computed here by an exact integer oracle
       (fraction-free Bareiss elimination over System.Numerics.BigInteger), so
       nothing is hard coded and a new case needs only a data file.  Cases with
       decimal coefficients are scaled to integers first, which leaves the rank
       and the consistency of the system unchanged.

    2. Printed classification: unique / infinitely many / no solution.

    3. Printed solution.
         unique    x must satisfy A x = b.  With rank(A) = n and a consistent
                   system that solution is unique, so the residual alone pins
                   it down.
         infinite  p must satisfy A p = b, every printed basis vector w must
                   satisfy A w = 0, there must be exactly n - rank of them, and
                   they must be linearly independent - together that means the
                   printed set really is the whole solution set.
         none      nothing further to check.

    Residuals are reported relative to the data magnitude, which makes the
    achieved numerical precision visible.

    Case file format
        line 1      <rows> <cols>
        following   <rows> lines of <cols> coefficients then the constant
        '#'         starts a comment line, blank lines are ignored

.PARAMETER Filter
    Filesystem wildcard (not a regular expression) selecting the case files,
    for example '1*' or '*unique*'.

.PARAMETER ScaleCheck
    Also re-run each case with A and b multiplied by each value in -Scales.
    That leaves the solution set unchanged, so the parsed solution must not
    change either; it is a cheap probe for tolerance/scale problems.  Off by
    default because an absolute eps is not scale invariant.

.EXAMPLE
    pwsh -NoProfile -File tests/run-tests.ps1

.EXAMPLE
    pwsh -NoProfile -File tests/run-tests.ps1 -ScaleCheck -Filter '1*'
#>
[CmdletBinding()]
param(
    [string]   $Compiler = 'C:\msys64\ucrt64\bin\g++.exe',
    [string]   $Source   = 'linear-equations-solver.cpp',
    [string]   $CaseDir  = 'tests',
    [string]   $Filter   = '*.txt',
    [double]   $Tol      = 1e-6,
    [switch]   $ScaleCheck,
    [double[]] $Scales   = @(1e-6, 1e6)
)

$ErrorActionPreference = 'Stop'
$Inv  = [Globalization.CultureInfo]::InvariantCulture
$BI   = [System.Numerics.BigInteger]
$Root = Split-Path -Parent $PSScriptRoot

function Resolve-UnderRoot([string]$path) {
    if ([IO.Path]::IsPathRooted($path)) { return $path }
    return Join-Path $Root $path
}

# ---------------------------------------------------------------- case files

function Read-CaseFile([string]$path) {
    $lines = @(Get-Content -LiteralPath $path |
               Where-Object { $_.Trim() -ne '' -and -not $_.TrimStart().StartsWith('#') })
    if ($lines.Count -lt 1) { throw "empty case file: $path" }

    $head  = $lines[0].Trim() -split '\s+'
    $nRows = [int]$head[0]
    $nCols = [int]$head[1]
    if ($lines.Count -ne $nRows + 1) {
        throw "$([IO.Path]::GetFileName($path)): declares $nRows row(s) but has $($lines.Count - 1)"
    }

    $matA = @(); $vecB = @(); $scale = 0.0
    for ($i = 1; $i -le $nRows; $i++) {
        $tok = $lines[$i].Trim() -split '\s+'
        if ($tok.Count -ne $nCols + 1) {
            throw "$([IO.Path]::GetFileName($path)) line $($i + 1): expected $($nCols + 1) numbers, got $($tok.Count)"
        }
        $vals = @()
        foreach ($t in $tok) {
            $v = 0.0
            if (-not [double]::TryParse($t, [Globalization.NumberStyles]::Float, $Inv, [ref]$v)) {
                throw "$([IO.Path]::GetFileName($path)) line $($i + 1): '$t' is not a number"
            }
            $vals += $v
            $scale = [Math]::Max($scale, [Math]::Abs($v))
        }
        $matA += ,@($vals[0..($nCols - 1)])
        $vecB += $vals[$nCols]
    }
    return @{ Name = [IO.Path]::GetFileName($path); Rows = $nRows; Cols = $nCols
              A = $matA; B = $vecB; Scale = $scale }
}

# ------------------------------------------------------------- exact oracle

function Test-Integral($case, [double]$factor) {
    for ($i = 0; $i -lt $case.Rows; $i++) {
        for ($j = 0; $j -lt $case.Cols; $j++) {
            $v = $case.A[$i][$j] * $factor
            if ([Math]::Abs($v - [Math]::Round($v)) -gt 1e-6) { return $false }
        }
        $v = $case.B[$i] * $factor
        if ([Math]::Abs($v - [Math]::Round($v)) -gt 1e-6) { return $false }
    }
    return $true
}

# Scale the system up until every entry is an integer (rank and consistency are
# invariant under a positive scaling), then eliminate [A|b] exactly.
function Get-ExactClass($case) {
    $dec = 0
    while ($dec -le 12) {
        $f = [Math]::Pow(10, $dec)
        if (Test-Integral $case $f) { break }
        $dec++
    }
    if ($dec -gt 12) { throw "$($case.Name): coefficients need more than 12 decimal places" }

    $f = [Math]::Pow(10, $dec)
    $mat = @()
    for ($i = 0; $i -lt $case.Rows; $i++) {
        $row = @()
        for ($j = 0; $j -lt $case.Cols; $j++) {
            $row += $BI::Parse([string][long][Math]::Round($case.A[$i][$j] * $f), $Inv)
        }
        $row += $BI::Parse([string][long][Math]::Round($case.B[$i] * $f), $Inv)
        $mat += ,$row
    }

    $prev = $BI::One
    $rank = 0
    for ($col = 0; $col -lt $case.Cols -and $rank -lt $case.Rows; $col++) {
        $piv = -1
        for ($r = $rank; $r -lt $case.Rows; $r++) {
            if ($mat[$r][$col] -ne $BI::Zero) { $piv = $r; break }
        }
        if ($piv -lt 0) { continue }                        # no pivot in this column
        $tmp = $mat[$rank]; $mat[$rank] = $mat[$piv]; $mat[$piv] = $tmp
        for ($r = $rank + 1; $r -lt $case.Rows; $r++) {
            for ($c = $col + 1; $c -le $case.Cols; $c++) {
                $num = $mat[$r][$c] * $mat[$rank][$col] - $mat[$r][$col] * $mat[$rank][$c]
                $mat[$r][$c] = $BI::Divide($num, $prev)     # exact by Bareiss
            }
            $mat[$r][$col] = $BI::Zero
        }
        $prev = $mat[$rank][$col]
        $rank++
    }

    $consistent = $true
    for ($r = $rank; $r -lt $case.Rows; $r++) {
        if ($mat[$r][$case.Cols] -ne $BI::Zero) { $consistent = $false; break }
    }
    if (-not $consistent)          { $type = 'none' }
    elseif ($rank -eq $case.Cols)  { $type = 'unique' }
    else                           { $type = 'infinite' }
    return @{ Type = $type; Rank = $rank }
}

# --------------------------------------------------------- solver interface

function Get-Verdict([string]$text) {
    if ($text -match 'unique')          { return 'unique' }
    if ($text -match 'infinit')         { return 'infinite' }
    if ($text -match '(?i)no solution') { return 'none' }
    return 'MISSING'
}

# Every "( ... )T" group: the first is the particular solution, the rest are
# null-space basis vectors.
function Get-SolutionGroups([string]$text) {
    $groups = @()
    foreach ($m in [regex]::Matches($text, '\(\s*([^)]*?)\s*\)\s*T')) {
        $tok = ($m.Groups[1].Value).Trim() -split '\s+'
        if ($tok.Count -eq 1 -and $tok[0] -eq '') { $groups += ,@(); continue }
        $vals = @(); $good = $true
        foreach ($t in $tok) {
            $v = 0.0
            if (-not [double]::TryParse($t, [Globalization.NumberStyles]::Float, $Inv, [ref]$v)) { $good = $false; break }
            $vals += $v
        }
        if ($good) { $groups += ,$vals } else { $groups += ,@() }
    }
    return ,$groups
}

function Invoke-Solver([string]$exe, [string[]]$lines) {
    $text = (@($lines) -join "`n") + "`n"
    return ($text | & $exe 2>&1 | Out-String)
}

# --------------------------------------------------------------- arithmetic

function Get-Residual($matA, $vecB, $x, [int]$nRows, [int]$nCols) {
    $worst = 0.0
    for ($i = 0; $i -lt $nRows; $i++) {
        $s = 0.0
        for ($j = 0; $j -lt $nCols; $j++) { $s += $matA[$i][$j] * $x[$j] }
        $worst = [Math]::Max($worst, [Math]::Abs($s - $vecB[$i]))
    }
    return $worst
}

function Get-Rank($rows, [int]$nCols) {
    $mat = @(); foreach ($r in $rows) { $mat += ,([double[]]$r) }
    $rank = 0
    for ($col = 0; $col -lt $nCols -and $rank -lt $mat.Count; $col++) {
        $sel = -1; $best = 1e-9
        for ($r = $rank; $r -lt $mat.Count; $r++) {
            if ([Math]::Abs($mat[$r][$col]) -gt $best) { $best = [Math]::Abs($mat[$r][$col]); $sel = $r }
        }
        if ($sel -lt 0) { continue }
        $tmp = $mat[$rank]; $mat[$rank] = $mat[$sel]; $mat[$sel] = $tmp
        for ($r = $rank + 1; $r -lt $mat.Count; $r++) {
            $f = $mat[$r][$col] / $mat[$rank][$col]
            if ($f -eq 0) { continue }
            for ($c = $col; $c -lt $nCols; $c++) { $mat[$r][$c] -= $f * $mat[$rank][$c] }
        }
        $rank++
    }
    return $rank
}

# ------------------------------------------------------------------- driver

$sourcePath = Resolve-UnderRoot $Source
$casePath   = Resolve-UnderRoot $CaseDir
$exePath    = Join-Path $env:TEMP 'linear-equations-solver-tests.exe'

Write-Host "Building $([IO.Path]::GetFileName($sourcePath)) ..."
$buildLog = & $Compiler '-std=c++17' '-Wall' '-Wextra' '-O2' $sourcePath '-o' $exePath 2>&1
if ($LASTEXITCODE -ne 0) {
    $buildLog | ForEach-Object { Write-Host $_ }
    Write-Error 'build failed'
    exit 2
}
$warnings = @($buildLog | Where-Object { $_ -match 'warning:' })
if ($warnings.Count) {
    Write-Host "  $($warnings.Count) compiler warning(s):" -ForegroundColor Yellow
    $warnings | ForEach-Object { Write-Host "    $_" -ForegroundColor Yellow }
}

$cases = @(Get-ChildItem -LiteralPath $casePath -Filter $Filter -File | Sort-Object Name)
if ($cases.Count -eq 0) { Write-Error "no case files matched $Filter in $casePath"; exit 2 }

$results = @()
foreach ($file in $cases) {
    $case   = Read-CaseFile $file.FullName
    $expect = Get-ExactClass $case
    $text   = Invoke-Solver $exePath (Get-Content -LiteralPath $file.FullName)
    $got    = Get-Verdict $text
    $groups = Get-SolutionGroups $text

    # tolerance scaled to the data, so it means the same thing for coefficients
    # of 1e-4 and for coefficients of 1e6
    $unit   = [Math]::Max(1.0, $case.Scale) * [Math]::Max(1, $case.Cols)
    $tolAbs = $Tol * $unit

    $faults = @()
    $relRes = $null

    if ($got -ne $expect.Type) {
        $faults += "classification: expected '$($expect.Type)', got '$got'"
    }
    elseif ($expect.Type -eq 'unique') {
        if ($groups.Count -lt 1)                  { $faults += 'no solution vector printed' }
        elseif ($groups[0].Count -ne $case.Cols)  { $faults += "solution has $($groups[0].Count) component(s), want $($case.Cols)" }
        else {
            $res    = Get-Residual $case.A $case.B $groups[0] $case.Rows $case.Cols
            $relRes = $res / $unit
            if ($res -gt $tolAbs) { $faults += "A x != b (relative residual $($relRes.ToString('E2', $Inv)))" }
        }
    }
    elseif ($expect.Type -eq 'infinite') {
        if ($groups.Count -lt 2) {
            $faults += "expected a particular solution plus $($case.Cols - $expect.Rank) basis vector(s), got $($groups.Count) group(s)"
        }
        else {
            $part  = $groups[0]
            $basis = @($groups[1..($groups.Count - 1)])
            if ($part.Count -ne $case.Cols) { $faults += "particular solution has $($part.Count) component(s), want $($case.Cols)" }
            else {
                $res    = Get-Residual $case.A $case.B $part $case.Rows $case.Cols
                $relRes = $res / $unit
                if ($res -gt $tolAbs) { $faults += "A p != b (relative residual $($relRes.ToString('E2', $Inv)))" }
            }
            foreach ($w in $basis) {
                if ($w.Count -ne $case.Cols) { $faults += "a basis vector has $($w.Count) component(s), want $($case.Cols)"; continue }
                $r0 = Get-Residual $case.A (New-Object double[] $case.Rows) $w $case.Rows $case.Cols
                if ($r0 -gt $tolAbs) { $faults += "basis vector outside the null space (residual $($r0.ToString('E2', $Inv)))" }
            }
            $want = $case.Cols - $expect.Rank
            if ($basis.Count -ne $want) {
                $faults += "solution set has the wrong dimension: $($basis.Count) basis vector(s), want $want"
            }
            elseif ((Get-Rank $basis $case.Cols) -ne $basis.Count) {
                $faults += 'basis vectors are linearly dependent'
            }
        }
    }

    $results += [pscustomobject]@{
        Name   = $case.Name
        Expect = $expect.Type
        Rank   = $expect.Rank
        Got    = $got
        RelRes = $relRes
        Faults = $faults
        Output = ($text.Trim() -replace "`r?`n", ' ~ ')
    }
}

# ------------------------------------------------------------------- report

$fail = @($results | Where-Object { $_.Faults.Count -gt 0 })

foreach ($r in $results) {
    $status = if ($r.Faults.Count -eq 0) { 'PASS' } else { 'FAIL' }
    $res    = if ($null -ne $r.RelRes) { "rel.resid $($r.RelRes.ToString('E1', $Inv))" } else { '' }
    Write-Host ("  {0}  {1,-30} {2,-9} rank {3,-3} {4}" -f $status, $r.Name, $r.Expect, $r.Rank, $res)
    if ($r.Faults.Count) {
        foreach ($f in $r.Faults) { Write-Host "          - $f" -ForegroundColor Red }
        Write-Host "          printed: $($r.Output)" -ForegroundColor DarkGray
    }
}

$measured = @($results | Where-Object { $null -ne $_.RelRes })
if ($measured.Count) {
    $worst = $measured | Sort-Object RelRes -Descending | Select-Object -First 1
    Write-Host ""
    Write-Host ("precision: worst relative residual {0} on {1} (tolerance {2})" -f `
        $worst.RelRes.ToString('E1', $Inv), $worst.Name, $Tol.ToString('E0', $Inv))
}

# ---------------------------------------------------- optional scale probe

$scaleFaults = 0
if ($ScaleCheck) {
    Write-Host ""
    Write-Host "scale probe (A,b multiplied by c; the solution must not change):"
    foreach ($file in $cases) {
        $case = Read-CaseFile $file.FullName
        $base = Get-SolutionGroups (Invoke-Solver $exePath (Get-Content -LiteralPath $file.FullName))
        foreach ($c in $Scales) {
            $lines = @("$($case.Rows) $($case.Cols)")
            for ($i = 0; $i -lt $case.Rows; $i++) {
                $vals = @()
                for ($j = 0; $j -lt $case.Cols; $j++) { $vals += ($case.A[$i][$j] * $c).ToString('R', $Inv) }
                $vals += ($case.B[$i] * $c).ToString('R', $Inv)
                $lines += ($vals -join ' ')
            }
            $scaled = Get-SolutionGroups (Invoke-Solver $exePath $lines)
            $why = $null
            if ($scaled.Count -ne $base.Count) {
                $why = "groups $($base.Count) -> $($scaled.Count)"
            }
            else {
                for ($k = 0; $k -lt $base.Count -and -not $why; $k++) {
                    if ($scaled[$k].Count -ne $base[$k].Count) { $why = "group $k length changed"; break }
                    for ($j = 0; $j -lt $base[$k].Count; $j++) {
                        $t = 1e-4 * [Math]::Max(1.0, [Math]::Abs($base[$k][$j]))
                        if ([Math]::Abs($scaled[$k][$j] - $base[$k][$j]) -gt $t) {
                            $why = "group $k component $($j + 1): $($base[$k][$j]) -> $($scaled[$k][$j])"
                            break
                        }
                    }
                }
            }
            if ($why) {
                $scaleFaults++
                Write-Host ("  FAIL  {0,-30} c={1,-7} {2}" -f $case.Name, $c, $why) -ForegroundColor Red
            }
        }
    }
    if ($scaleFaults -eq 0) { Write-Host "  every scale reproduced the baseline solution" }
}

# ------------------------------------------------------------------ summary

Write-Host ""
Write-Host ("{0} case(s): {1} passed, {2} failed" -f $results.Count, ($results.Count - $fail.Count), $fail.Count)
if ($ScaleCheck) {
    $probeTotal = $cases.Count * $Scales.Count
    Write-Host ("scale probe: {0} of {1} checks changed the solution" -f $scaleFaults, $probeTotal)
}
Remove-Item $exePath -ErrorAction SilentlyContinue
if ($fail.Count -gt 0 -or $scaleFaults -gt 0) { exit 1 }
exit 0
