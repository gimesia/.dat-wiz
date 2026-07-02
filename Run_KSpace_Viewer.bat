@echo off
REM ===================================================================
REM  Double-click launcher for the K-space Viewer / NUFFT recon tool.
REM  Opens a browse dialog for the .dat (the .seq is auto-located) and
REM  keeps this window open so you can read the startup log / errors.
REM
REM  Which Python it uses (first match wins):
REM    1. KSPACE_PYTHON, if you set it manually.
REM    2. env_python.txt next to this file, written by Setup_Conda.bat.
REM    3. `python` on your PATH.
REM  Run Setup_Conda.bat once and step 2 is handled for you.
REM ===================================================================
setlocal
set "HERE=%~dp0"

if not defined KSPACE_PYTHON (
    if exist "%HERE%env_python.txt" (
        set /p KSPACE_PYTHON=<"%HERE%env_python.txt"
    ) else (
        set "KSPACE_PYTHON=python"
    )
)

"%KSPACE_PYTHON%" "%HERE%main.py" %*

echo.
pause
