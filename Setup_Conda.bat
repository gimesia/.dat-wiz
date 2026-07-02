@echo off
REM ===================================================================
REM  One-time conda setup for the K-space Viewer.
REM  Double-click this ONCE. It:
REM    1. creates a conda env "kspace-viewer" (Python 3.11)
REM    2. installs the dependencies from requirements.txt
REM    3. records that env's python path so Run_KSpace_Viewer.bat
REM       uses it automatically (no manual configuration needed).
REM  After this finishes, just double-click Run_KSpace_Viewer.bat.
REM ===================================================================
setlocal
set "HERE=%~dp0"
set "ENV_NAME=kspace-viewer"

REM --- Locate conda ---------------------------------------------------
set "CONDA="
where conda >nul 2>nul && set "CONDA=conda"
if not defined CONDA (
    for %%C in (
        "%USERPROFILE%\miniconda3\Scripts\conda.exe"
        "%USERPROFILE%\anaconda3\Scripts\conda.exe"
        "%USERPROFILE%\AppData\Local\miniconda3\Scripts\conda.exe"
        "%USERPROFILE%\AppData\Local\Continuum\miniconda3\Scripts\conda.exe"
        "%ProgramData%\miniconda3\Scripts\conda.exe"
        "%ProgramData%\Anaconda3\Scripts\conda.exe"
    ) do (
        if exist %%C set "CONDA=%%~C"
    )
)
if not defined CONDA (
    echo.
    echo ERROR: Could not find conda.
    echo   - Install Miniconda: https://docs.conda.io/en/latest/miniconda.html
    echo   - Or run this script from an "Anaconda Prompt".
    echo.
    pause
    exit /b 1
)
echo Using conda: %CONDA%

REM --- Create the environment only if it doesn't already exist -------
call "%CONDA%" env list | findstr /I /C:"\envs\%ENV_NAME%" >nul
if not errorlevel 1 (
    echo.
    echo Conda env "%ENV_NAME%" already exists - reusing it.
) else (
    echo.
    echo Creating conda env "%ENV_NAME%" (Python 3.11)...
    call "%CONDA%" create -n %ENV_NAME% python=3.11 -y
    if errorlevel 1 (
        echo.
        echo ERROR: failed to create the conda env.
        pause
        exit /b 1
    )
)

REM --- Install dependencies ------------------------------------------
echo.
echo Installing dependencies (this can take a few minutes)...
call "%CONDA%" run -n %ENV_NAME% python -m pip install -r "%HERE%requirements.txt"
if errorlevel 1 (
    echo.
    echo ERROR: pip install failed.
    pause
    exit /b 1
)

REM --- Record the env's python path for the launcher -----------------
set "PYPATH="
for /f "delims=" %%p in ('"%CONDA%" run -n %ENV_NAME% python -c "import sys;print(sys.executable)"') do set "PYPATH=%%p"
if not defined PYPATH (
    echo.
    echo ERROR: could not determine the env's python path.
    pause
    exit /b 1
)
> "%HERE%env_python.txt" echo %PYPATH%

echo.
echo ===================================================================
echo  Setup complete.
echo  Environment: %ENV_NAME%
echo  Python:      %PYPATH%
echo.
echo  You can now double-click  Run_KSpace_Viewer.bat  to start the tool.
echo ===================================================================
echo.
pause
