@echo off
REM ============================================================================
REM CLIF-LungCx-Epi | Run all analysis scripts (Windows)
REM
REM Usage:
REM   run.bat
REM
REM Logs are saved to the output\run_<SITE>_<DATE>\ folder alongside results.
REM No patient-level PHI is printed - only aggregate counts, file paths, and
REM config info.
REM ============================================================================

setlocal enabledelayedexpansion

cd /d "%~dp0"

REM Read site name from config.json
set CONFIG_FILE=config\config.json
if not exist "%CONFIG_FILE%" (
    echo ERROR: %CONFIG_FILE% not found. Copy config_template.json and fill in your site details.
    exit /b 1
)

for /f "tokens=*" %%i in ('python -c "import json; print(json.load(open('%CONFIG_FILE%'))['site_name'])" 2^>nul') do set SITE_NAME=%%i

if "%SITE_NAME%"=="" (
    echo ERROR: Could not read site_name from %CONFIG_FILE%
    exit /b 1
)

for /f "tokens=2 delims==" %%I in ('wmic os get localdatetime /value') do set datetime=%%I
set DATE_STAMP=%datetime:~0,8%

set OUT_DIR=output\run_%SITE_NAME%_%DATE_STAMP%
if not exist "%OUT_DIR%" mkdir "%OUT_DIR%"

set LOG_FILE=%OUT_DIR%\pipeline_%SITE_NAME%_%DATE_STAMP%.log

call :log "=== CLIF-LungCx-Epi pipeline started ==="
call :log "Site: %SITE_NAME%"
call :log "Working directory: %cd%"
call :log "Output directory: %OUT_DIR%"
call :log "Log file: %LOG_FILE%"

REM Step 0: Install packages
call :log "--- Step 0: Checking packages ---"
Rscript code\00_renv_restore.R >> "%LOG_FILE%" 2>&1
if errorlevel 1 (
    call :log "ERROR: Step 0 failed. Check log for details."
    exit /b 1
)
call :log "--- Step 0: Complete ---"

REM Run full pipeline (02 sources 01 internally so all objects stay in one R process)
call :log "--- Running full pipeline (01 -> 02) ---"
Rscript code\02_federated_clusters.R >> "%LOG_FILE%" 2>&1
if errorlevel 1 (
    call :log "ERROR: Pipeline failed. Check log for details."
    exit /b 1
)
call :log "--- Pipeline complete ---"

call :log "=== Pipeline finished successfully ==="
call :log "Review log: %LOG_FILE%"
call :log "Review output: %OUT_DIR%\"

exit /b 0

:log
echo [%date% %time%] %~1
echo [%date% %time%] %~1 >> "%LOG_FILE%"
goto :eof
