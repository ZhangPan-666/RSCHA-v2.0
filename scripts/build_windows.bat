@echo off
setlocal EnableExtensions
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0prepare_table.ps1"
if errorlevel 1 exit /b 1
rem Public build wrapper for Intel oneAPI C++ on Windows.
if not defined ONEAPI_ROOT set "ONEAPI_ROOT=C:\Program Files (x86)\Intel\oneAPI"
if not exist "%ONEAPI_ROOT%\setvars.bat" (
  echo Set ONEAPI_ROOT to the Intel oneAPI installation directory.
  exit /b 2
)
call "%ONEAPI_ROOT%\setvars.bat" intel64 >nul
if errorlevel 1 exit /b 3
set "ROOT=%~dp0.."
set "OMP=%ROOT%\src\v2.0-openmp"
set "MPI=%ROOT%\src\v2.0-mpi"
icx /nologo /utf-8 /std:c++17 /O2 /fp:precise /Qopenmp /EHsc /W4 /clang:-Wall /clang:-Wextra /I"%OMP%" "%OMP%\Work1_RSCHAInversion.cpp" /Fe:"%OMP%\Work1_RSCHAInversion.exe"
if errorlevel 1 exit /b 4
icx /nologo /utf-8 /std:c++17 /O2 /fp:precise /Qopenmp /EHsc /W4 /clang:-Wall /clang:-Wextra /I"%OMP%" "%OMP%\Work2_RSCHAForward.cpp" /Fe:"%OMP%\Work2_RSCHAForward.exe"
if errorlevel 1 exit /b 5
icx /nologo /utf-8 /std:c++17 /O2 /fp:precise /Qopenmp /EHsc /W4 /clang:-Wall /clang:-Wextra /I"%OMP%" "%OMP%\Work3_LocalizeSHModel.cpp" /Fe:"%OMP%\Work3_LocalizeSHModel.exe"
if errorlevel 1 exit /b 6
if not defined I_MPI_ROOT (
  echo OpenMP build completed. Set I_MPI_ROOT to build MPI drivers.
  exit /b 0
)
for %%F in (Work1_RSCHAInversion Work2_RSCHAForward Work3_LocalizeSHModel) do (
  icx /nologo /utf-8 /std:c++17 /O2 /fp:precise /Qopenmp /EHsc /W4 /clang:-Wall /clang:-Wextra /I"%I_MPI_ROOT%\include" /I"%MPI%" "%MPI%\%%F.cpp" /Fe:"%MPI%\%%F.exe" /link /LIBPATH:"%I_MPI_ROOT%\lib" impi.lib
  if errorlevel 1 exit /b 10
)
echo OpenMP and MPI builds completed.
