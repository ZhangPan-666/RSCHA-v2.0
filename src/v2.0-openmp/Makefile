SHELL := /bin/bash
.SUFFIXES:
.DEFAULT_GOAL := all
# run.pbs initializes Intel before invoking make. Command-line/environment CXX is respected.
ifneq ($(filter default undefined,$(origin CXX)),)
CXX := $(shell for c in icpc icpx g++; do command -v $$c >/dev/null 2>&1 && { printf '%s' $$c; break; }; done)
endif
# Do not fail at parse time: clean and dry-run remain usable without the compiler module.
ifeq ($(strip $(CXX)),)
CXX := icpc
endif
ifneq ($(filter icpx icpc mpiicpx mpiicpc,$(notdir $(CXX))),)
OMPFLAGS = -qopenmp
FPFLAGS = -fp-model=precise
else
OMPFLAGS = -fopenmp
FPFLAGS = -fno-fast-math -ffp-contract=off
endif
CPPFLAGS += -I.
CXXFLAGS += -O2 -m64 -std=c++14 $(OMPFLAGS) $(FPFLAGS) -Wall -Wextra
LDFLAGS += $(OMPFLAGS)
PROGS = Work1_RSCHAInversion Work2_RSCHAForward Work3_LocalizeSHModel
HEADERS = Basis\ Function.h CanonicalResidualLegendre.h CoordinateTransformation.h \
  CRTable.h CRNearFallback.h RSCHACompatibility.h Inversion.h Legendre.h LocalizeSH2RSCHA.h \
  MPIParallel.h R_SCHA.h RSCHARoots.h RSCHAConsole.h
HEADERS += $(wildcard SHRSCHA/*.h) $(wildcard SHRSCHA/*.cpp)
.PHONY: all clean check-table check-compiler
all: $(PROGS)
check-compiler:
	@command -v "$(CXX)" >/dev/null 2>&1 || { echo "Cannot find C++ compiler $(CXX); submit with sbatch run.pbs, which loads the Intel environment."; exit 3; }
	@"$(CXX)" --version >/dev/null 2>&1 || { echo "C++ compiler/MPI wrapper backend is unavailable: $(CXX)"; exit 3; }
check-table:
	sha256sum -c cr.bin.sha256
$(PROGS): %: %.cpp $(HEADERS) | check-compiler check-table
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) "$<" -o "$@" $(LDFLAGS) $(LDLIBS)
clean:
	rm -f $(PROGS) *.o *.obj *.pdb *.dwf *.exe
