# =============================================================================
# Makefile for the Time-Dependent Fokker-Planck Solver (fkrplk)
#
# Compiles all Fortran 77 source files and links them into the executable.
# Uses gfortran by default. Override with: make FC=ifort FLAGS="..."
# =============================================================================

FC       = gfortran
FLAGS    = -ffixed-line-length-80 -c -O3
LDFLAGS  =
EXEC     = fkrplk.x
HEADER   = fkrplk.h

SOURCES  = main.f init.f inject.f bfield.f dands.f cntridag.f \
           tauup.f bcup.f eup.f muup.f

OBJECTS  = $(SOURCES:.f=.o)

# Default output directory (override with: make OUTDIR=../../output)
OUTDIR   = ../../output

# --- Build rules ---

.SUFFIXES: .o .f

.f.o:
	$(FC) $(FLAGS) $<

all: $(EXEC)

$(EXEC): $(OBJECTS)
	$(FC) -o $(EXEC) $(OBJECTS) $(LDFLAGS)
	@echo "Build successful: $(EXEC)"

$(OBJECTS): $(HEADER)

# --- Run the simulation ---
run: $(EXEC)
	@mkdir -p $(OUTDIR)
	./$(EXEC)
	@echo "Simulation complete."

# --- Clean build artifacts ---
clean:
	rm -f $(OBJECTS) $(EXEC)
	@echo "Cleaned."

.PHONY: all clean run
