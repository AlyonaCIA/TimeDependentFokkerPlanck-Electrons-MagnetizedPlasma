 FC=f90
FLAGS = -col80 -c -O3
.SUFFIXES: .o .f



EXEC1=fkrplk.x
INCLUDEF=fkrplk.h


OBJECTS= main.o\
        init.o\
        inject.o\
        bfield.o\
        dands.o\
        cntridag.o\
        tauup.o\
        bcup.o\
        eup.o\
        muup.o 
.f.o:
        $(FC) $(FLAGS) $<



$(EXEC1): $(OBJECTS)
        $(FC) -o $(EXEC1) $(OBJECTS)
$(OBJECTS): $(INCLUDEF)