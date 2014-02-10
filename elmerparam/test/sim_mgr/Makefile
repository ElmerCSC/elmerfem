all: runopt

include ../make.inc

FCFLAGS += $(MODSEARCH_OPT)$(SRC)
LIBS = -L$(SRC) -lelmerparamf -lelmerparam -L${PREFIX}/lib -lmatc

.SUFFIXES:
.SUFFIXES: .f .f90 .o

.f90.o:
	$(FC) -c $(FCFLAGS) $<


runopt: runopt.o 
	$(FC) $(FCFLAGS) -o $@ $? $(LIBS)

runopt.o: runopt.f90

clean:
	rm -rf *.o *.mod runopt *.dat mesh scalars.dat scalars.dat.names test.log
