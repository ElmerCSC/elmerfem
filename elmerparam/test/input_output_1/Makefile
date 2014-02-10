all: foobar

include ../make.inc

FCFLAGS += $(MODSEARCH_OPT)$(SRC)
LIBS = -L$(SRC) -lelmerparamf -lelmerparam -L${PREFIX}/lib -lmatc

.SUFFIXES:
.SUFFIXES: .f .f90 .o

.f90.o:
	$(FC) -c $(FCFLAGS) $<


foobar: foobar.o 
	$(FC) $(FCFLAGS) -o $@ $? $(LIBS)

foobar.o: foobar.f90

clean:
	rm -rf *.o *.mod foobar result save.dat test.log bar foo
