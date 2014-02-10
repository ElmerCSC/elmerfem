all: buuhuu

include ../make.inc

CFLAGS += -I$(SRC)
FCFLAGS += $(MODSEARCH_OPT)$(SRC)
LIBS = -L$(SRC) -lelmerparamf -lelmerparam -L${PREFIX}/lib -lmatc

.SUFFIXES:
.SUFFIXES: .f .f90 .o .c

.f90.o:
	$(FC) -c $(FCFLAGS) $<

.c.o:
	$(CC) -c $(CFLAGS) $<


buuhuu: buuhuu.o
	$(FC) $(FCFLAGS) -o $@ $? $(LIBS)

buuhuu.o: buuhuu.f90

clean:
	rm -rf *.o *.mod buuhuu buu result test.log
