all: buuhuu

include ../make.inc

CFLAGS += -I$(SRC)
FCFLAGS += $(MODSEARCH_OPT)$(SRC)
LIBS = -L$(SRC) -lelmerparam -L${PREFIX}/lib -lmatc -lm

.SUFFIXES:
.SUFFIXES: .f .f90 .o .c

.f90.o:
	$(FC) -c $(FCFLAGS) $<

.c.o:
	$(CC) -c $(CFLAGS) $<


buuhuu: buuhuu.o
	$(CC) $(CFLAGS) -o $@ $? $(LIBS)

buuhuu.o: buuhuu.c

clean:
	rm -rf *.o *.mod buuhuu result buu result test.log
