all: cmdfile_1

include ../make.inc

CFLAGS += -I${SRC}
LDFLAGS = -L${SRC} -lelmerparam -L${PREFIX}/lib -lmatc -lm

cmdfile_1: cmdfile_1.c
	${CC} ${CFLAGS} -o cmdfile_1 cmdfile_1.c $(LDFLAGS)

clean:
	rm -rf cmdfile_1 *.FOO result test.log blubb
