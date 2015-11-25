FC := mpif90
CC := mpicc

build: test

feti4i_mod.o: ../feti4i_mod.F90
	${FC} -c $^

feti4i.mod: 	../feti4i_mod.F90 feti4i_mod.o
	@true

test.o: test.F90 feti4i.mod
	${FC} -c -g test.F90 -o $@ 

feti4i_dummy.o: feti4i_dummy.c
	${CC} -c -g $^

test: test.o feti4i_dummy.o 
	${LINK.F} test.o feti4i_dummy.o -o $@ 
  
clean:
	${RM} test.o feti4i_dummy.o feti4i_mod.o feti4i.mod

.PHONY: build clean

