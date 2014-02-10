TESTS = run-writetest run-readtest
check_PROGRAMS = writetest readtest

run-writetest: writetest.sh
	cp writetest.sh run-writetest
	chmod a+x run-writetest

run-readtest: readtest.sh
	cp readtest.sh run-readtest
	chmod a+x run-readtest

AM_FCFLAGS = @INCLUDE_MODULE_FLAG@..

writetest_SOURCES = writetest.f90
writetest_LDADD = ../libbinio.a

readtest_SOURCES = readtest.f90
readtest_LDADD = ../libbinio.a

EXTRA_DIST = writetest.sh readtest.sh writetest.checkL writetest.checkB

clean-local:
	$(RM)  writetest.out
