TESTS = masterscript.sh

TESTDIRS = cmdfile_1 \
	input_output_1 \
	input_output_2 \
	input_output_3 \
	matc \
	shell_1 \
	shell_2 \
	sim_mgr

EXTRA_DIST = masterscript.sh $(TESTDIRS) make.inc.in make.inc

clean-local:
	for d in $(TESTDIRS) ; do \
		cd $$d; \
		make clean; \
		cd ..; \
	done

dist-hook: 
	cd $(distdir); \
	for d in $(TESTDIRS) ; do \
		cd $$d; \
		make clean; \
		cd ..; \
	done; \
	rm make.inc
