execute_process(COMMAND ${ELMERGRID_BIN} 14 2 sector.msh -autoclean)
execute_process(COMMAND ${ELMERGRID_BIN} sector_extrude.eg)
execute_process(COMMAND ${ELMERGRID_BIN} 2 2 sector3d -partdual -metis ${MPIEXEC_NTASKS} 3 -connect 5 6 -partrbc 2 -partlayers 0 -halor -nooverwrite)
RUN_ELMER_TEST()
