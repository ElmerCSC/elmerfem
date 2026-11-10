include(test_macros)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 winkel -partdual -metiskway ${MPIEXEC_NTASKS} -nooverwrite)

# Run the case twice, with the standard scheme and with total FETI, so both are
# covered. They differ in how the Dirichlet conditions are imposed -- eliminated
# from the local matrices, or constrained by Lagrange coefficients, in which case
# every subdomain floats -- and the solution has to come out the same either way.
# This is the case that caught total FETI leaving the copies of a constrained
# interface DOF untied, so it is worth running both ways here in particular.
#
# linsys.sif is pulled in by case.sif and is written here from the copy in the
# source tree, so the run does not depend on what a previous one left behind.
file(READ "${TEST_SOURCE}/linsys.sif" LINSYS_TEMPLATE)

string(REPLACE "Total Feti = Logical False" "Total Feti = Logical False"
       LINSYS "${LINSYS_TEMPLATE}")
file(WRITE linsys.sif "${LINSYS}")
RUN_ELMER_TEST()

string(REPLACE "Total Feti = Logical False" "Total Feti = Logical True"
       LINSYS "${LINSYS_TEMPLATE}")
file(WRITE linsys.sif "${LINSYS}")
RUN_ELMER_TEST()
