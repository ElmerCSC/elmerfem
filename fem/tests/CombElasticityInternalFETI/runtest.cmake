include(test_macros)

# Cut the comb lengthwise: one partition gets the base, the other all four
# teeth, which touch each other only through the base. That subdomain is then
# in four disconnected pieces and floats, so its kernel is four sets of rigid
# body modes and it is fixed with Lagrange coefficients.
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 comb.grd)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 comb -partcell 1 2 1 -nooverwrite)

file(READ "${TEST_SOURCE}/linsys.sif" LINSYS_TEMPLATE)

string(REPLACE "Total Feti = Logical False" "Total Feti = Logical False"
       LINSYS "${LINSYS_TEMPLATE}")
file(WRITE linsys.sif "${LINSYS}")
RUN_ELMER_TEST()

string(REPLACE "Total Feti = Logical False" "Total Feti = Logical True"
       LINSYS "${LINSYS_TEMPLATE}")
file(WRITE linsys.sif "${LINSYS}")
RUN_ELMER_TEST()

# And once forcing the Lagrange coefficients, so the run goes through the
# augmented matrix with three DOFs per node whatever the geometry decides.
string(REPLACE "Total Feti = Logical False" "Total Feti = Logical True"
       LINSYS "${LINSYS_TEMPLATE}")
file(WRITE linsys.sif "${LINSYS}\n  Feti Fixing Using L.C. = Logical True\n")
RUN_ELMER_TEST()
