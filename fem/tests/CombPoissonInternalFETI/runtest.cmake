include(test_macros)

# Cut the comb in two along its length, so that one partition gets the base and
# the other gets all four teeth. The teeth are joined to each other only through
# the base, so that partition is a subdomain in four disconnected pieces.
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 comb.grd)
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 comb -partcell 1 2 1 -nooverwrite)

# Run both ways, as for the other FETI cases: the schemes differ in how the
# Dirichlet conditions are imposed and have to agree on the solution.
file(READ "${TEST_SOURCE}/linsys.sif" LINSYS_TEMPLATE)

string(REPLACE "Total Feti = Logical False" "Total Feti = Logical False"
       LINSYS "${LINSYS_TEMPLATE}")
file(WRITE linsys.sif "${LINSYS}")
RUN_ELMER_TEST()

string(REPLACE "Total Feti = Logical False" "Total Feti = Logical True"
       LINSYS "${LINSYS_TEMPLATE}")
file(WRITE linsys.sif "${LINSYS}")
RUN_ELMER_TEST()

# And once more with the kernel removed by Lagrange coefficients instead of by
# pinning DOFs, which solves [A z^T; z 0] rather than a pinned A. Total FETI is
# on so that the subdomains actually float and the choice has something to act
# on. Nothing else covers that branch.
string(REPLACE "Total Feti = Logical False" "Total Feti = Logical True"
       LINSYS "${LINSYS_TEMPLATE}")
file(WRITE linsys.sif "${LINSYS}\n  Feti Fixing Using L.C. = Logical True\n")
RUN_ELMER_TEST()
