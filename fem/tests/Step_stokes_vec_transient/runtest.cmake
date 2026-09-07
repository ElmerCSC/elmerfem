include(test_macros)

# The Step mesh of Step_stokes_vec, run in time on "p:2 b:1" and with a slip
# condition on the walls. See transient.sif for what the slip is there for: it
# is what makes the solver assemble a boundary matrix at all, and the boundary
# path on a p-basis is the half of this case that the steady variants next door
# do not reach.
execute_process(COMMAND ${ELMERGRID_BIN} 1 2 Step)

RUN_ELMER_TEST()
