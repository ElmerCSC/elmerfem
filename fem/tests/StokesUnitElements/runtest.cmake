include(test_macros)

# One element of each 3D family -- tetrahedron, pyramid, prism, hexahedron --
# assembled by IncompressibleNSVec's equal-order pair. A smoke test that every
# family the solver can meet actually assembles; see case.sif for the two
# defects that made it worth having, and for why the MINI variant is not run
# alongside it.
RUN_ELMER_TEST()
