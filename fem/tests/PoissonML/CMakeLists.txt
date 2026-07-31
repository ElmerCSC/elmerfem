INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( ptest.sif ptest.sif COPYONLY)
ADD_ELMERTEST_MODULE(PoissonML PoissonML PoissonML.f90)

# LocalMatrix/LocalMatrixBC size their element matrices with runtime-adjustable
# automatic arrays (e.g. STIFF(nd+nb,nd+nb), ~80KB for the p:9 element here).
# gfortran's default at -O2/-O3 (unlike -Ofast) heap-allocates adjustable-size
# automatic arrays above a size threshold instead of using the stack; the
# resulting per-element malloc/free churn under OpenMP was found to cause a
# rare (~0.2-1%) intermittent SIGSEGV on Windows/MinGW, reproducible only at
# full native speed (not under gdb, Valgrind, or -O0 -- all of which slow
# execution enough to avoid it). -fstack-arrays forces true stack allocation;
# verified 1500/1500 clean (1/2/4 threads x 500) with it vs ~6/1500 without.
IF(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  TARGET_COMPILE_OPTIONS(PoissonML_PoissonML PRIVATE -fstack-arrays)
ENDIF()

file(COPY PoissonML.f90 square.grd ELMERSOLVER_STARTINFO DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(PoissonML LABELS quick)
