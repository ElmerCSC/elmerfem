INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/../test_macros.cmake)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE( interpolationtest.sif interpolationtest.sif COPYONLY)
ADD_ELMERTEST_MODULE(EdgeElementInterpolation EdgeElementSolve EdgeElementSolve.f90)

file(COPY ELMERSOLVER_STARTINFO toroidalcoil TEST.RESULT EdgeElementSolve.f90 DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(EdgeElementInterpolation)
