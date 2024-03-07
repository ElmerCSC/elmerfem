INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(IBC_circuit.sif IBC_circuit.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO wire.grd circuit.definition DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(mgdyn_harmonic_wire_impedanceBC_circuit LABELS whitney)
