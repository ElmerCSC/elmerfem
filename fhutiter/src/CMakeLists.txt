set(LIBSRC
    huti_aux.F90
    huti_bicgstab_2.F90
    huti_bicgstab.F90
    huti_cg.F90
    huti_cgs.F90
    huti_fdefs.h
    huti_gmres.F90
    huti_intdefs.h
    huti_qmr.F90
    huti_tfqmr.F90
    huti_interfaces.F90
    huti_sfe.F90)

add_library(fhuti SHARED ${LIBSRC})

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

if(UNIX)
  install(TARGETS fhuti LIBRARY DESTINATION ${ELMER_INSTALL_LIB_DIR})
endif()

if(APPLE)
  target_link_libraries(fhuti ${BLAS_LIBRARIES})
endif()

if(WIN32)
  install(
    TARGETS fhuti
    RUNTIME DESTINATION "bin"
    LIBRARY DESTINATION ${ELMER_INSTALL_LIB_DIR})
  target_link_libraries(fhuti ${BLAS_LIBRARIES})
endif()

# INSTALL(TARGETS fhuti LIBRARY DESTINATION lib)
