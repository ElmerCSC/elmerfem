SET(LIBSRC 
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

ADD_LIBRARY(fhuti SHARED ${LIBSRC})

IF(APPLE OR WIN32)
  TARGET_LINK_LIBRARIES(fhuti ${BLAS_LIBRARIES})
ENDIF()

INSTALL(TARGETS fhuti
  RUNTIME DESTINATION "bin"
  LIBRARY DESTINATION ${ELMER_INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${ELMER_INSTALL_LIB_DIR})
