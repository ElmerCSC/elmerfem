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

# message(STATUS "BLAS_LIBRARIES=${BLAS_LIBRARIES}")

find_package(BLAS REQUIRED)
if (BLAS_FOUND)
    message(STATUS "BLAS found: ${BLAS_LIBRARIES}")
else()
    message(STATUS "BLAS not found. Manually setting the path.")
    set(BLAS_LIBRARIES "/c/tools/msys64/mingw64/lib/libopenblas_64.a")
endif()


IF(UNIX)
INSTALL(TARGETS fhuti LIBRARY DESTINATION ${ELMER_INSTALL_LIB_DIR})
ENDIF()

# IF(APPLE)
#   TARGET_LINK_LIBRARIES(fhuti ${BLAS_LIBRARIES})
# ENDIF()

IF(WIN32)
  INSTALL(TARGETS fhuti RUNTIME DESTINATION "bin"
    LIBRARY DESTINATION ${ELMER_INSTALL_LIB_DIR})
#   TARGET_LINK_LIBRARIES(fhuti ${BLAS_LIBRARIES})
	# TARGET_LINK_LIBRARIES(fhuti "C:/tools/msys64/mingw64/lib/libblas.a" "C:/tools/msys64/mingw64/lib/liblapack.a")
	set_target_properties(fhuti PROPERTIES 
    	INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}"
	)
ENDIF()

# if(${CMAKE_HOST_SYSTEM_NAME} STREQUAL "MSYS")
# 	TARGET_LINK_LIBRARIES(fhuti ${BLAS_LIBRARIES})
# endif()

#INSTALL(TARGETS fhuti LIBRARY DESTINATION lib)

TARGET_LINK_LIBRARIES(fhuti ${BLAS_LIBRARIES})

RETURN()