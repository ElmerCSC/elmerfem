MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "Elmer/Ice package: Compiling PDAF coupling solver")



SET(CMAKE_Fortran_MODULE_DIRECTORY
  ${PROJECT_BINARY_DIR}/fmodules CACHE PATH "Directory for Fortran modules")

# ---------------------------------------------- #
# ---------------------------------------------- #
SET(ElmerIce_PDAF_SRC
	PDAF_COUPLER.F90
	mod_io.F90
	mod_assimilation.F90
	mod_observations.F90
	mod_statevector_pdaf.F90
	next_observation_pdaf.F90
	init_ens_pdaf.F90
	init_local_domains.F90
	mod_synobs.F90
	get_obs_f_pdaf.F90
	callback_obs_pdafomi.F90
	obs_AtNodes_pdafomi.F90
	prepoststep_pdaf.F90)

ADD_LIBRARY(ElmerIce_PDAFCoupler SHARED ${ElmerIce_PDAF_SRC})

INCLUDE_DIRECTORIES(${NETCDF_INCLUDE_DIR} ${PDAF_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(ElmerIce_PDAFCoupler Elmer::MPI_Fortran ${NETCDF_LIBRARIES} ${PDAF_LIBRARIES})

# Library object
SET_TARGET_PROPERTIES(ElmerIce_PDAFCoupler PROPERTIES PREFIX "")
SET_TARGET_PROPERTIES(ElmerIce_PDAFCoupler PROPERTIES LINKER_LANGUAGE Fortran)
IF(NOT(WIN32))
    SET_TARGET_PROPERTIES(ElmerIce_PDAFCoupler PROPERTIES 
                          INSTALL_RPATH ${ELMERSOLVER_RPATH_STRING})
ENDIF()

INSTALL(TARGETS ElmerIce_PDAFCoupler
  RUNTIME DESTINATION "share/elmersolver/lib"
  LIBRARY DESTINATION "share/elmersolver/lib"
  ARCHIVE DESTINATION "share/elmersolver/lib")
