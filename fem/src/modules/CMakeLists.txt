INCLUDE_DIRECTORIES("${CMAKE_BINARY_DIR}/fem/src")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/fhutiter/include")
INCLUDE(AddModules)

ADD_ELMER_MODULES(SKIP contrib)

IF(WITH_CONTRIB)
  ADD_SUBDIRECTORY(contrib)
ENDIF(WITH_CONTRIB)

ADD_DEPENDENCIES(EliminatePeriodic EliminateDirichlet)
ADD_DEPENDENCIES(FourierLoss MagnetoDynamics)
TARGET_LINK_LIBRARIES(EliminatePeriodic EliminateDirichlet)


IF(WITH_NETCDF) 
  # ---------------------- #
  # -- NetCDF libraries -- #
  #----------------------- #
  MESSAGE(STATUS "------------------------------------------------")
  MESSAGE(STATUS "Elmer package: Looking for [NetCDF] & [NetCDF Fortran] libraries")
  
  FIND_PACKAGE(NetCDF)
  
  IF(NetCDF_FOUND)
  
    SET(HAVE_NETCDF TRUE)
    MARK_AS_ADVANCED(HAVE_NETCDF)
    INCLUDE_DIRECTORIES(${NetCDF_INCLUDE_DIR})
    ADD_DEFINITIONS(-DHAVE_NETCDF)
    TARGET_LINK_LIBRARIES(SaveGridData ${NetCDF_LIBRARIES})
  
    MESSAGE(STATUS "  NetCDF:           " "${NetCDF_FOUND}")
    MESSAGE(STATUS "  NetCDF_INC:       " "${NetCDF_INCLUDE_DIR}")
    MESSAGE(STATUS "  NetCDF_LIBS:     " "${NetCDF_LIBRARIES}")
  
    # ---------------------- #
    # -- HDF5 libraries -- #
    #----------------------- #
    MESSAGE(STATUS "------------------------------------------------")
    MESSAGE(STATUS "Elmer package: Looking for [HDF5] & [HDF5 Fortran] libraries")
  
    FIND_PACKAGE(PHDF5)
  
    IF(PHDF5_FOUND)
  
      SET(HAVE_HDF5 TRUE)
      MARK_AS_ADVANCED(HAVE_HDF5)
      INCLUDE_DIRECTORIES(${PHDF5_INCLUDE_DIR})
      ADD_DEFINITIONS(-DHAVE_HDF5)
      TARGET_LINK_LIBRARIES(SaveGridData ${PHDF5_LIBRARIES})
  
      MESSAGE(STATUS "  HDF5:           " "${PHDF5_FOUND}")
      MESSAGE(STATUS "  HDF5_INC:       " "${PHDF5_INCLUDE_DIR}")
      MESSAGE(STATUS "  HDF5_LIBS:     " "${PHDF5_LIBRARIES}")
    ELSE()
      MESSAGE(STATUS "Library not found: >HDF5_FOUND< ")
      MESSAGE(WARNING "  \n Missing: >HDF5_INCLUDE_DIR< , >HDF5_LIBRARY<,  >HDF5F_LIBRARY< \n some functionalities will be disabled")
    ENDIF(PHDF5_FOUND)
  ELSE()
    MESSAGE(STATUS "Library not found: >NetCDF_FOUND< ")
    MESSAGE(WARNING "  \n Missing: >NetCDF_INCLUDE_DIR< , >NetCDF_LIBRARY<,  >NetCDFF_LIBRARY< \n some functionalities will be disabled")
  ENDIF(NetCDF_FOUND)
  
ENDIF(WITH_NETCDF) 
