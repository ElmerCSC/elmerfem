INCLUDE_DIRECTORIES("${CMAKE_BINARY_DIR}/fem/src")
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/fhutiter/include")
INCLUDE(AddModules)



if(INCLUDE_ELMER_MODULES)
  MAKE_SKIP_LIST(INCLUDE ${INCLUDE_ELMER_MODULES})
endif()

IF(NOT(HAVE_ADIOS2))
  LIST(APPEND SKIP_LIST ADIOS2OutputSolver)
ENDIF()

LIST(APPEND SKIP_LIST contrib)
LIST(APPEND SKIP_LIST TurbulenceModels)
ADD_ELMER_MODULES(SKIP ${SKIP_LIST})

# Turbulence model solvers live under TurbulenceModels/ for source tree
# tidiness, but each is still its own independent solver module/library.
FOREACH(_TurbMod KESolver Komega SSTKomega Spalart-Allmaras V2FSolver)
  IF(NOT INCLUDE_ELMER_MODULES OR _TurbMod IN_LIST INCLUDE_ELMER_MODULES)
    ADD_ELMER_MODULE(${_TurbMod} SOURCES TurbulenceModels/${_TurbMod}.F90)
  ENDIF()
ENDFOREACH()

IF(WITH_CONTRIB)
  ADD_SUBDIRECTORY(contrib)
ENDIF(WITH_CONTRIB)

if(TARGET EliminatePeriodic)
  ADD_DEPENDENCIES(EliminatePeriodic EliminateDirichlet)
  TARGET_LINK_LIBRARIES(EliminatePeriodic EliminateDirichlet)
endif()

if(TARGET FourierLoss)
  ADD_DEPENDENCIES(FourierLoss MagnetoDynamics)
  TARGET_LINK_LIBRARIES(FourierLoss MagnetoDynamics)
endif()


#IF(HAVE_NETCDF)
#  TARGET_LINK_LIBRARIES(SaveGridData ${NETCDF_LIBRARY})
#  TARGET_LINK_LIBRARIES(SaveGridData ${NETCDFF_LIBRARY})
#  #TODO - add to RPath
#ENDIF()
#
#IF(HAVE_HDF5)
#  TARGET_LINK_LIBRARIES(SaveGridData ${PHDF5_LIBRARY})
#  TARGET_LINK_LIBRARIES(SaveGridData ${PHDF5HL_LIBRARY})
#  #TODO - add to RPath
#ENDIF()
