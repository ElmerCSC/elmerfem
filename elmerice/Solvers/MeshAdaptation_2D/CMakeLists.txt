
SET(CMAKE_Fortran_MODULE_DIRECTORY
  ${PROJECT_BINARY_DIR}/fmodules CACHE PATH "Directory for Fortran modules")

# ---------------------------------------------- #
# -- set Solvers without external dependency  -- #
# ---------------------------------------------- #
SET(ElmerIce_MeshAdapt2D_SRC
    Compute2DNodalGradient.F90
    MMG2D_MetricAniso.F90
    MMG2D_MetricIntersect.F90)

# ---------------------- #
# -- LOOK FOR MMG     -- #
#----------------------- #
MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "Elmer/Ice package: Mesh adaptation 2D looking for [Mmg] tools ")

FIND_PACKAGE(MMG)

IF(MMG_FOUND)
  SET(HAVE_MMG TRUE CACHE BOOL "Has Mmg tools for Elmer/Ice Mesh adaptation" ) 
  MARK_AS_ADVANCED(HAVE_MMG)
  MESSAGE(STATUS "  Mmg:           " "${MMG_FOUND}")
  MESSAGE(STATUS "  Mmg_INC:       " "${MMG_INCLUDE_DIR}")
  MESSAGE(STATUS "  Mmg_LIB:      " "${MMG_LIBRARY}")
  MESSAGE(STATUS "Compile MMG2DSolver")
  LIST(APPEND  ElmerIce_MeshAdapt2D_SRC MMG2DSolver.F90)
ELSE()
  MESSAGE(STATUS "  Library not found: >MMG_FOUND< ")
  MESSAGE(STATUS "    Missing: >MMG_INCLUDE_DIR< , >MMG_LIBRARY<, to compile MMG2DSolver")
ENDIF(MMG_FOUND)

MESSAGE(STATUS "------------------------------------------------")

ADD_LIBRARY(ElmerIce_MeshAdapt2D SHARED ${ElmerIce_MeshAdapt2D_SRC})


# Library object
GET_TARGET_PROPERTY(CURRENT_LINK_FLAGS ElmerIce_MeshAdapt2D LINK_FLAGS) 
GET_TARGET_PROPERTY(CURRENT_COMPILE_FLAGS ElmerIce_MeshAdapt2D COMPILE_FLAGS)
IF(CURRENT_LINK_FLAGS)
  SET(CURR_LFLAGS "${CURRENT_LINK_FLAGS};${MPI_Fortran_LINK_FLAGS}")
ELSE()
  SET(CURR_LFLAGS "${MPI_Fortran_LINK_FLAGS}")
ENDIF()
IF(CURRENT_COMPILE_FLAGS)
  SET(CURR_CFLAGS "${CURRENT_COMPILE_FLAGS};${MPI_Fortran_COMPILE_FLAGS}")
ELSE()
  SET(CURR_CFLAGS "${MPI_Fortran_COMPILE_FLAGS}")
ENDIF()

SET_TARGET_PROPERTIES(ElmerIce_MeshAdapt2D PROPERTIES PREFIX "")  
SET_TARGET_PROPERTIES(ElmerIce_MeshAdapt2D PROPERTIES LINK_FLAGS "${CURR_LFLAGS}")
SET_TARGET_PROPERTIES(ElmerIce_MeshAdapt2D PROPERTIES COMPILE_FLAGS "${CURR_CFLAGS}")
SET_TARGET_PROPERTIES(ElmerIce_MeshAdapt2D PROPERTIES LINKER_LANGUAGE Fortran)
IF(NOT(WIN32))
    SET_TARGET_PROPERTIES(ElmerIce_MeshAdapt2D PROPERTIES 
                          INSTALL_RPATH ${ELMERSOLVER_RPATH_STRING})
ENDIF()

IF(MMG_FOUND)
  INCLUDE_DIRECTORIES(${MMG_INCLUDE_DIR}/mmg)
  TARGET_LINK_LIBRARIES(ElmerIce_MeshAdapt2D ${MMG_LIBRARY} elmersolver)
ELSE()
  TARGET_LINK_LIBRARIES(ElmerIce_MeshAdapt2D elmersolver)
ENDIF()

INSTALL(TARGETS ElmerIce_MeshAdapt2D LIBRARY DESTINATION "share/elmersolver/lib"
	    RUNTIME DESTINATION "share/elmersolver/lib")

