
ADD_SUBDIRECTORY(metis)

SET(elmergrid_SRCS common.h femdef.h femelmer.h femfilein.h
  femfileout.h feminfo.h femknot.h femmesh.h femsolve.h femtools.h
  femtypes.h nrutil.h common.c femelmer.c femfilein.c femfileout.c
  feminfo.c femknot.c femmesh.c fempre.c femsolve.c femtools.c
  nrutil.c)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/elmergrid/src/metis)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/elmergrid/src)

ADD_EXECUTABLE(ElmerGrid ${elmergrid_SRCS})

TARGET_LINK_LIBRARIES(ElmerGrid matc metis)

INSTALL(PROGRAMS ${CMAKE_BINARY_DIR}/elmergrid/src/ElmerGrid DESTINATION "bin")
