
SET(matc_SRCS c3d.c clip.c dri_ps.c eig.c error.c eval.c files.c
  funcs.c gra.c gra_com.c jacobi.c lists.c lu.c matc.c matrix.c oper.c
  optim.c parser.c str.c urand.c variable.c)

ADD_LIBRARY(matc SHARED ${matc_SRCS})

TARGET_LINK_LIBRARIES(matc m)

ADD_EXECUTABLE(Matc_bin main.c)
TARGET_LINK_LIBRARIES(Matc_bin matc)

SET_TARGET_PROPERTIES(Matc_bin PROPERTIES OUTPUT_NAME matc)
IF(NOT(WIN32))
  SET_TARGET_PROPERTIES(Matc_bin PROPERTIES INSTALL_RPATH ${ELMERSOLVER_RPATH_STRING})
ENDIF()

INSTALL(TARGETS matc Matc_bin RUNTIME DESTINATION "bin" LIBRARY DESTINATION ${ELMER_INSTALL_LIB_DIR})
