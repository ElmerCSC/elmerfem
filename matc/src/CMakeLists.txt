
SET(matc_SRCS c3d.c clip.c dri_ps.c eig.c error.c eval.c files.c
  funcs.c gra.c gra_com.c jacobi.c lists.c lu.c matc.c matrix.c oper.c
  optim.c parser.c str.c urand.c variable.c)

ADD_LIBRARY(matc SHARED ${matc_SRCS})

TARGET_LINK_LIBRARIES(matc m)

ADD_EXECUTABLE(Matc main.c)
TARGET_LINK_LIBRARIES(Matc matc)

SET_TARGET_PROPERTIES(Matc PROPERTIES OUTPUT_NAME matc)

INSTALL(TARGETS matc Matc RUNTIME DESTINATION "bin" LIBRARY DESTINATION "lib")
