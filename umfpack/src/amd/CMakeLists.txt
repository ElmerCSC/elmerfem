INCLUDE_DIRECTORIES(include)
SET(UMFPACK_AMD_ORIG_SOURCES
	amd_aat.c
	amd_1.c
	amd_2.c
	amd_dump.c
	amd_postorder.c
	amd_post_tree.c
	amd_defaults.c
	amd_order.c
	amd_control.c
	amd_info.c
	amd_valid.c
	amd_preprocess.c
)
SET(UMFPACK_AMD_CPP_SOURCES
	amd_i_aat.c
	amd_i_1.c
	amd_i_2.c
	amd_i_dump.c
	amd_i_postorder.c
	amd_i_post_tree.c
	amd_i_defaults.c
	amd_i_order.c
	amd_i_control.c
	amd_i_info.c
	amd_i_valid.c
	amd_i_preprocess.c
	amd_l_aat.c
	amd_l_1.c
	amd_l_2.c
	amd_l_dump.c
	amd_l_postorder.c
	amd_l_post_tree.c
	amd_l_defaults.c
	amd_l_order.c
	amd_l_control.c
	amd_l_info.c
	amd_l_valid.c
	amd_l_preprocess.c
)

SET(UMFMODIFIER_i_FLAGS -DDINT)
SET(UMFMODIFIER_l_FLAGS -DDLONG)

SET(UMFPACK_AMD_INCLUDES -I${CMAKE_CURRENT_SOURCE_DIR}/include -I${PROJECT_BINARY_DIR})
FOREACH(outfileName ${UMFPACK_AMD_CPP_SOURCES})
  STRING(REGEX REPLACE "^amd_([a-z]*)_.*" "\\1" umfModifier "${outfileName}")
  STRING(REGEX REPLACE "^amd_${umfModifier}_([a-zA-Z0-9_]*)\\.c" "\\1" stem "${outfileName}")
  SET(infileName "amd_${stem}.c")

  SET(flags ${UMFMODIFIER_${umfModifier}_FLAGS} -DNBLAS -E)

  ADD_CUSTOM_COMMAND(
    OUTPUT ${outfileName}
    COMMAND ${CMAKE_C_COMPILER}
    ARGS ${UMFPACK_AMD_INCLUDES} ${flags} ${CMAKE_CURRENT_SOURCE_DIR}/${infileName} > ${outfileName}
    DEPENDS ${infileName}
  )

ENDFOREACH(outfileName)
ADD_CUSTOM_TARGET(umfpack_amd_srcs DEPENDS ${UMFPACK_AMD_ORIG_SOURCES})

ADD_LIBRARY(amd STATIC amd_internal.h ${UMFPACK_AMD_CPP_SOURCES})
ADD_LIBRARY(amdf77 STATIC amd.f amdbar.f)

INSTALL(TARGETS amd amdf77
        DESTINATION "lib")
