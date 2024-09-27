set(UMFPACK_AMD_SOURCES
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

add_library(amd_i OBJECT ${UMFPACK_AMD_SOURCES})
target_include_directories(amd_i PRIVATE include)
target_compile_definitions(amd_i PRIVATE NBLAS DINT)

add_library(amd_l OBJECT ${UMFPACK_AMD_SOURCES})
target_include_directories(amd_l PRIVATE include)
target_compile_definitions(amd_l PRIVATE NBLAS DLONG)

add_library(amd STATIC $<TARGET_OBJECTS:amd_i> $<TARGET_OBJECTS:amd_l>)
add_library(amdf77 STATIC amd.f amdbar.f)

install(TARGETS amd amdf77
  ARCHIVE DESTINATION ${ELMER_INSTALL_LIB_DIR}
  RUNTIME DESTINATION ${ELMER_INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${ELMER_INSTALL_LIB_DIR})
