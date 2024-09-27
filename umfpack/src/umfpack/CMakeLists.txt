set(UMFPACK_INCLUDE_DIRS include ../amd ../amd/include ${CMAKE_CURRENT_BINARY_DIR})

# source files which only depend on size of indices
set(UMFPACK_INT_SOURCES
  umf_analyze.c
  umf_apply_order.c
  umf_colamd.c
  umf_free.c
  umf_fsize.c
  umf_is_permutation.c
  umf_malloc.c
  umf_realloc.c
  umf_report_perm.c
  umf_singletons.c
)
# source files which depend on size of indices and complex floating point
set(UMFPACK_INT_FLOAT_SOURCES
  umf_2by2.c
  umf_assemble.c
  umf_blas3_update.c
  umf_build_tuples.c
  umf_create_element.c
  umf_dump.c
  umf_extend_front.c
  umf_garbage_collection.c
  umf_get_memory.c
  umf_grow_front.c
  umf_init_front.c
  umf_kernel.c
  umf_kernel_init.c
  umf_kernel_wrapup.c
  umf_local_search.c
  umf_lsolve.c
  umf_ltsolve.c
  umf_mem_alloc_element.c
  umf_mem_alloc_head_block.c
  umf_mem_alloc_tail_block.c
  umf_mem_free_tail_block.c
  umf_mem_init_memoryspace.c
  umf_report_vector.c
  umf_row_search.c
  umf_scale.c
  umf_scale_column.c
  umf_set_stats.c
  umf_solve.c
  umf_start_front.c
  umf_store_lu.c
  umf_symbolic_usage.c
  umf_transpose.c
  umf_triplet.c
  umf_tuple_lengths.c
  umf_usolve.c
  umf_utsolve.c
  umf_valid_numeric.c
  umf_valid_symbolic.c
  umfpack_col_to_triplet.c
  umfpack_defaults.c
  umfpack_free_numeric.c
  umfpack_free_symbolic.c
  umfpack_get_determinant.c
  umfpack_get_lunz.c
  umfpack_get_numeric.c
  umfpack_get_symbolic.c
  umfpack_load_numeric.c
  umfpack_load_symbolic.c
  umfpack_numeric.c
  umfpack_qsymbolic.c
  umfpack_report_control.c
  umfpack_report_info.c
  umfpack_report_matrix.c
  umfpack_report_numeric.c
  umfpack_report_perm.c
  umfpack_report_status.c
  umfpack_report_symbolic.c
  umfpack_report_triplet.c
  umfpack_report_vector.c
  umfpack_save_numeric.c
  umfpack_save_symbolic.c
  umfpack_scale.c
  umfpack_solve.c
  umfpack_symbolic.c
  umfpack_transpose.c
  umfpack_triplet_to_col.c
)

function(add_object_library target_name)
  # Add an object library with specific preprocessor flags
  cmake_parse_arguments(OBJ "" "" "SOURCES;FLAGS" ${ARGN})
  add_library(${target_name} OBJECT ${OBJ_SOURCES})
  target_include_directories(${target_name} PRIVATE ${UMFPACK_INCLUDE_DIRS})
  target_compile_definitions(${target_name} PRIVATE NBLAS ${OBJ_FLAGS})
endfunction()

# add object libraries with the required preprocessor flags
add_object_library(umfpack_i SOURCES ${UMFPACK_INT_SOURCES} FLAGS DINT)
add_object_library(umfpack_l SOURCES ${UMFPACK_INT_SOURCES} FLAGS DLONG)
add_object_library(umfpack_di SOURCES ${UMFPACK_INT_FLOAT_SOURCES} FLAGS DINT)
add_object_library(umfpack_dl SOURCES ${UMFPACK_INT_FLOAT_SOURCES} FLAGS DLONG)
add_object_library(umfpack_zi SOURCES ${UMFPACK_INT_FLOAT_SOURCES} FLAGS ZINT)
add_object_library(umfpack_zl SOURCES ${UMFPACK_INT_FLOAT_SOURCES} FLAGS ZLONG)

# add targets for objects that require a special set of flags
add_object_library(hsolve_di SOURCES umf_ltsolve.c umf_utsolve.c FLAGS DINT CONJUGATE_SOLVE)
add_object_library(hsolve_dl SOURCES umf_ltsolve.c umf_utsolve.c FLAGS DLONG CONJUGATE_SOLVE)
add_object_library(hsolve_zi SOURCES umf_ltsolve.c umf_utsolve.c FLAGS ZINT CONJUGATE_SOLVE)
add_object_library(hsolve_zl SOURCES umf_ltsolve.c umf_utsolve.c FLAGS ZLONG CONJUGATE_SOLVE)

add_object_library(triplet_map_nox_di SOURCES umf_triplet.c FLAGS DINT MAP)
add_object_library(triplet_map_nox_dl SOURCES umf_triplet.c FLAGS DLONG MAP)
add_object_library(triplet_map_nox_zi SOURCES umf_triplet.c FLAGS ZINT MAP)
add_object_library(triplet_map_nox_zl SOURCES umf_triplet.c FLAGS ZLONG MAP)

add_object_library(triplet_nomap_x_di SOURCES umf_triplet.c FLAGS DINT DO_VALUES)
add_object_library(triplet_nomap_x_dl SOURCES umf_triplet.c FLAGS DLONG DO_VALUES)
add_object_library(triplet_nomap_x_zi SOURCES umf_triplet.c FLAGS ZINT DO_VALUES)
add_object_library(triplet_nomap_x_zl SOURCES umf_triplet.c FLAGS ZLONG DO_VALUES)

add_object_library(triplet_map_x_di SOURCES umf_triplet.c FLAGS DINT MAP DO_VALUES)
add_object_library(triplet_map_x_dl SOURCES umf_triplet.c FLAGS DLONG MAP DO_VALUES)
add_object_library(triplet_map_x_zi SOURCES umf_triplet.c FLAGS ZINT MAP DO_VALUES)
add_object_library(triplet_map_x_zl SOURCES umf_triplet.c FLAGS ZLONG MAP DO_VALUES)

add_object_library(assemble_fixq_di SOURCES umf_assemble.c FLAGS DINT FIXQ)
add_object_library(assemble_fixq_dl SOURCES umf_assemble.c FLAGS DLONG FIXQ)
add_object_library(assemble_fixq_zi SOURCES umf_assemble.c FLAGS ZINT FIXQ)
add_object_library(assemble_fixq_zl SOURCES umf_assemble.c FLAGS ZLONG FIXQ)

add_object_library(store_lu_drop_di SOURCES umf_store_lu.c FLAGS DINT DROP)
add_object_library(store_lu_drop_dl SOURCES umf_store_lu.c FLAGS DLONG DROP)
add_object_library(store_lu_drop_zi SOURCES umf_store_lu.c FLAGS ZINT DROP)
add_object_library(store_lu_drop_zl SOURCES umf_store_lu.c FLAGS ZLONG DROP)

add_object_library(wsolve_di SOURCES umfpack_solve.c FLAGS DINT WSOLVE)
add_object_library(wsolve_dl SOURCES umfpack_solve.c FLAGS DLONG WSOLVE)
add_object_library(wsolve_zi SOURCES umfpack_solve.c FLAGS ZINT WSOLVE)
add_object_library(wsolve_zl SOURCES umfpack_solve.c FLAGS ZLONG WSOLVE)


add_library(umfpack STATIC
  umfpack_timer.c umfpack_tictoc.c
  $<TARGET_OBJECTS:umfpack_i> $<TARGET_OBJECTS:umfpack_l>
  $<TARGET_OBJECTS:umfpack_di> $<TARGET_OBJECTS:umfpack_dl>
  $<TARGET_OBJECTS:umfpack_zi> $<TARGET_OBJECTS:umfpack_zl>
  $<TARGET_OBJECTS:hsolve_di> $<TARGET_OBJECTS:hsolve_dl>
  $<TARGET_OBJECTS:hsolve_zi> $<TARGET_OBJECTS:hsolve_zi>
  $<TARGET_OBJECTS:triplet_map_nox_di> $<TARGET_OBJECTS:triplet_map_nox_dl>
  $<TARGET_OBJECTS:triplet_map_nox_zi> $<TARGET_OBJECTS:triplet_map_nox_zl>
  $<TARGET_OBJECTS:triplet_nomap_x_di> $<TARGET_OBJECTS:triplet_nomap_x_dl>
  $<TARGET_OBJECTS:triplet_nomap_x_zi> $<TARGET_OBJECTS:triplet_nomap_x_zl>
  $<TARGET_OBJECTS:triplet_map_x_di> $<TARGET_OBJECTS:triplet_map_x_dl>
  $<TARGET_OBJECTS:triplet_map_x_zi> $<TARGET_OBJECTS:triplet_map_x_zl>
  $<TARGET_OBJECTS:assemble_fixq_di> $<TARGET_OBJECTS:assemble_fixq_dl>
  $<TARGET_OBJECTS:assemble_fixq_zi> $<TARGET_OBJECTS:assemble_fixq_zl>
  $<TARGET_OBJECTS:store_lu_drop_di> $<TARGET_OBJECTS:store_lu_drop_dl>
  $<TARGET_OBJECTS:store_lu_drop_zi> $<TARGET_OBJECTS:store_lu_drop_zl>
  $<TARGET_OBJECTS:wsolve_di> $<TARGET_OBJECTS:wsolve_dl>
  $<TARGET_OBJECTS:wsolve_zi> $<TARGET_OBJECTS:wsolve_zl>)

target_include_directories(umfpack PRIVATE ${UMFPACK_INCLUDE_DIRS})
target_compile_definitions(umfpack PRIVATE NBLAS)


install(TARGETS umfpack
  ARCHIVE DESTINATION ${ELMER_INSTALL_LIB_DIR}
  RUNTIME DESTINATION ${ELMER_INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${ELMER_INSTALL_LIB_DIR})
