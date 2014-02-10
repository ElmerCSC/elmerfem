INCLUDE_DIRECTORIES(include ../amd ../amd/include ${CMAKE_CURRENT_BINARY_DIR})

SET(UMFPACK_ORIG_SOURCES
	umf_usolve.c
	umf_triplet.c
	umfpack_get_determinant.c
	umfpack_get_lunz.c
	umfpack_get_symbolic.c
	umfpack_report_numeric.c
	umfpack_qsymbolic.c
	umf_assemble.c
	umf_local_search.c
	umfpack_report_matrix.c
	umfpack_defaults.c
	umf_scale.c
	umf_kernel_init.c
	umfpack_solve.c
	umf_start_front.c
	umfpack_save_numeric.c
	umfpack_numeric.c
	umf_2by2.c
	umfpack_col_to_triplet.c
	umf_free.c
	umf_lsolve.c
	umfpack_save_symbolic.c
	umf_grow_front.c
	umfpack_free_numeric.c
	umf_solve.c
	umf_ltsolve.c
	umf_kernel.c
	umfpack_symbolic.c
	umf_row_search.c
	umf_mem_alloc_tail_block.c
	umf_garbage_collection.c
	umf_colamd.c
	umf_fsize.c
	umf_is_permutation.c
	umf_get_memory.c
	umf_mem_alloc_element.c
	umf_mem_free_tail_block.c
	umfpack_report_status.c
	umfpack_report_control.c
	umf_utsolve.c
	umf_symbolic_usage.c
	umf_valid_symbolic.c
	umfpack_report_perm.c
	umfpack_report_triplet.c
	umfpack_report_symbolic.c
	umfpack_report_vector.c
	umf_dump.c
	umf_blas3_update.c
	umfpack_triplet_to_col.c
	umfpack_get_numeric.c
	umf_transpose.c
	umf_kernel_wrapup.c
	umfpack_free_symbolic.c
	umf_tuple_lengths.c
	umf_report_perm.c
	umf_realloc.c
	umf_mem_alloc_head_block.c
	umfpack_load_numeric.c
	umf_scale_column.c
	umf_create_element.c
	umf_apply_order.c
	umfpack_scale.c
	umf_set_stats.c
	umf_extend_front.c
	umf_mem_init_memoryspace.c
	umf_valid_numeric.c
	umf_init_front.c
	umfpack_transpose.c
	umfpack_report_info.c
	umf_analyze.c
	umfpack_load_symbolic.c
	umf_malloc.c
	umf_store_lu.c
	umf_singletons.c
	umf_report_vector.c
	umf_build_tuples.c
)

SET(UMFPACK_CPP_SOURCES
	umf_i_analyze.c
	umf_i_apply_order.c
	umf_i_colamd.c
	umf_i_free.c
	umf_i_fsize.c
	umf_i_is_permutation.c
	umf_i_malloc.c
	umf_i_realloc.c
	umf_i_report_perm.c
	umf_i_singletons.c
	umf_l_analyze.c
	umf_l_apply_order.c
	umf_l_colamd.c
	umf_l_free.c
	umf_l_fsize.c
	umf_l_is_permutation.c
	umf_l_malloc.c
	umf_l_realloc.c
	umf_l_report_perm.c
	umf_l_singletons.c
	umf_di_lhsolve.c
	umf_di_uhsolve.c
	umf_di_triplet_map_nox.c
	umf_di_triplet_nomap_x.c
	umf_di_triplet_nomap_nox.c
	umf_di_triplet_map_x.c
	umf_di_assemble_fixq.c
	umf_di_store_lu_drop.c
	umf_di_assemble.c
	umf_di_blas3_update.c
	umf_di_build_tuples.c
	umf_di_create_element.c
	umf_di_dump.c
	umf_di_extend_front.c
	umf_di_garbage_collection.c
	umf_di_get_memory.c
	umf_di_init_front.c
	umf_di_kernel.c
	umf_di_kernel_init.c
	umf_di_kernel_wrapup.c
	umf_di_local_search.c
	umf_di_lsolve.c
	umf_di_ltsolve.c
	umf_di_mem_alloc_element.c
	umf_di_mem_alloc_head_block.c
	umf_di_mem_alloc_tail_block.c
	umf_di_mem_free_tail_block.c
	umf_di_mem_init_memoryspace.c
	umf_di_report_vector.c
	umf_di_row_search.c
	umf_di_scale_column.c
	umf_di_set_stats.c
	umf_di_solve.c
	umf_di_symbolic_usage.c
	umf_di_transpose.c
	umf_di_tuple_lengths.c
	umf_di_usolve.c
	umf_di_utsolve.c
	umf_di_valid_numeric.c
	umf_di_valid_symbolic.c
	umf_di_grow_front.c
	umf_di_start_front.c
	umf_di_2by2.c
	umf_di_store_lu.c
	umf_di_scale.c
	umfpack_di_wsolve.c
	umfpack_di_col_to_triplet.c
	umfpack_di_defaults.c
	umfpack_di_free_numeric.c
	umfpack_di_free_symbolic.c
	umfpack_di_get_numeric.c
	umfpack_di_get_lunz.c
	umfpack_di_get_symbolic.c
	umfpack_di_get_determinant.c
	umfpack_di_numeric.c
	umfpack_di_qsymbolic.c
	umfpack_di_report_control.c
	umfpack_di_report_info.c
	umfpack_di_report_matrix.c
	umfpack_di_report_numeric.c
	umfpack_di_report_perm.c
	umfpack_di_report_status.c
	umfpack_di_report_symbolic.c
	umfpack_di_report_triplet.c
	umfpack_di_report_vector.c
	umfpack_di_solve.c
	umfpack_di_symbolic.c
	umfpack_di_transpose.c
	umfpack_di_triplet_to_col.c
	umfpack_di_scale.c
	umfpack_di_load_numeric.c
	umfpack_di_save_numeric.c
	umfpack_di_load_symbolic.c
	umfpack_di_save_symbolic.c
	umf_dl_lhsolve.c
	umf_dl_uhsolve.c
	umf_dl_triplet_map_nox.c
	umf_dl_triplet_nomap_x.c
	umf_dl_triplet_nomap_nox.c
	umf_dl_triplet_map_x.c
	umf_dl_assemble_fixq.c
	umf_dl_store_lu_drop.c
	umf_dl_assemble.c
	umf_dl_blas3_update.c
	umf_dl_build_tuples.c
	umf_dl_create_element.c
	umf_dl_dump.c
	umf_dl_extend_front.c
	umf_dl_garbage_collection.c
	umf_dl_get_memory.c
	umf_dl_init_front.c
	umf_dl_kernel.c
	umf_dl_kernel_init.c
	umf_dl_kernel_wrapup.c
	umf_dl_local_search.c
	umf_dl_lsolve.c
	umf_dl_ltsolve.c
	umf_dl_mem_alloc_element.c
	umf_dl_mem_alloc_head_block.c
	umf_dl_mem_alloc_tail_block.c
	umf_dl_mem_free_tail_block.c
	umf_dl_mem_init_memoryspace.c
	umf_dl_report_vector.c
	umf_dl_row_search.c
	umf_dl_scale_column.c
	umf_dl_set_stats.c
	umf_dl_solve.c
	umf_dl_symbolic_usage.c
	umf_dl_transpose.c
	umf_dl_tuple_lengths.c
	umf_dl_usolve.c
	umf_dl_utsolve.c
	umf_dl_valid_numeric.c
	umf_dl_valid_symbolic.c
	umf_dl_grow_front.c
	umf_dl_start_front.c
	umf_dl_2by2.c
	umf_dl_store_lu.c
	umf_dl_scale.c
	umfpack_dl_wsolve.c
	umfpack_dl_col_to_triplet.c
	umfpack_dl_defaults.c
	umfpack_dl_free_numeric.c
	umfpack_dl_free_symbolic.c
	umfpack_dl_get_numeric.c
	umfpack_dl_get_lunz.c
	umfpack_dl_get_symbolic.c
	umfpack_dl_get_determinant.c
	umfpack_dl_numeric.c
	umfpack_dl_qsymbolic.c
	umfpack_dl_report_control.c
	umfpack_dl_report_info.c
	umfpack_dl_report_matrix.c
	umfpack_dl_report_numeric.c
	umfpack_dl_report_perm.c
	umfpack_dl_report_status.c
	umfpack_dl_report_symbolic.c
	umfpack_dl_report_triplet.c
	umfpack_dl_report_vector.c
	umfpack_dl_solve.c
	umfpack_dl_symbolic.c
	umfpack_dl_transpose.c
	umfpack_dl_triplet_to_col.c
	umfpack_dl_scale.c
	umfpack_dl_load_numeric.c
	umfpack_dl_save_numeric.c
	umfpack_dl_load_symbolic.c
	umfpack_dl_save_symbolic.c
	umf_zi_lhsolve.c
	umf_zi_uhsolve.c
	umf_zi_triplet_map_nox.c
	umf_zi_triplet_nomap_x.c
	umf_zi_triplet_nomap_nox.c
	umf_zi_triplet_map_x.c
	umf_zi_assemble_fixq.c
	umf_zi_store_lu_drop.c
	umf_zi_assemble.c
	umf_zi_blas3_update.c
	umf_zi_build_tuples.c
	umf_zi_create_element.c
	umf_zi_dump.c
	umf_zi_extend_front.c
	umf_zi_garbage_collection.c
	umf_zi_get_memory.c
	umf_zi_init_front.c
	umf_zi_kernel.c
	umf_zi_kernel_init.c
	umf_zi_kernel_wrapup.c
	umf_zi_local_search.c
	umf_zi_lsolve.c
	umf_zi_ltsolve.c
	umf_zi_mem_alloc_element.c
	umf_zi_mem_alloc_head_block.c
	umf_zi_mem_alloc_tail_block.c
	umf_zi_mem_free_tail_block.c
	umf_zi_mem_init_memoryspace.c
	umf_zi_report_vector.c
	umf_zi_row_search.c
	umf_zi_scale_column.c
	umf_zi_set_stats.c
	umf_zi_solve.c
	umf_zi_symbolic_usage.c
	umf_zi_transpose.c
	umf_zi_tuple_lengths.c
	umf_zi_usolve.c
	umf_zi_utsolve.c
	umf_zi_valid_numeric.c
	umf_zi_valid_symbolic.c
	umf_zi_grow_front.c
	umf_zi_start_front.c
	umf_zi_2by2.c
	umf_zi_store_lu.c
	umf_zi_scale.c
	umfpack_zi_wsolve.c
	umfpack_zi_col_to_triplet.c
	umfpack_zi_defaults.c
	umfpack_zi_free_numeric.c
	umfpack_zi_free_symbolic.c
	umfpack_zi_get_numeric.c
	umfpack_zi_get_lunz.c
	umfpack_zi_get_symbolic.c
	umfpack_zi_get_determinant.c
	umfpack_zi_numeric.c
	umfpack_zi_qsymbolic.c
	umfpack_zi_report_control.c
	umfpack_zi_report_info.c
	umfpack_zi_report_matrix.c
	umfpack_zi_report_numeric.c
	umfpack_zi_report_perm.c
	umfpack_zi_report_status.c
	umfpack_zi_report_symbolic.c
	umfpack_zi_report_triplet.c
	umfpack_zi_report_vector.c
	umfpack_zi_solve.c
	umfpack_zi_symbolic.c
	umfpack_zi_transpose.c
	umfpack_zi_triplet_to_col.c
	umfpack_zi_scale.c
	umfpack_zi_load_numeric.c
	umfpack_zi_save_numeric.c
	umfpack_zi_load_symbolic.c
	umfpack_zi_save_symbolic.c
	umf_zl_lhsolve.c
	umf_zl_uhsolve.c
	umf_zl_triplet_map_nox.c
	umf_zl_triplet_nomap_x.c
	umf_zl_triplet_nomap_nox.c
	umf_zl_triplet_map_x.c
	umf_zl_assemble_fixq.c
	umf_zl_store_lu_drop.c
	umf_zl_assemble.c
	umf_zl_blas3_update.c
	umf_zl_build_tuples.c
	umf_zl_create_element.c
	umf_zl_dump.c
	umf_zl_extend_front.c
	umf_zl_garbage_collection.c
	umf_zl_get_memory.c
	umf_zl_init_front.c
	umf_zl_kernel.c
	umf_zl_kernel_init.c
	umf_zl_kernel_wrapup.c
	umf_zl_local_search.c
	umf_zl_lsolve.c
	umf_zl_ltsolve.c
	umf_zl_mem_alloc_element.c
	umf_zl_mem_alloc_head_block.c
	umf_zl_mem_alloc_tail_block.c
	umf_zl_mem_free_tail_block.c
	umf_zl_mem_init_memoryspace.c
	umf_zl_report_vector.c
	umf_zl_row_search.c
	umf_zl_scale_column.c
	umf_zl_set_stats.c
	umf_zl_solve.c
	umf_zl_symbolic_usage.c
	umf_zl_transpose.c
	umf_zl_tuple_lengths.c
	umf_zl_usolve.c
	umf_zl_utsolve.c
	umf_zl_valid_numeric.c
	umf_zl_valid_symbolic.c
	umf_zl_grow_front.c
	umf_zl_start_front.c
	umf_zl_2by2.c
	umf_zl_store_lu.c
	umf_zl_scale.c
	umfpack_zl_wsolve.c
	umfpack_zl_col_to_triplet.c
	umfpack_zl_defaults.c
	umfpack_zl_free_numeric.c
	umfpack_zl_free_symbolic.c
	umfpack_zl_get_numeric.c
	umfpack_zl_get_lunz.c
	umfpack_zl_get_symbolic.c
	umfpack_zl_get_determinant.c
	umfpack_zl_numeric.c
	umfpack_zl_qsymbolic.c
	umfpack_zl_report_control.c
	umfpack_zl_report_info.c
	umfpack_zl_report_matrix.c
	umfpack_zl_report_numeric.c
	umfpack_zl_report_perm.c
	umfpack_zl_report_status.c
	umfpack_zl_report_symbolic.c
	umfpack_zl_report_triplet.c
	umfpack_zl_report_vector.c
	umfpack_zl_solve.c
	umfpack_zl_symbolic.c
	umfpack_zl_transpose.c
	umfpack_zl_triplet_to_col.c
	umfpack_zl_scale.c
	umfpack_zl_load_numeric.c
	umfpack_zl_save_numeric.c
	umfpack_zl_load_symbolic.c
	umfpack_zl_save_symbolic.c
)

SET(UMFMODIFIER_i_FLAGS -DDINT)
SET(UMFMODIFIER_l_FLAGS -DDLONG)
SET(UMFMODIFIER_di_FLAGS -DDINT)
SET(UMFMODIFIER_dl_FLAGS -DDLONG)
SET(UMFMODIFIER_zi_FLAGS -DZINT)
SET(UMFMODIFIER_zl_FLAGS -DZLONG)

SET(UMFPACK_ODDITY_lhsolve ltsolve)
SET(UMFPACK_ODDITY_lhsolve_FLAGS -DCONJUGATE_SOLVE)
SET(UMFPACK_ODDITY_uhsolve utsolve)
SET(UMFPACK_ODDITY_uhsolve_FLAGS -DCONJUGATE_SOLVE)
SET(UMFPACK_ODDITY_triplet_map_nox triplet)
SET(UMFPACK_ODDITY_triplet_map_nox_FLAGS -DMAP)
SET(UMFPACK_ODDITY_triplet_nomap_x triplet)
SET(UMFPACK_ODDITY_triplet_nomap_x_FLAGS -DDO_VALUES)
SET(UMFPACK_ODDITY_triplet_nomap_nox triplet)
SET(UMFPACK_ODDITY_triplet_nomap_nox_FLAGS)
SET(UMFPACK_ODDITY_triplet_map_x triplet)
SET(UMFPACK_ODDITY_triplet_map_x_FLAGS -DMAP -DDO_VALUES)
SET(UMFPACK_ODDITY_assemble_fixq assemble)
SET(UMFPACK_ODDITY_assemble_fixq_FLAGS -DFIXQ)
SET(UMFPACK_ODDITY_store_lu_drop store_lu)
SET(UMFPACK_ODDITY_store_lu_drop_FLAGS -DDROP)
SET(UMFPACK_ODDITY_wsolve solve)
SET(UMFPACK_ODDITY_wsolve_FLAGS -DWSOLVE)

SET(UMFPACK_INCLUDES -I${CMAKE_CURRENT_SOURCE_DIR}/include 
  -I${CMAKE_CURRENT_SOURCE_DIR}/../amd 
  -I${CMAKE_CURRENT_SOURCE_DIR}/../amd/include 
  -I${PROJECT_BINARY_DIR} 
  -I${CMAKE_CURRENT_BINARY_DIR})

FOREACH(outfileName ${UMFPACK_CPP_SOURCES})
  STRING(REGEX REPLACE "^(umf|umfpack)_([a-z]*)_.*" "\\1" prefix "${outfileName}")
  STRING(REGEX REPLACE "^(umf|umfpack)_([a-z]*)_.*" "\\2" umfModifier "${outfileName}")
  STRING(REGEX REPLACE "^(umf|umfpack)_${umfModifier}_([a-zA-Z0-9_]*)\\.c" "\\2" stem "${outfileName}")
  SET(infileName "${prefix}_${stem}.c")

  SET(flags ${UMFMODIFIER_${umfModifier}_FLAGS} -DNBLAS -E)
  IF(DEFINED UMFPACK_ODDITY_${stem}_FLAGS)
    SET(flags ${flags} ${UMFPACK_ODDITY_${stem}_FLAGS})
#    MESSAGE(${flags})
  ENDIF(DEFINED UMFPACK_ODDITY_${stem}_FLAGS)

  IF(DEFINED UMFPACK_ODDITY_${stem})
    SET(stem ${UMFPACK_ODDITY_${stem}})
    SET(infileName "${prefix}_${stem}.c")
  ENDIF(DEFINED UMFPACK_ODDITY_${stem})

  ADD_CUSTOM_COMMAND(
    OUTPUT ${outfileName}
    COMMAND ${CMAKE_C_COMPILER}
    ARGS ${UMFPACK_INCLUDES} ${flags} ${CMAKE_CURRENT_SOURCE_DIR}/${infileName} > ${outfileName}
    DEPENDS ${infileName}
  )

ENDFOREACH(outfileName)
ADD_CUSTOM_TARGET(umfpack_srcs DEPENDS ${UMFPACK_ORIG_SOURCES})

ADD_LIBRARY(umfpack STATIC umfpack_timer.c umfpack_tictoc.c ${UMFPACK_CPP_SOURCES})
ADD_DEPENDENCIES(umfpack umfpack_srcs)

INSTALL(TARGETS umfpack
        DESTINATION "lib")
