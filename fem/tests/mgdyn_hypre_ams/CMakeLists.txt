if(WITH_Hypre)
  include(test_macros)

  configure_file(case.sif case.sif COPYONLY)

  file(COPY ELMERSOLVER_STARTINFO cylinders
       DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

  add_elmer_test(mgdyn_hypre_ams NPROCS 2 LABELS quick hypre)
endif()
