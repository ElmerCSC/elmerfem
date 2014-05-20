
# Quick hack to test if the Fortran compiler supports the EXECUTE_COMMAND_LINE

message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports EXECUTE_COMMAND_LINE")
file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompilerExecCommand.f90
  "
      PROGRAM TESTFortranExecCommand
      CALL EXECUTE_COMMAND_LINE('echo Hello World!',.TRUE.)
      END PROGRAM TESTFortranExecCommand
  ")
try_compile(FC_HAS_EXECUTECOMMANDLINE ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testFortranCompilerExecCommand.f90
  OUTPUT_VARIABLE OUTPUT)
if(FC_HAS_EXECUTECOMMANDLINE)
  message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports EXECUTE_COMMAND_LINE -- yes")
  set(CMAKE_Fortran_COMPILER_SUPPORTS_EXECUTECOMMANDLINE 1 CACHE INTERNAL "")
else(FC_HAS_EXECOMMAND)
  message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports EXECUTE_COMMAND_LINE -- no")
  set(CMAKE_Fortran_COMPILER_SUPPORTS_EXECUTECOMMANDLINE 0 CACHE INTERNAL "")
endif(FC_HAS_EXECUTECOMMANDLINE)