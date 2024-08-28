
# Some versions of VTK don't correctly instantiate extern template classes
# for some platforms (notably Windows). Check whether that is an issue with the
# used version of VTK.

message(STATUS "Checking whether GetValue works for vtkIntArray")

file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_vtkIntArray_GetValue.cpp
  "
  #include <vtkIntArray.h>
  #include <vtkNew.h>

  int main(void)
  {
    vtkNew<vtkIntArray> vtk_int_array;
    vtk_int_array->InsertNextValue(42);
    vtk_int_array->GetValue(0);

    return 0;
  }
  ")
try_compile(CHECK_VTKINTARRAY_GETVALUE_BUILD_SUCCESS ${CMAKE_BINARY_DIR}
  SOURCES ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_vtkIntArray_GetValue.cpp
  CMAKE_FLAGS -DINCLUDE_DIRECTORIES=${VTK_INCLUDE_DIRS}
  LINK_LIBRARIES ${VTK_LIBRARIES})

if(CHECK_VTKINTARRAY_GETVALUE_BUILD_SUCCESS)
  message(STATUS "Checking whether GetValue works for vtkIntArray -- yes")
  set(ELMER_INSTANTIATE_VTK_ARRAY_TEMPLATE OFF CACHE INTERNAL "")
else()
  message(STATUS "Checking whether GetValue works for vtkIntArray -- no")
  set(ELMER_INSTANTIATE_VTK_ARRAY_TEMPLATE ON CACHE INTERNAL "")
endif()
