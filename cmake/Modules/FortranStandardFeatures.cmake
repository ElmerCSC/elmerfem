# fem/include/elmer_fortran_features.h decides the ELMER_HAVE_F20xx_* macros
# from compiler version macros. This only checks, with the real build flags,
# that everything the header enables actually compiles; a mismatch usually
# means a -std flag passed by hand or an untested compiler build.
#
# Sets ELMER_FORTRAN_FEATURES_ENABLED, and a variable named after each macro,
# for the features the header enables.

INCLUDE(CheckFortranSourceCompiles)

SET(ELMER_FORTRAN_FEATURES_INCLUDE_DIR "${CMAKE_CURRENT_LIST_DIR}/../../fem/include")
GET_FILENAME_COMPONENT(ELMER_FORTRAN_FEATURES_INCLUDE_DIR
  "${ELMER_FORTRAN_FEATURES_INCLUDE_DIR}" ABSOLUTE)

SET(_elmer_saved_required_flags "${CMAKE_REQUIRED_FLAGS}")
SET(_elmer_saved_required_includes "${CMAKE_REQUIRED_INCLUDES}")
SET(_elmer_saved_required_definitions "${CMAKE_REQUIRED_DEFINITIONS}")
SET(CMAKE_REQUIRED_FLAGS
  "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_UPCASE}}")
SET(CMAKE_REQUIRED_INCLUDES "${ELMER_FORTRAN_FEATURES_INCLUDE_DIR}")
IF(DEFINED ELMER_FORTRAN_STD_DEFINITION)
  SET(CMAKE_REQUIRED_DEFINITIONS "${ELMER_FORTRAN_STD_DEFINITION}")
ENDIF()

SET(ELMER_FORTRAN_FEATURES_ENABLED "")

# Results are not cached: the answer changes with the compiler flags.
MACRO(ELMER_CHECK_FORTRAN_FEATURE _feat _source)
  UNSET(${_feat}_CLAIMED CACHE)
  UNSET(${_feat}_WORKS CACHE)
  CHECK_FORTRAN_SOURCE_COMPILES("
#include \"elmer_fortran_features.h\"
#ifndef ${_feat}
#error not enabled
#endif
program probe
end program probe"
    ${_feat}_CLAIMED SRC_EXT F90)
  IF(${_feat}_CLAIMED)
    CHECK_FORTRAN_SOURCE_COMPILES("${_source}" ${_feat}_WORKS SRC_EXT F90)
    IF(NOT ${_feat}_WORKS)
      MESSAGE(FATAL_ERROR "elmer_fortran_features.h enables ${_feat}, but "
        "${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} rejects it "
        "with the flags '${CMAKE_REQUIRED_FLAGS}'. If you pass -std by hand, use "
        "ELMER_FORTRAN_STANDARD instead.")
    ENDIF()
    LIST(APPEND ELMER_FORTRAN_FEATURES_ENABLED ${_feat})
    SET(${_feat} TRUE)
  ELSE()
    SET(${_feat} FALSE)
  ENDIF()
ENDMACRO()

# The internal-procedure call through a dummy procedure pointer is what
# gfortran 16 gets wrong.
ELMER_CHECK_FORTRAN_FEATURE(ELMER_HAVE_F2018_IMPLICIT_NONE_EXTERNAL "
module probe_m
  implicit none (type, external)
  abstract interface
    real function f_iface(x)
      real, intent(in) :: x
    end function f_iface
  end interface
contains
  real function outer(f, x)
    implicit none (type, external)
    procedure(f_iface), pointer :: f
    real, intent(in) :: x
    outer = inner()
  contains
    real function inner()
      implicit none (type, external)
      inner = f(x)
    end function inner
  end function outer
end module probe_m
program probe
  implicit none (type, external)
  external :: probe_s
  call probe_s()
end program probe
subroutine probe_s()
end subroutine probe_s")

ELMER_CHECK_FORTRAN_FEATURE(ELMER_HAVE_F2018_DO_CONCURRENT_LOCALITY "
program probe
  implicit none
  integer :: i
  real :: tmp, a(4)
  do concurrent (i = 1:4) local(tmp) shared(a)
    tmp = real(i)
    a(i) = tmp
  end do
end program probe")

ELMER_CHECK_FORTRAN_FEATURE(ELMER_HAVE_F2023_DO_CONCURRENT_REDUCE "
program probe
  implicit none
  integer :: i
  real :: s
  s = 0.0
  do concurrent (i = 1:4) reduce(+:s)
    s = s + real(i)
  end do
end program probe")

SET(CMAKE_REQUIRED_FLAGS "${_elmer_saved_required_flags}")
SET(CMAKE_REQUIRED_INCLUDES "${_elmer_saved_required_includes}")
SET(CMAKE_REQUIRED_DEFINITIONS "${_elmer_saved_required_definitions}")

MESSAGE(STATUS "------------------------------------------------")
MESSAGE(STATUS "Optional Fortran features (F2008 is the required floor):")
FOREACH(_feat
    ELMER_HAVE_F2018_IMPLICIT_NONE_EXTERNAL
    ELMER_HAVE_F2018_DO_CONCURRENT_LOCALITY
    ELMER_HAVE_F2023_DO_CONCURRENT_REDUCE)
  IF(_feat IN_LIST ELMER_FORTRAN_FEATURES_ENABLED)
    MESSAGE(STATUS "  ${_feat} -- yes")
  ELSE()
    MESSAGE(STATUS "  ${_feat} -- no")
  ENDIF()
ENDFOREACH()
UNSET(_feat)
MESSAGE(STATUS "------------------------------------------------")
