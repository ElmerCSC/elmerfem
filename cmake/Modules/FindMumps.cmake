# cmake script for finding MUMPS sparse direct solver
# Supports both parallel MUMPS (MPI + ScaLAPACK) and sequential MUMPS
# (compiled with MUMPS's own libmpiseq stub, no MPI or ScaLAPACK needed).
#
# On Debian/Ubuntu with libmumps-seq-dev installed, sequential libraries are
# found automatically. Without the dev package, set MUMPS_ROOT or MUMPSROOT to
# point to a custom MUMPS sequential build.
INCLUDE(FindPackageHandleStandardArgs)

# If Mumps libraries are already defined, do nothing
IF(Mumps_LIBRARIES AND Mumps_INCLUDE_DIR)
   SET(Mumps_FOUND TRUE)
   RETURN()
ENDIF()

# Try to find with CMake config file of upstream mumps-superbuild.
FIND_PACKAGE(MUMPS CONFIG)

IF(MUMPS_FOUND)
  SET(Mumps_FOUND TRUE)
  SET(Mumps_LIBRARIES ${MUMPS_LIBRARIES})
  SET(Mumps_INCLUDE_DIR ${MUMPS_INCLUDE_DIRS})
  IF(MUMPS_METIS_FOUND)
    SET(Metis_FOUND TRUE)
    SET(Metis_LIBRARIES ${METIS_LIBRARIES})
    SET(Metis_INCLUDE_DIR ${METIS_INCLUDE_DIRS})
  ENDIF()

  RETURN()
ENDIF()

# Fall back to manual search
SET(Mumps_FOUND FALSE)
MESSAGE(STATUS "Finding Mumps")

# ── Search hints ──────────────────────────────────────────────────────────────
SET(MUMPSINCLUDE
  "${MUMPSROOT}/include"
  "$ENV{MUMPSROOT}/include"
  "${MUMPS_ROOT}/include"
  "$ENV{MUMPS_ROOT}/include"
  "${CMAKE_SOURCE_DIR}/mumps/include"
  INTERNAL)

SET(MUMPSLIB
  "${MUMPSROOT}/lib"
  "$ENV{MUMPSROOT}/lib"
  "${MUMPS_ROOT}/lib"
  "$ENV{MUMPS_ROOT}/lib"
  "${CMAKE_SOURCE_DIR}/mumps/lib"
  INTERNAL)

# Debian/Ubuntu system library directories
SET(MUMPS_SYSTEM_LIBDIRS
  /usr/lib/x86_64-linux-gnu
  /usr/lib/aarch64-linux-gnu
  /usr/lib/arm-linux-gnueabihf
  /usr/lib64
  /usr/lib)

# ── Find headers ──────────────────────────────────────────────────────────────
FIND_PATH(Mumps_INCLUDE_DIR
  NAMES dmumps_struc.h zmumps_struc.h smumps_struc.h cmumps_struc.h
  HINTS ${MUMPSINCLUDE})

# ── Helper macro: find a library by name, with Debian versioned-name fallback ─
# On Debian/Ubuntu, sequential MUMPS uses libname-X.Y.so (version before .so),
# which FIND_LIBRARY does not recognise without an unversioned symlink.
MACRO(FIND_MUMPS_LIB outvar libname)
  FIND_LIBRARY(${outvar} NAMES ${libname} HINTS ${MUMPSLIB})
  IF(NOT ${outvar})
    FOREACH(_dir ${MUMPSLIB} ${MUMPS_SYSTEM_LIBDIRS})
      IF(NOT ${outvar})
        FILE(GLOB _cands
          "${_dir}/lib${libname}-*.so"
          "${_dir}/lib${libname}-*.so.*"
          "${_dir}/lib${libname}-*.a")
        IF(_cands)
          LIST(SORT _cands)
          LIST(GET _cands 0 ${outvar})
          SET(${outvar} "${${outvar}}" CACHE FILEPATH
              "Path to ${libname} library" FORCE)
          BREAK()
        ENDIF()
      ENDIF()
    ENDFOREACH()
    UNSET(_cands)
    UNSET(_dir)
  ENDIF()
ENDMACRO()

# ── Sequential vs parallel path ───────────────────────────────────────────────
IF(MPI_FOUND)
  # ── PARALLEL MUMPS ──────────────────────────────────────────────────────────
  FIND_MUMPS_LIB(MUMPS_D_LIB      dmumps)
  FIND_MUMPS_LIB(MUMPS_Z_LIB      zmumps)
  FIND_MUMPS_LIB(MUMPS_S_LIB      smumps)
  FIND_MUMPS_LIB(MUMPS_C_LIB      cmumps)
  FIND_MUMPS_LIB(MUMPS_COMMON_LIB mumps_common)
  FIND_MUMPS_LIB(MUMPS_PORD_LIB   pord)

  IF(Mumps_INCLUDE_DIR AND MUMPS_D_LIB AND MUMPS_Z_LIB AND
     MUMPS_COMMON_LIB AND MUMPS_PORD_LIB)
    SET(Mumps_LIBRARIES
      ${MUMPS_D_LIB} ${MUMPS_Z_LIB} ${MUMPS_S_LIB} ${MUMPS_C_LIB}
      ${MUMPS_COMMON_LIB} ${MUMPS_PORD_LIB})

    # Parallel MUMPS always needs ScaLAPACK
    FIND_PACKAGE(SCALAPACK QUIET)
    IF(SCALAPACK_FOUND)
      LIST(APPEND Mumps_LIBRARIES ${SCALAPACK_LIBRARIES})

      # Check for Metis.
      #
      # Use "nm -D" first: a stripped shared library has no static symbol
      # table, so plain "nm" prints nothing at all and reports "no symbols" on
      # stderr. Requiring stderr to be empty, as this test used to, therefore
      # made the "not needed" branch unreachable for any shared library, and
      # Metis and ParMetis were always declared mandatory. On Debian/Ubuntu
      # that leaves parallel Mumps unconfigurable, even though the packaged
      # Mumps is built against PT-Scotch and references neither.
      MESSAGE(STATUS "Checking if Metis library is needed by Mumps")
      EXECUTE_PROCESS(COMMAND ${CMAKE_NM} -D ${MUMPS_D_LIB}
        OUTPUT_VARIABLE _nm_out ERROR_VARIABLE _nm_err)
      IF("${_nm_out}" STREQUAL "")
        # Static archive, or an object format with no dynamic symbol table.
        EXECUTE_PROCESS(COMMAND ${CMAKE_NM} ${MUMPS_D_LIB}
          OUTPUT_VARIABLE _nm_out ERROR_VARIABLE _nm_err)
      ENDIF()
      STRING(FIND "${_nm_out}" "metis_nodend" _metis_pos)
      IF("${_metis_pos}" STREQUAL "-1")
        MESSAGE(STATUS "Checking if Metis library is needed by Mumps -- no")
      ELSE()
        MESSAGE(STATUS "Checking if Metis library is needed by Mumps -- yes")
        IF(EXTERNAL_METIS)
          FIND_PACKAGE(Metis QUIET)
          IF(Metis_FOUND)
            LIST(APPEND Mumps_LIBRARIES ${Metis_LIBRARIES})
            LIST(APPEND Mumps_INCLUDE_DIR ${Metis_INCLUDE_DIR})
          ELSE()
            SET(MUMPS_FAILMSG "Metis not found, needed by Mumps.")
          ENDIF()
        ELSE()
          MESSAGE(STATUS "Using bundled Metis")
        ENDIF()
      ENDIF()

      # Check for ParMetis
      IF(NOT MUMPS_FAILMSG)
        # See the Metis check above for why "nm -D" is used here.
        MESSAGE(STATUS "Checking if ParMetis library is needed by Mumps")
        EXECUTE_PROCESS(COMMAND ${CMAKE_NM} -D ${MUMPS_COMMON_LIB}
          OUTPUT_VARIABLE _nm_out ERROR_VARIABLE _nm_err)
        IF("${_nm_out}" STREQUAL "")
          EXECUTE_PROCESS(COMMAND ${CMAKE_NM} ${MUMPS_COMMON_LIB}
            OUTPUT_VARIABLE _nm_out ERROR_VARIABLE _nm_err)
        ENDIF()
        STRING(FIND "${_nm_out}" "ParMETIS_V3_NodeND" _parmetis_pos)
        IF("${_parmetis_pos}" STREQUAL "-1")
          MESSAGE(STATUS "Checking if ParMetis library is needed by Mumps -- no")
        ELSE()
          MESSAGE(STATUS "Checking if ParMetis library is needed by Mumps -- yes")
          FIND_PACKAGE(ParMetis QUIET)
          IF(ParMetis_FOUND)
            LIST(APPEND Mumps_LIBRARIES ${ParMetis_LIBRARIES})
            LIST(APPEND Mumps_INCLUDE_DIR ${ParMetis_INCLUDE_DIR})
          ELSE()
            SET(MUMPS_FAILMSG "ParMetis not found, needed by Mumps.")
          ENDIF()
        ENDIF()
      ENDIF()

    ELSE()
      SET(MUMPS_FAILMSG "ScaLAPACK not found, required by parallel Mumps.")
    ENDIF()
  ELSE()
    SET(MUMPS_FAILMSG "Parallel Mumps libraries not found.")
  ENDIF()

ELSE()
  # ── SEQUENTIAL MUMPS (no real MPI) ─────────────────────────────────────────
  # On Debian/Ubuntu, sequential libraries carry a _seq suffix:
  #   libdmumps_seq-5.8.so, libmumps_common_seq-5.8.so, libmpiseq_seq-5.8.so
  # Standard MUMPS builds put libmpiseq in lib/ alongside the main libs.

  FIND_MUMPS_LIB(MUMPS_D_SEQ_LIB      dmumps_seq)
  FIND_MUMPS_LIB(MUMPS_Z_SEQ_LIB      zmumps_seq)
  FIND_MUMPS_LIB(MUMPS_S_SEQ_LIB      smumps_seq)
  FIND_MUMPS_LIB(MUMPS_C_SEQ_LIB      cmumps_seq)
  FIND_MUMPS_LIB(MUMPS_COMMON_SEQ_LIB mumps_common_seq)
  FIND_MUMPS_LIB(MUMPS_PORD_SEQ_LIB   pord_seq)

  # The sequential MPI stub: try mpiseq_seq (Debian) and mpiseq (upstream)
  FIND_MUMPS_LIB(MUMPS_MPISEQ_LIB mpiseq_seq)
  IF(NOT MUMPS_MPISEQ_LIB)
    FIND_MUMPS_LIB(MUMPS_MPISEQ_LIB mpiseq)
  ENDIF()

  IF(Mumps_INCLUDE_DIR AND MUMPS_D_SEQ_LIB AND MUMPS_Z_SEQ_LIB AND
     MUMPS_COMMON_SEQ_LIB AND MUMPS_PORD_SEQ_LIB)

    SET(Mumps_LIBRARIES
      ${MUMPS_D_SEQ_LIB} ${MUMPS_Z_SEQ_LIB} ${MUMPS_S_SEQ_LIB}
      ${MUMPS_C_SEQ_LIB} ${MUMPS_COMMON_SEQ_LIB} ${MUMPS_PORD_SEQ_LIB})

    IF(MUMPS_MPISEQ_LIB)
      LIST(APPEND Mumps_LIBRARIES ${MUMPS_MPISEQ_LIB})
      MESSAGE(STATUS "Found sequential MUMPS mpiseq stub: ${MUMPS_MPISEQ_LIB}")
    ELSE()
      MESSAGE(STATUS "mpiseq stub not found; Elmer's MPI stubs will provide "
                     "the Fortran-level MPI symbols for sequential MUMPS.")
    ENDIF()

    # Sequential mpif.h stub (needed so MUMPS headers compile cleanly)
    FIND_PATH(MUMPS_SEQ_INCLUDE_DIR mpif.h
      HINTS
        /usr/include/mumps_seq
        "${MUMPSROOT}/libseq"
        "$ENV{MUMPSROOT}/libseq"
        "${MUMPS_ROOT}/libseq"
        "$ENV{MUMPS_ROOT}/libseq")
    IF(MUMPS_SEQ_INCLUDE_DIR)
      LIST(APPEND Mumps_INCLUDE_DIR ${MUMPS_SEQ_INCLUDE_DIR})
    ENDIF()

    LIST(REMOVE_DUPLICATES Mumps_LIBRARIES)
    LIST(REMOVE_DUPLICATES Mumps_INCLUDE_DIR)

  ELSE()
    SET(MUMPS_FAILMSG
      "Sequential MUMPS (_seq) libraries not found. "
      "Install libmumps-seq-dev (Debian/Ubuntu) or set MUMPS_ROOT to a "
      "sequential MUMPS build directory.")
  ENDIF()
ENDIF()

# ── Report ─────────────────────────────────────────────────────────────────────
IF(NOT MUMPS_FAILMSG)
  SET(Mumps_FOUND TRUE)
ENDIF()

IF(Mumps_FOUND)
  IF(NOT Mumps_FIND_QUIETLY)
    MESSAGE(STATUS "Found Mumps:")
    MESSAGE(STATUS "  Include dirs: ${Mumps_INCLUDE_DIR}")
    MESSAGE(STATUS "  Libraries:    ${Mumps_LIBRARIES}")
  ENDIF()
ELSE()
  IF(Mumps_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR ${MUMPS_FAILMSG})
  ELSE()
    MESSAGE(STATUS "Mumps not found: ${MUMPS_FAILMSG}")
  ENDIF()
ENDIF()

MARK_AS_ADVANCED(
  MUMPSINCLUDE MUMPSLIB
  MUMPS_FAILMSG
  Mumps_INCLUDE_DIR Mumps_LIBRARIES
  MUMPS_D_LIB MUMPS_Z_LIB MUMPS_S_LIB MUMPS_C_LIB
  MUMPS_COMMON_LIB MUMPS_PORD_LIB
  MUMPS_D_SEQ_LIB MUMPS_Z_SEQ_LIB MUMPS_S_SEQ_LIB MUMPS_C_SEQ_LIB
  MUMPS_COMMON_SEQ_LIB MUMPS_PORD_SEQ_LIB
  MUMPS_MPISEQ_LIB MUMPS_SEQ_INCLUDE_DIR
  SCALAPACK_LIBRARIES)
