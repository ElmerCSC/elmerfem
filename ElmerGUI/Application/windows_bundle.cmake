IF(CPACK_BUNDLE_EXTRA_WINDOWS_DLLS)
  # Collect list of dependent DLLs using bundledlls
  # This assumes that ElmerGUI is built using MSYS2.
  FIND_PROGRAM(MSYS2_BASH NAMES bash REQUIRED)
  ADD_CUSTOM_COMMAND(
    TARGET ElmerGUI POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bundle/bin
    COMMAND ${MSYS2_BASH} -lc "${CMAKE_CURRENT_SOURCE_DIR}/../../contrib/msys2-bundledlls/bundledlls $<TARGET_FILE:ElmerGUI> ${CMAKE_BINARY_DIR}/bundle/bin"
  )

  INSTALL(DIRECTORY ${CMAKE_BINARY_DIR}/bundle/bin/
    DESTINATION "bin"
    COMPONENT "ELMERGUI")

  # Manually copy Qt platform plugins because they are loaded at runtime
  IF(WITH_QT6)
    EXECUTE_PROCESS(
      COMMAND ${MSYS2_BASH} -lc "cygpath -m \"$\{MINGW_PREFIX}/share/qt6/plugins/platforms\""
      OUTPUT_VARIABLE QT6_PLUGIN_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    INSTALL(DIRECTORY "${QT6_PLUGIN_DIR}/"
      DESTINATION "bin/platforms"
      COMPONENT "ELMERGUI")
  ELSEIF(WITH_QT5)
    EXECUTE_PROCESS(
      COMMAND ${MSYS2_BASH} -lc "cygpath -m \"$\{MINGW_PREFIX}/share/qt5/plugins/platforms\""
      OUTPUT_VARIABLE QT5_PLUGIN_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    INSTALL(DIRECTORY "${QT5_PLUGIN_DIR}/"
      DESTINATION "bin/platforms"
      COMPONENT "ELMERGUI")
  ELSE()
    MESSAGE(WARNING "Qt platform plugins must be bundled manually")
  ENDIF()

  # Additional files (licenses) to be bundled
  SET(LICENSE_FILES "")
  IF(WITH_OCC)
    LIST(APPEND LICENSE_FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/../license_texts/FTL.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/../license_texts/OCE.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/../license_texts/OCE.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/../license_texts/OCCT_LGPL_EXCEPTION.txt
      ${CMAKE_CURRENT_SOURCE_DIR}/../license_texts/LICENSE_LGPL_21.txt)
  ENDIF(WITH_OCC)

  IF(WITH_QWT)
    LIST(APPEND LICENSE_FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/../license_texts/Qwt.txt)
  ENDIF(WITH_QWT)

  IF(WITH_VTK)
    LIST(APPEND LICENSE_FILES
      ${CMAKE_CURRENT_SOURCE_DIR}/../license_texts/vtk-copyright.txt)
  ENDIF(WITH_VTK)

  INSTALL(FILES ${LICENSE_FILES}
    DESTINATION "share/ElmerGUI/license_texts"
    COMPONENT "ELMERGUI")

ENDIF(CPACK_BUNDLE_EXTRA_WINDOWS_DLLS)
