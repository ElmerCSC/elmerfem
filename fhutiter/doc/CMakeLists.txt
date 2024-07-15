FIND_PACKAGE(LATEX COMPONENTS PDFLATEX)
IF(LATEX_PDFLATEX_FOUND)
    ADD_CUSTOM_COMMAND(
        OUTPUT hutidoc.pdf
        COMMAND "${PDFLATEX_COMPILER}" "${CMAKE_CURRENT_SOURCE_DIR}/hutidoc.tex"
        DEPENDS hutidoc.tex
    )
    ADD_CUSTOM_TARGET(fhutiter_docs DEPENDS hutidoc.pdf)
    MESSAGE(STATUS "PdfLaTeX found, use fhutiter_docs target to generate documentation")
ELSE()
    MESSAGE(STATUS "No PdfLaTeX found, hutidoc.pdf cannot be generated.")
ENDIF()
