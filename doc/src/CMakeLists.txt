set(arageli_latex_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/arageli_overview_announce_preprint.tex
)

find_package(LATEX)

if(LATEX_FOUND)
    message(STATUS "latex is found")

    set(DVI_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/../dvi)

    add_custom_target(dvi
        COMMAND ${CMAKE_COMMAND} -E make_directory ${DVI_OUTPUT_DIR}
        COMMAND ${LATEX_COMPILER} -output-directory ${DVI_OUTPUT_DIR} ${arageli_latex_sources}
        COMMAND ${LATEX_COMPILER} -output-directory ${DVI_OUTPUT_DIR} ${arageli_latex_sources}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/build
        VERBATIM
    )

    set(PDF_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/../pdf)
    add_custom_target(pdf
        COMMAND ${CMAKE_COMMAND} -E make_directory ${PDF_OUTPUT_DIR}
        COMMAND ${PDFLATEX_COMPILER} -output-directory ${PDF_OUTPUT_DIR} ${arageli_latex_sources}
        COMMAND ${PDFLATEX_COMPILER} -output-directory ${PDF_OUTPUT_DIR} ${arageli_latex_sources}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/build
        VERBATIM
    )
endif(LATEX_FOUND)
