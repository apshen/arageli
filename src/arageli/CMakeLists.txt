set(arageli_lib_sources
    advsimplmeth.cpp
    algebraic.cpp
    algebrslt.cpp
    bareiss.cpp
    berlekamp_hensel.cpp
    bigar.cpp
    big_const.cpp
    big_float.cpp
    big_int.cpp
    cantor_zassenhaus.cpp
    cone.cpp
    exception.cpp
    gauss.cpp
    gcd.cpp
    generic.cpp
    gram_schmidt.cpp
    groebner_basis.cpp
    hermite.cpp
    hermite/classic.cpp
    hermite/domich.cpp
    hermite/hafner.cpp
    hermite/storjohann.cpp
    intalg.cpp
    intconvex.cpp
    intcount_barvinok.cpp
    io.cpp
    karatsuba.cpp
    linopt.cpp
    lll.cpp
    logarithm.cpp
    matrix.cpp
    motzkin_burger.cpp
    orthogonal.cpp
    pollard_mult.cpp
    polyhedron.cpp
    polyhedron_relation.cpp
    polynom.cpp
    powerest.cpp
    prime.cpp
    rational.cpp
    resultant.cpp
    sideset.cpp
    smith.cpp
    smith/classic.cpp
    smith/near_optimal.cpp
    smithpoly.cpp
    smith/storjohann.cpp
    solve_linsys.cpp
    sparse_multipolynom.cpp
    sparse_polynom.cpp
    sturm.cpp
    tex_cartesian2d.cpp
    texout.cpp
    timer.cpp
    triangulation.cpp
    type_opers.cpp
    type_traits.cpp
    _utility.cpp
    vecalg.cpp
    vector.cpp
    winograd.cpp
)

# arageli lib
add_library(arageli STATIC ${arageli_lib_sources})
target_compile_definitions(arageli PRIVATE
    ARAGELI_SOURCE
    ARAGELI_DEBUG_LEVEL=1
    ARAGELI_ASSERT_THROW_EXCEPTION
)
target_compile_options(arageli PRIVATE
    -O2
    -Wall
)

install(TARGETS arageli)
install(DIRECTORY .
    DESTINATION include/arageli
    FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "*.cpp"
)

# testarageli lib
add_library(testarageli STATIC EXCLUDE_FROM_ALL ${arageli_lib_sources})
target_compile_definitions(testarageli PRIVATE
    ARAGELI_SOURCE
    ARAGELI_DEBUG_LEVEL=3
    ARAGELI_ASSERT_THROW_EXCEPTION
)
target_compile_options(testarageli PRIVATE
    -g
    -Wall
)

# perfarageli lib
add_library(perfarageli STATIC EXCLUDE_FROM_ALL ${arageli_lib_sources})
target_compile_definitions(perfarageli PRIVATE
    ARAGELI_SOURCE
    ARAGELI_DEBUG_LEVEL=0
    ARAGELI_ASSERT_THROW_EXCEPTION
)
target_compile_options(perfarageli PRIVATE
    -O3
    -pg
    -g
    -Wall
)

# Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
    message(STATUS "Doxygen is found")

    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(HTML
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif(DOXYGEN_FOUND)
