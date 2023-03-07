set(arageli_karatsuba_vs_classic_sources
    karatsuba_vs_classic.cpp
)

add_executable(karatsuba_vs_classic EXCLUDE_FROM_ALL ${arageli_karatsuba_vs_classic_sources})

target_compile_definitions(karatsuba_vs_classic PRIVATE
    ARAGELI_SOURCE
)

target_compile_options(karatsuba_vs_classic PRIVATE
    -pg
    -O3
)

target_include_directories(karatsuba_vs_classic PRIVATE
    ../src
)

target_link_libraries(karatsuba_vs_classic perfarageli)
