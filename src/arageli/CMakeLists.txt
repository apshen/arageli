# Common compiler flags for all configurations
add_definitions (-DARAGELI_SOURCE -MD)

# Target-specific compilers flags
set (arageli_cxx_flags "-O2 -DARAGELI_DEBUG_LEVEL=1")
set (aragelitest_cxx_flags "-g -DARAGELI_DEBUG_LEVEL=3 -DARAGELI_ASSERT_THROW_EXCEPTION")

# Generate list of source files
aux_source_directory (. root_src)
aux_source_directory (hermite hermite_src)
aux_source_directory (smith smith_src)
set (arageli_src ${root_src} ${smith_src} ${hermite_src})

# Targets
add_library (arageli STATIC ${arageli_src})
add_library (aragelitest STATIC ${arageli_src})

# Associate targets with compiler flags
set_target_properties (arageli PROPERTIES COMPILE_FLAGS ${arageli_cxx_flags})
set_target_properties (aragelitest PROPERTIES COMPILE_FLAGS ${aragelitest_cxx_flags})

