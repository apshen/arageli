# Common compiler flags for all configurations
add_definitions (-MD -g -DARAGELI_SOURCE -O2 -DARAGELI_DEBUG_LEVEL=3 -DARAGELI_ASSERT_THROW_EXCEPTION)

# Generate list of source files
aux_source_directory (. root_src)
aux_source_directory (ITLab-produced ITLab_src)
set (test_src ${root_src} ${ITLab_src})

include_directories (../src ../tool/ts/src)
link_directories (../lib)

# Targets
add_executable (test ${test_src})
target_link_libraries(test ts aragelitest)
