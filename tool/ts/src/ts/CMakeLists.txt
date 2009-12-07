# Common compiler flags for all configurations
add_definitions (-MD -g -O2)

# Generate list of source files
aux_source_directory (. ts_src)

# Targets
add_library (ts STATIC ${ts_src})

