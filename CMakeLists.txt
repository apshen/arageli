set (CMAKE_CONFIGURATION_TYPES Release Debug Test Fast)

# Project name affects name of solution file.
project (arageli)
set (PROJECT_BINARY_DIR bin)

add_subdirectory (src/arageli)
# ????
# Here I'm trying to point out directory with generated solution and project files.
if (WIN32)
    set (CMAKE_BINARY_DIR "build/mcvs2005")
endif (WIN32)

# Current version of Cygwin has a certain problems with precompiled headers.
if (CYGWIN)
# TODO: Disable precompiled headers!!!
endif (CYGWIN)

# add_subdirectory (test)
set (EXECUTABLE_OUTPUT_PATH bin)
set (LIBRARY_OUTPUT_PATH lib)
set (DEBUG_POSTFIX "d")
set (RELEASE_POSTFIX "r")

