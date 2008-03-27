add_subdirectory (src/arageli)
# ????
# Here I'm trying to point out directory with generated solution and project files.
if (DEFINED WIN32)
    set (CMAKE_BINARY_DIR "build/mcvs2005")
endif (DEFINED WIN32)

# Current version of Cygwin has a certain problems with precompiled headers.
if (DEFINED CYGWIN)
# TODO: Disable precompiled headers!!!
endif (DEFINED CYGWIN)

set (CMAKE_CONFIGURATION_TYPES "Release Debug Test Fast")
set (EXECUTABLE_OUTPUT_PATH bin)
set (LIBRARY_OUTPUT_PATH lib)
set (DEBUG_POSTFIX "d")
set (RELEASE_POSTFIX "r")
# add_subdirectory (test)
