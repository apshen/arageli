cmake_minimum_required(VERSION 3.16)

project(arageli CXX)

add_subdirectory(src/arageli)
add_subdirectory(tool/ts/src/ts)
add_subdirectory(perf)
add_subdirectory(test)
add_subdirectory(doc/src)
add_subdirectory(example)

# run script to satisfy Coding Guidelines
add_custom_target(script
	COMMAND python ${CMAKE_SOURCE_DIR}/tool/script/tab_finding.py ${${CMAKE_SOURCE_DIR}}
	COMMAND python ${CMAKE_SOURCE_DIR}/tool/script/crlf_finding.py ${${CMAKE_SOURCE_DIR}}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    #COMMAND python ${CMAKE_SOURCE_DIR}/tool/script/end_spaces_del.py ${${CMAKE_SOURCE_DIR}}
    #COMMAND python ${CMAKE_SOURCE_DIR}./tool/script/crlf2lf.py ${${CMAKE_SOURCE_DIR}}
    #COMMAND python ${CMAKE_SOURCE_DIR}/tool/script/tab_replacing.py ${${CMAKE_SOURCE_DIR}}
)

# update TAGS for source files
add_custom_target(TAGS
	COMMAND ctags -R -h .cpp.hpp
    WORKING_DIRECTORY  ${CMAKE_SOURCE_DIR}
)

# run unit tests
add_custom_target(runtests
    DEPENDS arageli_test
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test/arageli_test
)
