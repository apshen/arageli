file(GLOB arageli_examples_list
    "*.cpp"
)

add_custom_target(example)

foreach(example_source ${arageli_examples_list})
    get_filename_component(bin_name ${example_source} NAME_WE)
    add_executable(${bin_name} EXCLUDE_FROM_ALL ${example_source})
    target_include_directories(${bin_name} PRIVATE ../src)
    target_link_libraries(${bin_name} arageli)

    add_dependencies(example ${bin_name})
endforeach()

