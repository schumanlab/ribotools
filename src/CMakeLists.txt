project(main)

add_subdirectory(argparser)
add_subdirectory(bedrecord)
add_subdirectory(metagene)

# main program
add_executable(main neuroribotools.cpp)
target_link_libraries(main PRIVATE argparser metagene)