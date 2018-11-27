project(argparser)

# ArgParser
add_library(argparser STATIC ArgParser.cpp ArgParser.h)
target_include_directories(argparser PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})