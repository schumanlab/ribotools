project(argparser)
add_library(argparser STATIC argparser.cpp argparser.h)
target_include_directories(argparser PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})