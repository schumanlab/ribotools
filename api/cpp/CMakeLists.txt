cmake_minimum_required(VERSION 3.0)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall")
set (CMAKE_CXX_STANDARD 17)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++17 -W -Wall -pedantic")
endif()

project(codonNFC)
set(SOURCE codonNFC.cpp
           bedrecord.cpp bedrecord.hpp)
add_executable(${PROJECT_NAME} ${SOURCE})

find_package(PkgConfig REQUIRED)

pkg_check_modules(ZLIB REQUIRED zlib)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${ZLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${ZLIB_LIBRARIES})

pkg_check_modules(HTSLIB REQUIRED htslib)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${HTSLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${HTSLIB_LIBRARIES})
