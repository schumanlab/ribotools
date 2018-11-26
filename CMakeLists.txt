cmake_minimum_required(VERSION 3.2)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++17 -W -Wall -pedantic")
endif()


project(neuroribotools)
add_subdirectory(src)
add_subdirectory(test)

# external dependeincies
find_package(PkgConfig REQUIRED)

pkg_check_modules(ZLIB REQUIRED zlib)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${ZLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${ZLIB_LIBRARIES})

pkg_check_modules(HTSLIB REQUIRED htslib)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${HTSLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${HTSLIB_LIBRARIES})