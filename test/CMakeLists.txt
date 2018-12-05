project(ribotoolsTests VERSION 1.0.0 LANGUAGES CXX)
set(SOURCE  test_tabix.cpp)
add_executable(${PROJECT_NAME} ${SOURCE})

find_package(PkgConfig REQUIRED)

pkg_check_modules(ZLIB REQUIRED zlib)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${ZLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${ZLIB_LIBRARIES})

pkg_check_modules(HTSLIB REQUIRED htslib)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${HTSLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${HTSLIB_LIBRARIES})