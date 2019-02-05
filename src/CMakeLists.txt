project(ribotools VERSION 1.0.0 LANGUAGES CXX)
set(SOURCE  ribotools.cpp
            version.cpp version.hpp
            bedrecord.cpp bedrecord.hpp
            gbedrecord.cpp gbedrecord.hpp
            parserbed.cpp parserbed.hpp
            parsergbed.cpp parsergbed.hpp
            metageneMain.cpp metagene.cpp metagene.hpp
            extractumi.cpp)
add_executable(${PROJECT_NAME} ${SOURCE})

find_package(PkgConfig REQUIRED)

pkg_check_modules(ZLIB REQUIRED zlib)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${ZLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${ZLIB_LIBRARIES})

pkg_check_modules(HTSLIB REQUIRED htslib)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${HTSLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${HTSLIB_LIBRARIES})