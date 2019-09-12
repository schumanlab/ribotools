cmake_minimum_required(VERSION 3.0)

if(APPLE)
    # Fix linking on 10.14+ See https://stackoverflow.com/questions/54068035
    link_directories(/usr/local/lib)
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall")
set (CMAKE_CXX_STANDARD 14)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++14 -W -Wall -pedantic")
endif()

project(ribotools)
set(SOURCE 
    ribotools.cpp
    aminoacids.cpp aminoacids.h
    bamhandle.cpp bamhandle.h
    bedrecord.cpp bedrecord.h
    parserargv.cpp parserargv.h
    main_codonfreq.cpp
    main_codonrate.cpp
    main_count.cpp
    main_features.cpp
    main_gcratio.cpp
    main_irate.cpp
    main_length.cpp
    main_metagene.cpp
    main_mtdr.cpp
    main_pausing.cpp
    main_poffset.cpp
    main_uorfs.cpp
    main_utrseq.cpp
    )
add_executable(${PROJECT_NAME} ${SOURCE})



find_package(PkgConfig REQUIRED)

pkg_check_modules(ZLIB REQUIRED zlib)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${ZLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${ZLIB_LIBRARIES})

pkg_check_modules(HTSLIB REQUIRED htslib)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${HTSLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${HTSLIB_LIBRARIES})
