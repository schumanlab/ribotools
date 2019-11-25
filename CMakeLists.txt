cmake_minimum_required(VERSION 3.0)

if(APPLE)
    # Fix linking on 10.14+ See https://stackoverflow.com/questions/54068035
    link_directories(/usr/local/lib)
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall")
set (CMAKE_CXX_STANDARD 17)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++17 -W -Wall -pedantic")
endif()

project(ribotools)
set(SOURCE 
    ribotools.cpp
    argumentparser.h
    aminoacidtable.cpp aminoacidtable.h
    bamhandle.cpp bamhandle.h
    bamio.cpp bamio.h
    bedio.cpp bedio.h
    bedrecord.cpp bedrecord.h
    seqio.cpp seqio.h
    version.cpp version.h
    main_basefreq.cpp
    main_cai.cpp
    main_codonfreq.cpp
    main_codonrate.cpp
    main_count.cpp
    main_depth.cpp
    main_features.cpp
    main_gcratio.cpp
    main_gcref.cpp
    main_irate.cpp
    main_length.cpp
    main_metagene.cpp
    main_mtdr.cpp
    main_pausing.cpp
    main_periodicity.cpp
    main_poffset.cpp
    main_runoff.cpp
    main_structure.cpp
    main_translate.cpp
    main_uorfs.cpp
    main_utrseq.cpp
    parserargv.cpp parserargv.h
    )
add_executable(${PROJECT_NAME} ${SOURCE})



find_package(PkgConfig REQUIRED)

pkg_check_modules(ZLIB REQUIRED zlib)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${ZLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${ZLIB_LIBRARIES})

pkg_check_modules(HTSLIB REQUIRED htslib)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${HTSLIB_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${HTSLIB_LIBRARIES})
