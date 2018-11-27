# Metagene
add_executable(metagene MetaGeneMain.cpp MetaGene.cpp MetaGene.h)
target_link_libraries(metagene argparser bedrecord)

# PkgConfig to resolve external dependencies
find_package(PkgConfig REQUIRED)

# ZLIB
pkg_check_modules(ZLIB REQUIRED zlib)
target_include_directories(metagene SYSTEM PUBLIC ${ZLIB_INCLUDE_DIRS})
target_link_libraries(metagene ${ZLIB_LIBRARIES})

# HTSLIB
pkg_check_modules(HTSLIB REQUIRED htslib)
target_include_directories(metagene SYSTEM PUBLIC ${HTSLIB_INCLUDE_DIRS})
target_link_libraries(metagene ${HTSLIB_LIBRARIES})
