project(bedrecord)

# BedRecord
add_library(bedrecord STATIC BedRecord.cpp BedRecord.h)
target_include_directories(bedrecord PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})