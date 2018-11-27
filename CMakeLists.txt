cmake_minimum_required(VERSION 3.12)

project(neuroribotools VERSION 1.0.0 LANGUAGES CXX)
target_compile_options(neuroribotools PRIVATE -Wall)
target_compile_options(neuroribotools PRIVATE cxx_std_17)
add_subdirectory(src)
