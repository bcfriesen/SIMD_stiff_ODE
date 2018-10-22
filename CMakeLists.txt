cmake_minimum_required (VERSION 2.8.11)
project (SIMD_stiff_ODE)

enable_language(Fortran)

set (CMAKE_VERBOSE_MAKEFILE "TRUE")

add_subdirectory (src)
