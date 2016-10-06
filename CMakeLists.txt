cmake_minimum_required (VERSION 2.8.11)
project (SIMD_stiff_ODE)

enable_language(Fortran)

set (CMAKE_Fortran_FLAGS "-qopt-report=5 -qopt-report-phase=vec -qopenmp")

add_subdirectory (src)
