add_executable (main
                main.f90
                rhs_simd.f90
                rhs_scalar.f90
                rkf45_simd.f90
                rkf45_scalar.f90)
