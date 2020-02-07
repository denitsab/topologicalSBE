#!/bin/bash
#gfortran -O3 -g -fbounds-check -Wall -ftracer -L/usr/local/Cellar/lapack/3.8.0_1/lib -lblas -llapack -L/Users/denitsab/Programs/Haldane\ Model\ SBE/zvode_adapted_step -llinpack_z zgbfa.f zgesl.f zgefa.f zvode.f sbe_bi2se3_B.f90 -o sbe.x lapack.a blas.a liblinpack_z.a liblapack.a libblas.a
#mpif90 -fopenmp -O3 -g -fbounds-check -Wall -ftracer -lopenblas zgbfa.f zgesl.f zgefa.f zvode.f sbe_bi2se3_B.f90 -o sbe.x
mpif90 -fopenmp -O3 -g -fbounds-check -Wall -ftracer -lblas -llapack -L/cygdrive/d/Haldane\ Model\ SBE/Bi2Se3\ Bulk zgbfa.f zgesl.f zgefa.f zvode.f math.f90 legendre.f90 phys.f90 iostream.f90 grid.f90 field.f90 solid.f90 bloch.f90 info.f90 main.f90 -o unisbe.x lapack.a blas.a liblinpack_z.a liblapack.a libblas.a

