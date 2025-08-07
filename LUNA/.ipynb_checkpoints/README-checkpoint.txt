gfortran -c test_mandel.f90
gfortran occultquad.o mandelmod.o test_mandel.o -o test_mandel.exe
./test_mandel.exe