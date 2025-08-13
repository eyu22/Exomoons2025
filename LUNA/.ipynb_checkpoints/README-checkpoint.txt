TO RUN LUNA
Outputs: lightcurve.txt

To edit parameters: change the script test_mandel.f90 and regenerate test_mandel.exe

gfortran -c test_mandel.f90
gfortran occultquad.o mandelmod.o test_mandel.o -o test_mandel.exe
./test_mandel.exe