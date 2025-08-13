TO RUN LUNA
Outputs: lightcurve.txt

To edit parameters: change the script test_mandel.f90 and regenerate test_mandel.exe
You can change the range of impact parameters that we generate flux for 


PLANET ONLY (TEST)
gfortran -c test_mandel.f90
gfortran occultquad.o mandelmod.o test_mandel.o -o test_mandel.exe
./test_mandel.exe


PLANET + MOON
gfortran -O3 -o gen_pm generate_pm_lightcurve.f90 mandelmod.o luna.o occultquad.o
./gen_pm