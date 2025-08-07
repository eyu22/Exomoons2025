# Clean first
rm *.o *.mod test_mandel

# Compile old Fortran code
gfortran -c occultquad.f

# Compile modern Fortran module
gfortran -c mandelmod.f90

# Compile and link the main test file
gfortran test_mandel.f90 mandelmod.o occultquad.o -o test_mandel
