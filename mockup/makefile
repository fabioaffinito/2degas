#
# make file for mockup
#

FC=ifort
SRC=ewald.f90 utils.f90 input.f90 mockup.f90
OBJ=ewald.o utils.o input.o mockup.o



mockup:$(OBJ)
	$(FC) -qopenmp -o $@ $(OBJ)

%.o:%.f90
	$(FC) -c -qopenmp $<

.PHONY:clean
	rm -f *.o *.mod

