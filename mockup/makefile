#
# make file for mockup
#

FC=ifort
SRCDIR=../src
SRC=ewald.f90 utils.f90 input.f90 mockup.f90 sonaseppia.f90
OBJ=ewald.o utils.o input.o mockup.o sonaseppia.o
FLAGS=-O0 -g


mockup:$(OBJ)
	$(FC) $(FLAGS) -qopenmp -o $@ $(OBJ)

mockup.o:mockup.f90
	 $(FC) $(FLAGS) -c -qopenmp -o $@ $<

sonaseppia.o:sonaseppia.f90
	 $(FC) $(FLAGS) -c -qopenmp -o $@ $<

%.o:$(SRCDIR)/%.f90
	$(FC) $(FLAGS) -c -qopenmp -o $@ $<

clean:
	rm -f *.o *.mod mockup

