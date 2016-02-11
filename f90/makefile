# main makefile
# arch dependent flags
#
include make.defs


#F90 parallel version
F90_SRC =  tools.f90 ewald.f90 main.f90 reptation.f90
F90_OBJS =  tools.o ewald.o main.o reptation.o 

all : $(EXE) 
 

#
#  F90
#
$(EXE):$(F90_OBJS)
	$(F90) $(F90FLAGS) -o $@ $(F90_OBJS) $(LDFLAGS)

%.o:%.f90
	$(F90) $(FPP) $(F90FLAGS) -c $< 

#
#

clean:
	-rm *.o *.mod $(EXE) 


