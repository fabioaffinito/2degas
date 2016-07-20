# main makefile
# arch dependent flags
#
include make.defs


#F90 parallel version
SRCDIR=src
OBJDIR=obj
F90_SRC =  tools.f90 ewald.f90 utils.f90 input.f90 main.f90 reptation.f90 three_body.f90 oldblas.f90 algorithms.f90 mpif_mod.f90 mpif.f90  
OBJS =  tools.o ewald.o utils.o input.o main.o reptation.o three_body.o oldblas.o algorithms.o mpif_mod.o mpif.o 

F90_OBJS = $(patsubst %,$(OBJDIR)/%,$(OBJS))

all : $(EXE) 
 

#
#  F90
#
$(EXE):$(F90_OBJS)
	$(F90) $(F90FLAGS) -o $@ $(F90_OBJS) $(LDFLAGS)

$(OBJDIR)/%.o:$(SRCDIR)/%.f90
	$(F90) $(FPP) $(F90FLAGS) -c -o $@ $< 

#
#

clean:
	-rm $(F90_OBJS) *.mod $(EXE) 


