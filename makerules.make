# Use gfortran unless already defined
F90 ?= mpif90

ifeq ($(F90), mpif90)
	FFLAGS	?= -O2 -g -std=f2008 -Wall -Wextra
	ifeq ($(DEBUG), 1)
		FFLAGS += -fcheck=all -g -ffpe-trap=invalid,zero,overflow
	endif
else ifeq ($(F90), mpiifort)
	FFLAGS	:= -O2 -stand f08 -warn all
endif

%.o: 	%.f90
	$(F90) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

../%:	%.o
	$(F90) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))
