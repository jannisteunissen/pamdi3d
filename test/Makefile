FC 	:= gfortran
FFLAGS	:= -Wall -ffpe-trap=invalid,zero,overflow -fcheck=all -g -O2 -fopenmp
TESTS	:= test_m_particle_core test_m_cross_sec

# Need these libraries for the tests
BASEDIR	:= ../../
LIBS	:= particle_core fosito
INCDIRS := $(addprefix $(BASEDIR), $(LIBS))
LIBDIRS := ../ $(INCDIRS)

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all clean

all: 	$(TESTS)

test:	$(TESTS)
	$(foreach test, $(TESTS), ./$(test);)

clean:
	$(RM) *.o *.mod $(TESTS)
