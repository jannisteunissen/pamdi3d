FC 	:= gfortran
FFLAGS	:= -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow -g -O3
OBJS	:= m_gas.o m_particle_core.o m_cross_sec.o kdtree2.o
TESTS	:= test_m_particle_core test_m_cross_sec

INCDIRS	:= ../fosito
LIBDIRS := . ../fosito
LIBS	:= part_core fosito

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

.PHONY: all test clean

all: 	libpart_core.a

libpart_core.a: $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

$(TESTS): libpart_core.a

test: 	$(TESTS)
	$(foreach test, $(TESTS), ./$(test);)

clean:
	$(RM) -f *.o *.mod $(TESTS)

# Dependency information
m_particle_core.o:	m_cross_sec.o kdtree2.o
