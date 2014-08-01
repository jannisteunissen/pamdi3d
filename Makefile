FC 	:= mpif90
FFLAGS	:= -Wall -ffpe-trap=invalid,zero,overflow -fcheck=all -g -O2 -fopenmp
OBJS	:= m_gas.o m_particle_core.o m_particle_par.o m_cross_sec.o kdtree2.o

INCDIRS	:= ../fosito

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

.PHONY: all test clean

all: 	libparticle_core.a

libparticle_core.a: $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^


test: 	libparticle_core.a
	$(MAKE) -C test

clean:
	$(RM) -f *.o *.mod
	$(MAKE) -C test clean

# Dependency information
m_particle_core.o:	m_cross_sec.o kdtree2.o
m_particle_par.o:	m_particle_core.o
