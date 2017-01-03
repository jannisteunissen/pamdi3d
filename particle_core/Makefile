OBJS	:= m_gas.o m_particle_core.o m_particle_par.o m_cross_sec.o kdtree2.o

INCDIRS	:= ../fosito

include ../makerules.make

.PHONY: all test clean

all: 	libparticle_core.a

libparticle_core.a: $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^


test: 	libparticle_core.a
	$(MAKE) -C test

clean:
	$(RM) -f *.o *.mod libparticle_core.a
	$(MAKE) -C test clean

# Dependency information
m_particle_core.o:	m_cross_sec.o kdtree2.o
m_particle_par.o:	m_particle_core.o
