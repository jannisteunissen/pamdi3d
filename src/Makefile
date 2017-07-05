FC := gfortran
FFLAGS := -Wall -O2 -fopenmp -Wno-unused-dummy-argument
VPATH := kdtree2:lookup_table_fortran:rng_fortran
OBJS := m_units_constants.o m_gas.o m_cross_sec.o kdtree2.o m_lookup_table.o	\
m_linked_list.o m_find_index.o m_random.o m_mrgrnk.o m_particle_core.o		\
m_photoi.o m_pc_all.o #m_particle_par.o

ifeq ($(DEBUG), 1)
	FFLAGS += -fcheck=array-temps,bounds,do,mem,pointer\
	-g -ffpe-trap=invalid,zero,overflow
endif

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

.PHONY: all clean

all: 	libparticle_core.a

libparticle_core.a: $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

test: 	libparticle_core.a
	$(MAKE) -C test

clean:
	$(RM) *.o *.mod libparticle_core.a

# Dependency information
m_particle_core.o:	m_cross_sec.o kdtree2.o m_lookup_table.o \
			m_linked_list.o m_random.o m_mrgrnk.o
# m_particle_par.o:	m_particle_core.o
m_gas.o:		m_units_constants.o
m_lookup_table.o:	m_find_index.o
m_photoi.o: m_random.o m_particle_core.o m_units_constants.o
m_pc_all.o:		m_particle_core.o m_gas.o m_photoi.o
