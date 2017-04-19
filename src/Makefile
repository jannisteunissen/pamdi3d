OBJS	:= m_particle.o m_efield_amr.o m_init_cond.o m_electrode.o\
	 m_write_silo.o m_misc_process.o m_phys_domain.o m_photoi.o

INCDIRS := ../fosito ../particle_core/src ../pamdi3d_libs/silo_lib/include
LIBDIRS := ../fosito ../particle_core/src ../pamdi3d_libs/fishpack4.1/lib	\
 ../pamdi3d_libs/silo_lib/lib
LIBS	:= particle_core fosito fishpack silo

include ../makerules.make

.PHONY: all clean

all: ../pamdi3d

clean:
	$(RM) *.o *.mod ../pamdi3d

# Dependency information
../pamdi3d:		$(OBJS)
m_particle.o:		m_efield_amr.o m_photoi.o
m_efield_amr.o:		m_electrode.o m_write_silo.o m_phys_domain.o
m_init_cond.o:		m_efield_amr.o
m_misc_process.o:	m_particle.o
m_phys_domain.o:	m_electrode.o
