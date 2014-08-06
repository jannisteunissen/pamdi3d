FC 	:= mpif90
FFLAGS	:= -Wall -O2 #-fcheck=all -ffpe-trap=invalid,zero,overflow -g
OBJS	:= m_particle.o m_efield_amr.o m_init_cond.o m_electrode.o\
	 m_write_silo.o m_misc_process.o m_phys_domain.o m_photoi.o

INCDIRS	:= ../../fosito ../../particle_core\
 ../pamdi3d_libs/silo_lib/include
LIBDIRS := ../../fosito ../../particle_core  ../pamdi3d_libs/fishpack4.1/lib\
 ../pamdi3d_libs/silo_lib/lib
LIBS	:= particle_core fosito fishpack silo

%.o: 	%.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS))

../%:	%.o
	$(FC) -o $@ $^ $(FFLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

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

