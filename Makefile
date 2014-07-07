# ##############################################################
#
# Streamer simulation makefile. Please read the installation instructions in INSTALL.txt
#
# ##############################################################


# ###############################
# Define the objects, source and module, include and library directories,
# and set the libraries that have to be linked.
# ###############################

ODIR		= objects
SDIR		= src
MDIR		= modules
UDIR	= ./pampi3d_libs
INCDIRS	= $(UDIR)/silo_lib/include $(UDIR)/kdtree2/src-f90
LIBDIRS	= $(UDIR)/fishpack4.1/lib $(UDIR)/silo_lib/lib $(UDIR)/kdtree2/src-f90
LIBS	= fishpack silo kdtree2


# ###############################
# Define the compilers to be used
# ###############################

F90C			= mpif90
F90L			= mpif90

# ###############################
# Set up flags for the compilation process
# ###############################

# Flags for warnings
F90_FWARNING 		= -Wall -finit-integer=-9999999 -finit-real=nan\
			 -fcheck=all -ffpe-trap=invalid,zero,overflow -g
# Flags for optimizations
F90_FOPTIM		= -O2
# Flags for compatability related settings
F90_FCOMPAT		= -fdefault-real-8 -fdefault-double-8 -cpp
F77_FCOMPAT		= -fdefault-real-8 -fdefault-double-8 -std=legacy
# Flags used for debugging
F90_FDEBUG		= -g -pg -O0

# Now create the full flags from the partial flags
F90FLAGS	= $(F90_FWARNING) $(F90_FOPTIM) $(F90_FCOMPAT)
FFLAGS		= $(F90_FWARNING) $(F90_FOPTIM) $(F77_FCOMPAT)

# The debug flag can set by typing make 'debug=1'
debug = 0
ifeq ($(debug), 1)
	F90FLAGS	+= $(F90_FDEBUG)
	FFLAGS		+= $(F90_FDEBUG)
endif


# ###############################
# Commands to compile and link Fortran 90 code.
# ###############################

# Compile standard source files
CompileF90	= $(F90C) -c -o $@ $< $(F90FLAGS) $(addprefix -I,$(INCDIRS)) -J $(MDIR)

# Compile F77 source files
Compilef77	= $(F90C) -c -o $@ $< $(FFLAGS) $(addprefix -I,$(INCDIRS)) -J $(MDIR)

# Command to link objects into executable
LinkAll	= $(F90L) -o $@ $^ $(F90FLAGS) $(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))


# ###############################
# Implicit rules for compilation.
# ###############################

# F95/90 files should have the extension .f90
$(ODIR)/%.o:	$(SDIR)/%.f90
		$(CompileF90)

# F77 files should have the extension .f
$(ODIR)/%.o:	$(SDIR)/%.f
		$(Compilef77)

# ###############################
# The source and object files
# ###############################

_SRC_F90	= module_globals.f90 module_parameters.f90 module_kiss.f90\
		generalUtilities.f90 module_crossSections.f90 module_constants.f90 \
		module_particle.f90 module_photoionization.f90 module_initialCond.f90 \
		m_efield_amr.f90 module_config.f90 module_electrode.f90 \
		module_IO.f90 module_gas.f90 m_write_silo.f90

_SRC_F77	=

_OBJS		= $(_SRC_F90:.f90=.o) $(_SRC_F77:.f=.o)
OBJS		= $(patsubst %,$(ODIR)/%,$(_OBJS))

# ###############################
# The targets
# ###############################

.PHONY:		all cleanlibs clean tarball

all: | pampi3d_libs
	@echo
	@echo "~~~ BUILDING pampi3d  ~~~"
	@echo " Checking whether directories exist... "
	@test -d $(ODIR) || mkdir $(ODIR)
	@test -d $(MDIR) || mkdir $(MDIR)
	@test -d output/ || mkdir output/
	@test -d $(UDIR) || (echo "Error: user directory does not exists! This is needed for libraries."; exit 1)
	@test -d input || (echo "Error: input directory does not exists! This is needed for cross sections."; exit 1)
	@test -d $(SDIR) || (echo "Error: Cannot find the source files."; exit 1)
	make pampi3d
	@echo
	@echo "~~~ Compilation complete   ~~~"
	@echo

pampi3d_libs:
	tar -xzf pampi3d_libs.tar.gz
	make -C pampi3d_libs all

cleanlibs:
	rm -rf pampi3d_libs

clean:
	rm -f $(ODIR)/*.o $(MDIR)/*.mod pampi3d

tarball:
	tar -czf pampi3d.tar.gz $(patsubst %,$(SDIR)/%,$(_SRC_F90)) $(SDIR)/pampi3d.f90 Makefile config.txt Doxyfile \
	input/sigloNew.txt COPYING.txt INSTALL.txt

# ###############################
# Now list the dependency information so
# that make knows how to compile the program.
# ###############################

pampi3d: 	$(ODIR)/pampi3d.o $(OBJS)
		$(LinkAll)

$(ODIR)/pampi3d.o :	$(ODIR)/module_crossSections.o \
				$(ODIR)/module_initialCond.o \
				$(ODIR)/module_particle.o \
				$(ODIR)/m_efield_amr.o \
				$(ODIR)/module_parameters.o \
				$(ODIR)/module_config.o \
				$(ODIR)/module_globals.o \
				$(ODIR)/module_constants.o \
				$(ODIR)/module_IO.o \
				$(ODIR)/module_photoionization.o \
				$(ODIR)/module_electrode.o \
				$(ODIR)/module_kiss.o \

$(ODIR)/module_crossSections.o :$(ODIR)/module_globals.o \
				$(ODIR)/generalUtilities.o \
				$(ODIR)/module_config.o \
				$(ODIR)/module_gas.o \

$(ODIR)/module_gas.o:		$(ODIR)/module_config.o \
				$(ODIR)/module_constants.o

$(ODIR)/module_particle.o :	$(ODIR)/module_globals.o \
				$(ODIR)/generalUtilities.o \
				$(ODIR)/m_efield_amr.o \
				$(ODIR)/module_constants.o \
				$(ODIR)/module_kiss.o \
				$(ODIR)/module_config.o \
				$(ODIR)/module_crossSections.o \
				$(ODIR)/module_electrode.o \

$(ODIR)/module_initialCond.o :	$(ODIR)/module_globals.o \
				$(ODIR)/generalUtilities.o \
				$(ODIR)/module_particle.o \
				$(ODIR)/module_kiss.o \
				$(ODIR)/m_efield_amr.o \
				$(ODIR)/module_config.o \
				$(ODIR)/module_constants.o \
				$(ODIR)/module_electrode.o \

$(ODIR)/m_efield_amr.o :	$(ODIR)/module_globals.o \
				$(ODIR)/m_write_silo.o \
				$(ODIR)/generalUtilities.o \
				$(ODIR)/module_electrode.o \
				$(ODIR)/module_constants.o \
				$(ODIR)/module_config.o \
				$(ODIR)/module_gas.o \

$(ODIR)/module_parameters.o :	$(ODIR)/module_globals.o \
				$(ODIR)/module_config.o \
				$(ODIR)/module_constants.o \

$(ODIR)/module_config.o :

$(ODIR)/module_globals.o :	$(ODIR)/module_constants.o \
				$(ODIR)/module_config.o \
				$(ODIR)/generalUtilities.o \

$(ODIR)/generalUtilities.o :	$(ODIR)/module_constants.o

$(ODIR)/module_electrode.o :	$(ODIR)/module_config.o \
				$(ODIR)/module_constants.o \
				$(ODIR)/generalUtilities.o \
				$(ODIR)/module_kiss.o \

$(ODIR)/module_photoionization.o : 	$(ODIR)/module_globals.o \
				$(ODIR)/module_kiss.o \
				$(ODIR)/generalUtilities.o \
				$(ODIR)/module_particle.o \
				$(ODIR)/m_efield_amr.o \
				$(ODIR)/module_config.o \
				$(ODIR)/module_constants.o \

$(ODIR)/module_IO.o :		$(ODIR)/module_particle.o \
				$(ODIR)/m_efield_amr.o \
				$(ODIR)/module_config.o \
				$(ODIR)/module_globals.o \
				$(ODIR)/m_write_silo.o \
