SRC_DIRS	:= src fosito particle_core
CREATE_DIRS	:= output

# Directories with altered names (useful for cleaning)
CLEANSRC	:= $(SRC_DIRS:%=clean-%)

.PHONY:	all clean $(SRC_DIRS) $(CLEANSRC)

all: 		$(SRC_DIRS) | $(CREATE_DIRS)

clean: 		$(CLEANSRC)

$(SRC_DIRS):
		@echo -e "\n*********** Build information ***********"
		@echo -e "  Debug is set to: [$(DEBUG)],"
		@echo -e "  Set it to 1 to enable a debug build."
		@echo -e "  For example: make clean; make DEBUG=1"
		@echo -e "*****************************************\n"
		$(MAKE) -C $@
$(CREATE_DIRS):
		mkdir -p $@
$(CLEANSRC):
		$(MAKE) -C $(@:clean-%=%) clean

pamdi3d_libs/ready:
	make -C pamdi3d_libs

# Dependecy information
$(SRC_DIRS):	| pamdi3d_libs/ready
src:		fosito particle_core
particle_core: 	fosito
