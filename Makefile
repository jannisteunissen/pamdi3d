SRC_DIRS	:= src fosito particle_core
CREATE_DIRS	:= output

# Directories with altered names (useful for cleaning)
CLEANSRC	:= $(SRC_DIRS:%=clean-%)

.PHONY:	all clean $(SRC_DIRS) $(CLEANSRC) $(SUBMODULES)

all: 		$(SRC_DIRS) | $(CREATE_DIRS)

clean: 		$(CLEANSRC)

$(SRC_DIRS):
		$(MAKE) -C $@
$(CREATE_DIRS):
		mkdir -p $@
$(CLEANSRC):
		$(MAKE) -C $(@:clean-%=%) clean

pamdi3d_libs:
	tar -xzf pamdi3d_libs.tar.gz
	make -C pamdi3d_libs all

# Dependecy information
$(SRC_DIRS):	| pamdi3d_libs
src:		fosito particle_core
particle_core: 	fosito
