SRC_DIRS	:= src

# Directories with altered names (useful for cleaning)
CLEANSRC	:= $(SRC_DIRS:%=clean-%)

.PHONY:	all clean $(SRC_DIRS) $(CLEANSRC)

all:	$(SRC_DIRS)

clean:	$(CLEANSRC)

$(SRC_DIRS):
	@$(MAKE) -C $@

$(CLEANSRC):
	$(MAKE) -C $(@:clean-%=%) clean

# Dependency information
tests: src
