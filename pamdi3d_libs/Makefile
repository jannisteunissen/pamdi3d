include ../makerules.make

SILO_DIR := silo_lib

.PHONY: all clean fishpack4.1 kdtree2
all: $(SILO_DIR) fishpack4.1 kdtree2 ready

fishpack4.1:
	$(MAKE) -C fishpack4.1

kdtree2:
	$(MAKE) -C kdtree2/src-f90

$(SILO_DIR):
	@echo "Locally installing the Silo library"
	@./build_silo.sh

ready:
	touch ready

clean:
	$(MAKE) -C fishpack4.1 clean
	$(MAKE) -C kdtree2/src-f90 clean
	$(MAKE) -C $(SILO_DIR) clean
	rm -rf silo_lib
