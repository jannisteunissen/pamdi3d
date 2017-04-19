FFLAGS	:= -O2 -Wall -Wextra -fopenmp
FC 	:= gfortran

all:	kdtree2_test

kdtree2_test.o: kdtree2.o
kdtree2_test: kdtree2.o

clean:
	$(RM) *.o *.mod kdtree2_test

%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS)

%: %.o
	$(FC) -o $@ $^ $(FFLAGS)
