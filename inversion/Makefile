FC = gfortran
 
# make        -> version rapide, pour les vrais runs
# make debug  -> version de verification (5 a 20x plus lente)
PROG = b-Bayesian_binned
 
FFLAGS_OPT   = -O3 -march=native
FFLAGS_DEBUG = -g -O0 -Wall -Wextra -Wno-tabs -fcheck=all -fbacktrace \
               -ffpe-trap=invalid,zero,overflow \
               -finit-real=snan -finit-integer=-99999999
 
FFLAGS ?= $(FFLAGS_OPT)
 
SRC = src_likelihood/sort_functions.f90 \
      src_likelihood/likelihoods_FMD_gridsearch.f90 \
      src_likelihood/likelihoods_FMD_random_complete.f90
OBJ = sort_functions.o likelihoods_FMD_gridsearch.o likelihoods_FMD_random_complete.o
 
all: run
 
debug:
	$(MAKE) clean
	$(MAKE) run FFLAGS="$(FFLAGS_DEBUG)"
 
run: libroutines.a $(PROG).f90
	$(FC) $(FFLAGS) -o run $(PROG).f90 -L./ -lroutines -lm
 
libroutines.a: $(SRC)
	for f in $(SRC); do $(FC) $(FFLAGS) -c $$f || exit 1; done
	ar -r libroutines.a $(OBJ)
	rm -f $(OBJ)
 
clean:
	rm -f *.a *.o *.mod run
 
.PHONY: all debug clean

