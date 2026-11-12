FC = gfortran
 
# make        -> version rapide, pour les vrais runs
# make debug  -> version de verification (5 a 20x plus lente)
FFLAGS_OPT   = -O3 -march=native
FFLAGS_DEBUG = -g -O0 -Wall -Wextra -Wno-tabs -fcheck=all -fbacktrace \
               -ffpe-trap=invalid,zero,overflow \
               -finit-real=snan -finit-integer=-99999999
 
FFLAGS ?= $(FFLAGS_OPT)
 
SRC = func_likelihood/sort_functions.f90 \
      func_likelihood/loglikelihood_errf_bin.f90
OBJ = sort_functions.o loglikelihood_errf_bin.o
 
all: run
 
debug:
	$(MAKE) clean
	$(MAKE) run FFLAGS="$(FFLAGS_DEBUG)"
 
run: libroutines.a b-Bayesian_binned_McMC.f90
	$(FC) $(FFLAGS) -o run b-Bayesian_binned_McMC.f90 -L./ -lroutines -lm
 
libroutines.a: $(SRC)
	for f in $(SRC); do $(FC) $(FFLAGS) -c $$f || exit 1; done
	ar -r libroutines.a $(OBJ)
	rm -f $(OBJ)
 
clean:
	rm -f *.a *.o *.mod run
 
.PHONY: all debug clean

	
