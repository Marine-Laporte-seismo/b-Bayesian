F90 =  gfortran -fbounds-check # change above if needs some checks

all  :  obj global 

	
global: 
	$(F90) -o run b-Bayesian.f90 \
	-L./ -lroutines -lm
obj:
	$(F90) -c src/sort_functions.f90
	$(F90) -c src/likelihoods_random_marginals.f90
	$(F90) -c src/likelihoods_full.f90
	
	ar -r libroutines.a sort_functions.o likelihoods_random_marginals.o likelihoods_full.o
	\rm ./*.o
clean:
	/bin/rm *.a
	
