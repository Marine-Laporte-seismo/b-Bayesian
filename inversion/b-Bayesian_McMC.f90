program b_Bayesian_binned_McMC_detect
  implicit none
  real , EXTERNAL    ::   ran3, gasdev
  
  !            ----------------------------------
  !              # BEGINING OF THE USER-DEFINED PARAMETERS
  !            ----------------------------------
!-----------------------------------------
!-----------------------------------------
! # Data 
!----------------------------------------- 
  character(LEN=*), parameter :: Magnitude_file = '../inputs/TB_6DISC_UNCLEAR_2710_ERROR_ETAS_erf.in' !path of the time-magnitude catalog
  integer, parameter :: N_data_max = 800000		! Number maximum of data
  real*8, parameter :: delta_bin = 0.01	  		! Size of the Discrete Time Vector
  integer, parameter :: T_bins = 200		    	! Binning of the catalog	
  real*8, parameter :: min_bin = 0		      	! Min bin < min magnitudes for integration lower bound
  real*8, parameter :: max_bin = 7.0		      	! Max bin > max magnitudes for integration upper bound
  
!-----------------------------------------
! # Parameters of the Markov chain
!-----------------------------------------
  integer, parameter :: it_max = 10000000	    	! Total number of McMC iterations
  integer, parameter :: it_burnin = 4000000		! Burn-in iterations
  integer, parameter :: it_thin = 100	      		! Thinning for independency of accepted models
  integer, parameter :: it_std = 1000000			! Number of iterations between two perturbation standard-deviation changes
  integer, parameter :: size_thin = (it_max/it_thin)	! For storing models every thinned iterations -> less heavy outputs
  integer, parameter :: size_std = (it_max/it_std)	! For changing the std of all parameters -> accelerating convergence
  
  integer, parameter :: window_AR = 10000  		! for monitoring acceptance rate every window_AR iterations

!-----------------------------------------
! # Definition of priors
!-----------------------------------------
  real*8, parameter :: beta_min = 0.5d0*2.3d0, beta_max = 2.0d0*2.3d0     ! Bounds of the uniform prior of Beta (b_value*ln(10))
  real*8, parameter :: mu_min = 0.0d0, mu_max = 4.0d0                     ! Bounds of the Mu uniform prior
  real*8, parameter :: sigma_min = 0.0001d0, sigma_max = 0.5d0            ! Bounds of the Sigma uniform prior

  !-----------------------------------------
  ! # Definition de l'initialisation GR
  !-----------------------------------------  
  real*8, parameter :: beta_init = 1.0d0*2.3d0, beta_stdi = 0.25d0, beta_std_birth = 0.5d0
  real*8, parameter :: mu_init = 0.5d0, mu_stdi = 0.15d0, mu_std_birth = 0.3d0
  real*8, parameter :: sigma_init = 0.3d0, sigma_stdi = 0.04d0, sigma_std_birth = 0.1d0

  real*8, parameter :: beta_store = 0.025d0, mu_store = 0.025d0, sigma_store = 0.01d0

!-----------------------------------------

!-----------------------------------------
  ! # Definition of initialization T
  !-----------------------------------------	
  integer, parameter :: max_layers = 50		! Threshold on maximum layers	
  integer, parameter :: Nc_start_min = 10		! Initial number of layers min
  integer, parameter :: Nc_start_max = 20		! Initial number of layers max
  real*8, parameter :: Disc_stdi = 1.5			! Initial standard deviation for layer perturbation 
 
  !            ----------------------------------
  !              # END OF THE USER-DEFINED PARAMETERS
  !            ----------------------------------
	


!-----------------------------------------
! # Functions declaration
!-----------------------------------------			
interface  
    ! A function that computes the loglikelihood of the model ("beta", "mu", "sigma") for a subset of data
    ! Considering a GR law described by Beta modulated by a detection function described by Mu/Sigma
    ! No binning 
    subroutine LogLikelihood_function_nobin_detect(mag_obs, N_data_subset, min_mags, &
    						 New_Beta, New_Mu, New_Sigma, LogLikelihood_prop)
  	real*8, intent(in) :: mag_obs(:)
  	integer, intent(in) :: N_data_subset
  	real*8, intent(in) :: min_mags
  	real*8, intent(in) :: New_Beta
  	real*8, intent(in) :: New_Mu
  	real*8, intent(in) :: New_Sigma
  	real*8, intent(out) :: LogLikelihood_prop
    end subroutine LogLikelihood_function_nobin_detect
    
    subroutine sort_shell(arr)						! Fonction pour trier des listes d'entiers
  	integer, DIMENSION(:), INTENT(INOUT) :: arr
	integer :: i,j,inc,n
	REAL*8 :: v
    end subroutine sort_shell
    subroutine sort_shell_real(arr)					! Fonction pour trier des listes de reels
  	real*8, DIMENSION(:), INTENT(INOUT) :: arr
	integer :: i,j,inc,n
	REAL*8 :: v
    end subroutine sort_shell_real 
end interface  
   
	!-----------------------------------------
	! # DECLARATIONS OF VARIABLES 
	!-----------------------------------------
integer :: ra							! seed of our run (change if we want a different run)

!-----------------------------------------
! # DECLARATIONS OF VARIABLES DATA
!-----------------------------------------
real*8 :: mag_obs0(N_data_max), times0(N_data_max)   	! max data	
real*8, allocatable :: mag_obs(:), times(:) 		! allocated data
integer :: N_data_mag, N_data_subset			! number of data in the full dataset / by subsets
real*8 :: time_lim, time_max, time_min			! min/max times, bounds for random draw
real*8 :: max_mags, min_mags				! min/max magnitudes in the full dataset
! magnitude bins
integer :: bin_index, N_bins

integer :: Nc, New_Nc				! nombre de couches 
!-----------------------------------------
! # DECLARATIONS OF VARIABLES BMUSIGMA
!-----------------------------------------
! Beta 
integer :: beta_bins ! number of bins in the beta grid 
integer :: bin_b     ! bin in the beta grid for a beta value 
real*8 :: Acc_Beta   ! accepted beta
real*8 :: New_Beta   ! proposed beta
real*8 :: beta_std   ! for the beta moves 
real*8 :: beta_delta ! for the birth/death term of acceptance alpha


! Mu
integer :: mu_bins    ! number of bins in the mu grid 
integer :: bin_mu     ! bin in the mu grid for a mu value 
real*8 :: Acc_Mu      ! accepted mu
real*8 :: New_Mu      ! proposed mu
real*8 :: mu_std      ! for the mu moves 
real*8 :: mu_delta    ! for the birth/death term of acceptance alpha

! Sigma
integer :: sigma_bins 	! number of bins in the sigma grid 
integer :: bin_sigma    ! bin in the sigma grid for a sigma value 
real*8 :: Acc_Sigma  	! accepted sigma
real*8 :: New_Sigma   	! proposed sigma
real*8 :: sigma_std   	! for the sigma moves 
real*8 :: sigma_delta 	! for the birth/death term of acceptance alpha


real*8 :: Acc_Loglikelihood, loglikelihood, New_Loglikelihood 
real*8 :: Acc_loglikelihood_layers(max_layers+2),New_loglikelihood_layers(max_layers+2)
real*8 :: Sum_LogLikelihood, Acc_Sum_Loglikelihood, New_Sum_Loglikelihood


!-----------------------------------------
! # DECLARATIONS OF VARIABLES TEMPORAL MODEL 
!-----------------------------------------
real*8 :: model_times(max_layers+1)	! vecteur temps des perturbations
integer :: model_ibins(max_layers+1)	! vecteur indices des perturbations dans le vecteur prop
integer :: model_idata(max_layers+1)	! vecteur indices des perturbations dans le vecteur data 

real*8 :: New_model_times(max_layers+1)		! [Proposed] Array of the times of discontinuities (different of catalog times)
integer :: New_model_ibins(max_layers+1)	! [Proposed] Array of indexed of these times in the T_bins vector
integer :: New_model_idata(max_layers+1)	! [Proposed] Array of indexed of these times in the data 

integer :: a_data, b_data, a_ibins, b_ibins	! Indices pour decoupe de catalogues entre deux discontinuites temporelles
integer ::  index_data, index_ibins_old, index_ibins

integer :: index_pert					! Indice de la discontinuite a perturber
real*8 :: t_old, t_new, t_init, t_std				! Temps de discontinuite uen fois perturbe selon tstd, temps minimal des perturbations, standard deviation temporelle

integer :: model_T(T_bins), Sum_accepted_T(T_bins)		! Matrice des discontinuites acc/ Stack des discontinuites acc
integer :: New_model_T(T_bins)				! Matrice temporaire des discontinuites proposees

!-----------------------------------------
! DECLARATIONS OF VARIABLES GR MODEL
!-----------------------------------------
! Beta
real*8 :: model_beta(max_layers+2)		! vector of beta values for each segment of the current temporal model
real*8 :: New_model_beta(max_layers+2)		! vector of beta values for each segment of the proposed (TMP) temporal model
real*8, allocatable :: Beta_T_posterior(:,:)	! matrix of sum of accepted beta on a beta_grid for each temporal bins 

! Mu
real*8 :: model_mu(max_layers+2)		! vector of mu values for each segment of the current temporal model
real*8 :: New_model_mu(max_layers+2)		! vector of mu values for each segment of the proposed (TMP) temporal model
real*8, allocatable :: Mu_T_posterior(:,:)	! matrix of sum of accepted mu on a mu_grid for each temporal bins 

! Sigma
real*8 :: model_sigma(max_layers+2)		! vector of sigma values for each segment of the current temporal model
real*8 :: New_model_sigma(max_layers+2)		! vector of sigma values for each segment of the proposed (TMP) temporal model
real*8, allocatable :: Sigma_T_posterior(:,:)	! matrix of sum of accepted sigma on a sigma_grid for each temporal bins 

!-----------------------------------------
! Store each change at each iteration McMC/thinning
!-----------------------------------------

integer :: Nc_total(size_thin+1)	
real*8 :: Store_times_proposed(max_layers+1, size_std+1)
	
real*8 :: Store_Beta_proposed(max_layers+2, size_std+1)
real*8 :: Store_Mu_proposed(max_layers+2, size_std+1)
real*8 :: Store_Sigma_proposed(max_layers+2, size_std+1)

!-----------------------------------------
! # DECLARATIONS OF VARIABLES Acceptance Rates
!-----------------------------------------
integer :: N_GEO, N_GR						! Count number of change in geometry/gutenberg-richter
integer :: N_GEO_accept, N_GR_accept				! Count number of change in geometry/gutenberg-richter accepted
integer :: N_move, N_birth, N_death , store_death		! Count number of move/birth /death 
integer :: N_move_accept, N_birth_accept, N_death_accept	! Count number of move/birth /death accepted
integer :: N_beta, N_mu, N_sigma				! Count number of change of beta/mu/sigma
integer :: N_beta_accept, N_mu_accept, N_sigma_accept		! Count number of change of beta/mu/sigma accepted

real*8 ::  AR_GEO, AR_GR			! Acceptance Rate geometry/gutenberg-richter
real*8 ::  AR_move, AR_birth, AR_death		! Acceptance Rate move/birth/death
real*8 ::  AR_beta, AR_mu, AR_sigma		! Acceptance Rate beta/mu/sigma

integer :: flag_prior_beta
integer :: flag_prior_mu
integer :: flag_prior_sigma
integer ::  c, flag, flag_std


integer            :: n_seg, i_seg
! Compteurs par segment (remis à zéro à chaque nouveau segment)
integer :: seg_N_GEO,   seg_N_GEO_accept
integer :: seg_N_move,  seg_N_move_accept
integer :: seg_N_death, seg_N_death_accept
integer :: seg_N_birth, seg_N_birth_accept
integer :: seg_N_GR,    seg_N_GR_accept
integer :: seg_N_beta,  seg_N_beta_accept
integer :: seg_N_mu,    seg_N_mu_accept
integer :: seg_N_sigma, seg_N_sigma_accept
!-----------------------------------------
! # DECLARATIONS OF VARIABLES Metropolis Hastings
!-----------------------------------------
integer ::  flag_unaccepted, flag_ndata(max_layers+1)			! some models are ill constructed - reject directly
real*8 :: alpha_A, alpha_T, alpha_GR					! conditions A : geometry/GR ; T: birth/death/pert ; GR : beta/mu/sigma
real*8 :: u, condition, diff_log, term_added							! condition acceptation Metropolis Hastings
real*8 :: terme_beta_1, terme_beta_2, terme_mu_1, terme_mu_2, terme_sigma_1, terme_sigma_2
real*8:: New_std						! Methode T.Santos : adaptation du std en fonction de Acc_pert
integer :: N_thin, flag_thin						! Pour independance des parametres : garder les modeles acceptes tous les dix iterations
!-----------------------------------------
! # If Init == True
!-----------------------------------------
real*8 :: t_true(2)=(/40,60/)
real*8 :: b_true(3)=(/0.8*2.3,1.1*2.3,0.8*2.3/)
!-----------------------------------------
! # Others
!-----------------------------------------
integer :: io, j, i, k, it
integer :: line, column
real*8 :: t1, t2, comp_time, t3, t4, progress
real*8 :: t_it(it_max),mean_t_it

real*8 :: t_offset                       ! origine des temps soustraite au catalogue
real*8 :: AR_move_seg                    ! taux d'acceptation des moves sur la dernière fenêtre
integer :: j_start                       ! premier bin temporel d'un segment

real*8, parameter :: pi = acos(-1.0d0)
logical, parameter :: prior_only = .false.


	!-----------------------------------------
	! # END OF VARIABLES DECLARATION 
	!-----------------------------------------









	!-----------------------------------------
	! # VARIABLES INIT
	!-----------------------------------------
!open(unit=90,file='seed.txt',status='old')
!read(90,*) ra
!close(90)
ra=14	! seed of our run (change if we want a different run)

model_times =  9999999
model_ibins =  9999999
model_idata =  9999999

Nc_total = 0
Store_times_proposed = 99
Store_Beta_proposed = 99
Store_Mu_proposed = 99
Store_Sigma_proposed = 99

Acc_Loglikelihood = 0
loglikelihood = 0
sum_logLikelihood =0 
New_model_T =0
model_T=0
Sum_accepted_T=0

model_beta = 99
New_model_beta = 99
model_mu = 99
New_model_mu = 99
model_sigma = 99
New_model_sigma = 99

N_thin=1
n_seg = it_max / window_AR
N_GEO = 0
N_move = 0
N_death = 0
N_birth = 0
N_GR = 0
N_beta = 0
N_mu = 0
N_sigma = 0

AR_GEO = 0; AR_move = 0; AR_death = 0; AR_birth = 0
AR_GR = 0;  AR_beta = 0; AR_mu = 0;    AR_sigma = 0
AR_move_seg = 40.0d0
condition = 0; diff_log = 0; term_added = 0
New_Sum_Loglikelihood = 0
Acc_loglikelihood_layers = 0.0d0
New_loglikelihood_layers = 0.0d0

seg_N_GEO=0;   seg_N_GEO_accept=0
seg_N_move=0;  seg_N_move_accept=0
seg_N_death=0; seg_N_death_accept=0
seg_N_birth=0; seg_N_birth_accept=0
seg_N_GR=0;    seg_N_GR_accept=0
seg_N_beta=0;  seg_N_beta_accept=0
seg_N_mu=0;    seg_N_mu_accept=0
seg_N_sigma=0; seg_N_sigma_accept=0





!-----------------------------------------
! Etape 0.B : Read the magnitude data
! The magnitude catalog must only have two colums
! C1 : Time between 0 and TMAX in minuts or hour or day 
! C2 : Event magnitude
!-----------------------------------------
open(50,file= Magnitude_file ,status='old')
    N_data_mag = 0
    do line =1, N_data_max
	      read(50,*,IOSTAT=io) times0(line),mag_obs0(line)	
	      if (io > 0) then				
			    stop "Check input.  Something was wrong"
	      elseif (io < 0) then
			    exit
	      else
	      	N_data_mag = N_data_mag+1
	      endif
	      if (line == N_data_max) stop "number of magnitudes data >= N_data_magmax"
    enddo
 close(50)! close the file

!-----------------------------------------
! allocate de tous les tableaux qui dependent de ndata 
allocate(mag_obs(N_data_mag))
allocate(times(N_data_mag))
mag_obs = mag_obs0(1:N_data_mag)
times =  times0(1:N_data_mag)

do i = 2, N_data_mag
   if (times(i) < times(i-1)) stop "Catalogue non trie en temps"
end do

write(*,*) 'Size of the catalog :', N_data_mag
write(*,*) 'Minimum magnitude :', minval(mag_obs)
write(*,*) 'Maximum magnitude :', maxval(mag_obs)
open(unit=61, file="../outputs/input.txt",  status="replace", action="write")
	do i = 1, N_data_mag
		write(61, *) times(i),mag_obs(i)
	end do
close(61)
 
!-----------------------------------------
! Etape 0.b : Create the list of mag_bins (bins =! uncertainty)2828
!-----------------------------------------
min_mags = minval(mag_obs) 
max_mags = maxval(mag_obs) 
!-----------------------------------------
!-----------------------------------------
! Etape 0.c : Initialisation time variables
!-----------------------------------------
t_offset = minval(times)
times = times - t_offset          ! le temps commence maintenant à 0
time_max = maxval(times)
time_min = minval(times)          ! = 0


if (time_max .EQ. time_min ) then 
	if (time_min .NE. 0) then
		time_max = 2 * time_min 
	else 
		time_max = 1
	end if
else if (time_max == 0) then
	time_max = 1
end if

time_lim =  0		
t_std = Disc_stdi			! Disc stdi correspond au vrai std


!-----------------------------------------
! Etape 0.d : Create the list of GR-variables bins
!-----------------------------------------

beta_bins = nint((beta_max - beta_min)/(beta_store)) + 1
mu_bins = ((mu_max - mu_min)/(mu_store))+1
sigma_bins = nint((sigma_max - sigma_min)/(sigma_store)) + 1
write(*,*) "Beta Bins", beta_bins

allocate(Beta_T_posterior(beta_bins, T_bins))
allocate(Mu_T_posterior(mu_bins, T_bins))
allocate(Sigma_T_posterior(sigma_bins, T_bins))
Beta_T_posterior = 0
Mu_T_posterior = 0
Sigma_T_posterior = 0









!------------------------------------------------------------------------
! 1- INITIALIZATION OF THE TRANSDIMENSIONAL INVERSION
!-----------------------------------------------------------------------

  CALL cpu_time(t1)  !Tic. start counting time 
	
  !!1.1 Random draw of the number and location of temporal discontinuities !!
  !-----------------------------------------
  write(*,*) "> Begin __init__... "
  write(*,*) "------------------- "
  write(*,*) "-- init temporal discontinuities :"
  ! Random draw of the number of discontinuities
  Nc = Nc_start_min + int(ran3(ra)*(Nc_start_max-Nc_start_min+1))
  Nc = min(Nc, Nc_start_max)
  write(*,*) ". Initial number of discontinuities :", Nc
   
 
  do k=1, Nc
  	! Random draw of discontinuity time (not necessarily a catalog times = real between time_min time_max)
  	!t_init = t_true(k)
	t_init = time_lim + (ran3(ra)*(time_max-(2*time_lim)))

	! a - we find the index associated with this random observation time in the T_bins vector
    	index_ibins = ceiling(T_bins * t_init/(time_max))
    	! b - we find the index associated with this random observation time in the data 	
    	index_data = 0
    	do i = 1, N_data_mag 
    		if (times(i) .le. t_init) then 	
    			index_data = i
    		end if 
    	end do				
    	model_times(k) = t_init			! store the discontinuities times				
    	model_ibins(k) = index_ibins		! store its index in the T_bins vector
    	model_idata(k) = index_data		! store its index in the data vector
    	model_T(index_ibins) = model_T(index_ibins) + 1		! assigns 1 to the location of a discontinuity [0,0,0,1,0,...] in the vector model_T(T_bins)
    	Store_times_proposed(k,1) = t_init	! store the discontinuities times at each iteration
    	Nc_total(1) = Nc			! store the number of layers accepted Nc at each iteration 
   end do
   call sort_shell_real(Store_times_proposed(:,1))
   call sort_shell_real(model_times)
   call sort_shell(model_ibins)
   call sort_shell(model_idata)
   

   write(*,*) ". Initial model times :", model_times(1:Nc)
   write(*,*) ". Initial model prop :", model_ibins(1:Nc)
   write(*,*) ". Initial model data :", model_idata(1:Nc)
  
  
  
  !! 1.2 Computing conditional P(W|T,D) for each segment !! 
  !------------------------------------------------------
   write(*,*) ". Now, computing first likelihood..."
   write(*,*)".. Likelihood :  id | begin | end | n_data | Beta | Mu | Sigma | LogLikelihood Cumul_Log | "
   do k = 1, Nc+1
         if (k==1) then
    		flag_ndata(1) = model_idata(1) 
    		a_data = 1 
      		b_data = a_data + flag_ndata(1) - 1
      		
      		a_ibins = 1
      		b_ibins = model_ibins(1)	
    	
    	else if (k == Nc+1) then
    		flag_ndata(k) = N_data_mag - model_idata(k-1)  
    		a_data = a_data + flag_ndata(k-1)
	    	b_data = a_data + flag_ndata(k) - 1
    		
    		a_ibins = model_ibins(k-1)
      		b_ibins = T_bins
      		
    	else
    		flag_ndata(k) = model_idata(k) - model_idata(k-1) 
    		a_data = a_data + flag_ndata(k-1)
	    	b_data = a_data + flag_ndata(k) - 1
	    	
	      	a_ibins = model_ibins(k-1)
	      	b_ibins = model_ibins(k)
    	end if

    	
    	! a1 - random draw of initial beta
    	Acc_Beta = beta_init + gasdev(ra) * beta_stdi
  	!Acc_Beta = b_true(k)
	
	! beta init
    	if ((Acc_Beta.lt.beta_min) .or. (Acc_Beta.gt.beta_max)) then
        	Acc_Beta = beta_init + gasdev(ra) * beta_stdi
    		if ((Acc_Beta.lt.beta_min) .or. (Acc_Beta.gt.beta_max)) then
    			write(*,*) " ... New Beta random draw ! : problem in the hyperparameters  "
    			stop
    		end if
    	end if
    	! beta store
    	bin_b = int((Acc_Beta - beta_min)/beta_store)+1 
	model_beta(k) = Acc_Beta
	Store_Beta_proposed(k, 1) = Acc_Beta
	
	! a2 - random draw of initial mu
    	Acc_Mu = mu_init + gasdev(ra) * mu_stdi
	! mu init
    	if ((Acc_Mu.lt.mu_min) .or. (Acc_Mu.gt.mu_max)) then
       		Acc_Mu = mu_init + gasdev(ra) * mu_stdi
    		if ((Acc_Mu.lt.mu_min) .or. (Acc_Mu.gt.mu_max)) then
    			write(*,*) " ... New Mu random draw ! : problem in the hyperparameters  "
    			stop
    		end if
    	end if
    	! mustore
    	bin_mu = int((Acc_Mu - mu_min)/mu_store)+1 
	model_mu(k) = Acc_Mu
	Store_Mu_proposed(k, 1) = Acc_Mu
	
	
	! a3 - random draw of initial sigma
    	Acc_Sigma = sigma_init + gasdev(ra) * sigma_stdi
	! sigma init
    	if ((Acc_Sigma.lt.sigma_min) .or. (Acc_Sigma.gt.sigma_max)) then
      		Acc_Sigma = sigma_init + gasdev(ra) * sigma_stdi
    		if ((Acc_Sigma.lt.sigma_min) .or. (Acc_Sigma.gt.sigma_max)) then
    			write(*,*) " ... New Sigma random draw ! : problem in the hyperparameters  "
    			stop 	
    		end if
    	end if
    	! sigma store
    	bin_sigma = int((Acc_Sigma - sigma_min)/sigma_store)+1 
	model_sigma(k) = Acc_Sigma
	Store_Sigma_proposed(k, 1) = Acc_Sigma
	


	
	
	! b - compute loglikelihood if possible
    	if (flag_ndata(k)==0) then 
    		!write(*,*) '.....No events in partition'
    		Acc_Loglikelihood = 0
    	
    	else 	
		!write(*,*) ".....Events in partition :", flag_ndata(k)							! Au moins un evenement dans le segment								

		N_data_subset = size(mag_obs(a_data:b_data))



		call LogLikelihood_function_nobin_detect(mag_obs(a_data:b_data),N_data_subset, min_mags, & 
							Acc_Beta, Acc_Mu, Acc_Sigma, Acc_Loglikelihood)
	end if 	
	
	Acc_loglikelihood_layers(k) = Acc_Loglikelihood
	
	write(*,*) ".. Likelihood /segt",k,a_data,b_data,flag_ndata(k), Acc_Beta, Acc_Mu, Acc_Sigma, &
					Acc_Loglikelihood	
    end do 
   Sum_loglikelihood = sum(Acc_loglikelihood_layers(1:Nc+1))
   write(*,*) "------------------------------------------------"
   write(*,*) ". Initial model beta :", model_beta(1:Nc+1)
   write(*,*) ". Initial model mu :", model_mu(1:Nc+1)
   write(*,*) ". Initial model sigma :", model_sigma(1:Nc+1)
   write(*,*) ". Initial Log-likelihood :", Acc_loglikelihood_layers(1:Nc+1)
   write(*,*) "------------------------------------------------"	
    	
    	



  !!-------------------------------------------
  !! 1.3 First randow draw new model : 
  !!  MOVE/BIRTH/DEATH - BETA/MU/SIGMA!!
  !!------------------- ---------------------- 
  flag_unaccepted = 0 
  index_ibins = 0
  index_ibins_old = 0
  
  New_Nc = Nc
  
  New_model_T = model_T
  New_model_times =  model_times
  New_model_ibins =  model_ibins
  New_model_idata =  model_idata
  
  New_model_beta = model_beta
  
  
  
  

  flag_prior_beta = 0
  flag_thin = 1
  
  
   
   
!----------------------------------------------------------------------------------------------------
! 2- BEGIN TRANSDIMENSIONAL ITERATIONS
!------------------------------------------------------------------------------------------------------
  !write(*,*) ". Initial Log-likelihood :", Acc_loglikelihood_layers(1:Nc+1)
  write(*,*)
  write(*,*)
  write(*,*)
  write(*,*) "> Begin reversible-jump..."
  write(*,*) "--------------------------"
  write(*,*) ". Total number of iterations : ", it_max
  write(*,*) ". it_thin : ", it_std
  write(*,*)
  write(*,*) "> First iteration : 1"
  write(*,*) ". Computing the likelihood of the new model (for comparison)..."
  
  flag=0
  ! Delta GR
  beta_delta = beta_max - beta_min
  beta_std = beta_stdi
  
  mu_delta = mu_max - mu_min
  mu_std = mu_stdi
  
  sigma_delta = sigma_max - sigma_min
  sigma_std = sigma_stdi
  
  term_added = 0
  ! Acceptance counts
  N_GEO_accept = 0
  N_GR_accept = 0
  N_move_accept = 0
  N_death_accept = 0
  N_birth_accept = 0
  N_beta_accept = 0
  N_mu_accept = 0
  N_sigma_accept = 0
  
  Acc_Sum_Loglikelihood = Sum_loglikelihood
   ! case random exploration pour condition MCMC
  New_std = Disc_stdi
  
  
  open(unit=75, file="../outputs/AR_evolution.txt", status="replace", action="write")
  do it = 2, it_max
  	CALL cpu_time(t4) 
  	
  	if (mod(it, window_AR) == 0) then
	    i_seg = it / window_AR

	    ! Calcul AR segment (éviter division par zéro)
	    if (seg_N_GEO   > 0) AR_GEO   = real(seg_N_GEO_accept)   / real(seg_N_GEO)   * 100
	    if (seg_N_death > 0) AR_death = real(seg_N_death_accept) / real(seg_N_death) * 100
	    if (seg_N_birth > 0) AR_birth = real(seg_N_birth_accept) / real(seg_N_birth) * 100
	    if (seg_N_GR    > 0) AR_GR    = real(seg_N_GR_accept)    / real(seg_N_GR)    * 100
	    if (seg_N_beta  > 0) AR_beta  = real(seg_N_beta_accept)  / real(seg_N_beta)  * 100
	    if (seg_N_mu    > 0) AR_mu    = real(seg_N_mu_accept)    / real(seg_N_mu)    * 100
	    if (seg_N_sigma > 0) AR_sigma = real(seg_N_sigma_accept) / real(seg_N_sigma) * 100
	    
	    if (seg_N_move  > 0) AR_move  = real(seg_N_move_accept)  / real(seg_N_move)  * 100
	    if (seg_N_move  > 0) AR_move_seg = dble(seg_N_move_accept) / dble(seg_N_move) * 100.0d0

	    ! Sauvegarde dans le fichier d'évolution
	    write(75, '(I10, 8F10.3)') it, AR_GEO, AR_move, AR_death, AR_birth, &
		                            AR_GR, AR_beta, AR_mu, AR_sigma

	    ! Remise à zéro des compteurs du segment
	    seg_N_GEO=0;   seg_N_GEO_accept=0
	    seg_N_move=0;  seg_N_move_accept=0
	    seg_N_death=0; seg_N_death_accept=0
	    seg_N_birth=0; seg_N_birth_accept=0
	    seg_N_GR=0;    seg_N_GR_accept=0
	    seg_N_beta=0;  seg_N_beta_accept=0
	    seg_N_mu=0;    seg_N_mu_accept=0
	    seg_N_sigma=0; seg_N_sigma_accept=0

		  ! Petit bilan des taux d'acceptance 
		  AR_GEO = real(N_GEO_accept)/real(N_GEO)*100
		  AR_move = real(N_move_accept)/real(N_move)*100
		  AR_death = real(N_death_accept)/real(N_death)*100
		  AR_birth = real(N_birth_accept)/real(N_birth)*100
		  AR_GR = real(N_GR_accept)/real(N_GR)*100
		  AR_beta = real(N_beta_accept)/real(N_beta)*100
		  AR_mu = real(N_mu_accept)/real(N_mu)*100
		  AR_sigma = real(N_sigma_accept)/real(N_sigma)*100
		  write(*,*) ". Review of acceptance rates : "
		  write(*,*) ". Name : number_accepted  number_proposed    rate(pct)"
		  write(*,*) ".. # GEOMETRY :", N_GEO_accept, N_GEO, AR_GEO
		  write(*,*) ".. Births : ", N_birth_accept, N_birth, AR_birth
		  write(*,*) ".. Deaths : ", N_death_accept, N_death, AR_death
		  write(*,*) ".. Moves  : ", N_move_accept, N_move, AR_move
		  write(*,*) ".. # GR parameters :", N_GR_accept, N_GR, AR_GR
		  write(*,*) ".. Beta  : ", N_beta_accept, N_beta, AR_beta
		  write(*,*) ".. Mu  : ", N_mu_accept, N_mu, AR_mu
		  write(*,*) ".. Sigma  : ", N_sigma_accept, N_sigma, AR_sigma

		  

		if ((mod(it,it_std)==0).and.(it .le. it_burnin)) then
		    write(*,*) ". Old t_std :", t_std, "  AR_move fenetre (%) :", AR_move_seg
		    t_std = t_std * max(0.5d0, min(2.0d0, AR_move_seg/40.0d0))
		    t_std = max(t_std, 1.0d-3*(time_max-2*time_lim))
		    t_std = min(t_std, time_max-2*time_lim)
		    write(*,*) ". New t_std :", t_std
		end if
		  
		  
		  open(unit=62, file="../outputs/Acceptances.txt",  status="replace", action="write")
			  write(62,*) " Acceptance rates at iteration : ", it
			  write(62,*) " Name : number_accepted  number_proposed    rate(pct)"
			  write(62,*) " # GEOMETRY :", N_GEO_accept, N_GEO, AR_GEO
			  write(62,*) "Births : ", N_birth_accept, N_birth, AR_birth
			  write(62,*) "Deaths : ", N_death_accept, N_death, AR_death
			  write(62,*) "Moves  : ", N_move_accept, N_move, AR_move
			  write(62,*) " # GR parameters :", N_GR_accept, N_GR, AR_GR
			  write(62,*) "Beta  : ", N_beta_accept, N_beta, AR_beta
			  write(62,*) "Mu  : ", N_mu_accept, N_mu, AR_mu
			  write(62,*) "Sigma  : ", N_sigma_accept, N_sigma, AR_sigma
		close(62)
	end if

	
	! INITIALISATION DU PROPOSAL 
	index_ibins = 0
	index_ibins_old = 0

	New_Nc = Nc
	New_model_T = model_T
	New_model_times =  model_times   ! proposal = current
	New_model_ibins =  model_ibins   ! proposal = current
	New_model_idata =  model_idata   ! proposal = current

	New_model_beta = model_beta
	New_model_mu = model_mu
	New_model_sigma = model_sigma

	flag_unaccepted = 0
	flag_prior_beta = 0
	flag_prior_mu = 0
	flag_prior_sigma = 0
	
	
	! UPDATE DU PROPOSAL
	!write(*,*) ". New random draw : Geometry vs Gutenberg-Richter parameters"
	alpha_A = ran3(ra) 
	
	
	if (alpha_A.lt.1./2) then 							! Case we perturb the geometry of T
		N_GEO = N_GEO + 1
		seg_N_geo = seg_N_geo + 1
		!write(*,*) ".. Case Geometry : second randow draw move/birth/death "
		alpha_T = ran3(ra)   							! tir du type de perturbation alpha_T
		!write(*,*) alpha_T
		if (alpha_T.lt.1./3) then 									! CASE MOVE
			if (Nc/=0) then 
			  	N_move = N_move +1
			  	seg_N_move =  seg_N_move +1 
			  	index_pert = min(Nc, 1 + int(ran3(ra)*Nc))			   ! on tire aleatoire une des Nc perturbations
			  	!write(*,*) "... Case Move : selecting discontinuity... ", index_pert
			  	
			  	t_old = model_times(index_pert)							! on recupere le temps associe a cette perturbation
			  	t_new = t_old + gasdev(ra) * t_std						! on perturbe selon une gaussienne

			  	if ((t_new.ge.time_lim) .and. (t_new.le.(time_max-time_lim))) then 		! sinon do nothing
			  		New_model_times(index_pert) = t_new	
			  		index_ibins_old = ceiling(T_bins * t_old/time_max)			! on retrouve l'index associe a ce temps d'observation random dans le vecteur prop		
			  		index_ibins = ceiling(T_bins * t_new/time_max)				! on retrouve l'index associe a ce temps d'observation random dans le vecteur prop
			  		New_model_T(index_ibins_old) = 0
			  		New_model_T(index_ibins) = 1
			  		New_model_ibins(index_pert) = index_ibins
			  		index_data = 0
			    		do i = 1, N_data_mag 
			    			if (times(i) .le. t_new) then 					! on retrouve l'index associe a ce temps d'observation random dans les donnees
			    				index_data = i
			    			end if 
			    		end do
			    		New_model_idata(index_pert) = index_data
			    		New_Nc = Nc 
			  		!write(*,*) "...  Move accepted"
			  		do k = 1, New_Nc+1
			  			if ((k.lt.index_pert).or.(k.gt.index_pert+1)) then
			  				New_loglikelihood_layers(k) = Acc_loglikelihood_layers(k)

			  			else 
			  				New_loglikelihood_layers(k) = 0
			  			end if
			  		end do
			  	else 
			  		!write(*,*) "... Move rejected : bounds condition"
			  		flag_unaccepted = 1
			  	end if 
			else 
			  	!write(*,*) "... Move rejected : bounds condition"
			  	flag_unaccepted = 1
			end if 
		  	call sort_shell_real(New_model_times)
		  	call sort_shell(New_model_ibins)
		  	call sort_shell(New_model_idata)
		  	
		  	if (mod(it,10000000)==0) then ! talk
			  	write(*,*) "........................................."
			  	write(*,*) "... Current model ibins Move (new Timegrid)", model_ibins(1:Nc)
				write(*,*) "... Current model  idata Move (new Timegrid)", model_idata(1:Nc)
				write(*,*) "... Current model times Move (new Timegrid)", model_times(1:Nc)
				write(*,*) "... Current model beta ",model_beta(1:Nc+1)
				write(*,*) "... Current model loglike ",Acc_loglikelihood_layers(1:Nc+1) ! Attention changer 
				write(*,*)
				write(*,*) "........................................."
			  	write(*,*) "... Proposed model ibins Move (new Timegrid)", New_model_ibins(1:New_Nc)
				write(*,*) "... Proposed model  idata Move (new Timegrid)", New_model_idata(1:New_Nc)
				write(*,*) "... Proposed model times Move (new Timegrid)", New_model_times(1:New_Nc)
				write(*,*) "... Proposed model beta ",New_model_beta(1:New_Nc+1)
				write(*,*) "... Proposed model loglike ",New_loglikelihood_layers(1:New_Nc+1)
				write(*,*) "........................................."
			end if
			

				
		else if ((alpha_T.ge.1./3) .and. (alpha_T.lt.2./3))  then  					 !Case death 
		!else if (alpha_T.lt.1./2)  then 
		  	!write(*,*) "... Case Death : selecting discontinuity..."
		  	N_death = N_death +1
		  	seg_N_death = seg_N_death + 1
			if (Nc/=0) then 
				New_Nc = Nc-1
		  		index_pert = min(Nc, 1 + int(ran3(ra)*Nc))		
		  		t_new = model_times(index_pert)
		  		
		  		Acc_Beta = model_beta(index_pert)
		  		Acc_Mu = model_mu(index_pert)
		  		Acc_Sigma = model_sigma(index_pert)
		  		!write(*,*) "Layer removed :", t_new, index_pert, Acc_Beta, Acc_Mu, Acc_Sigma
		  		
		  		New_model_times(index_pert) = 9999999
		  		New_model_idata(index_pert) = 9999999
		  		New_model_ibins(index_pert) = 9999999
		  		
		  		index_ibins = ceiling(T_bins * t_new/time_max)
		  		New_model_T(index_ibins) = 0	! on supprime la discontinuite
		  		
		  		call sort_shell_real(New_model_times)
		  		call sort_shell(New_model_ibins)
		  		call sort_shell(New_model_idata)
		  		
		  		New_Beta = model_beta(index_pert+1)
		  		New_Mu = model_mu(index_pert+1)
		  		New_sigma = model_sigma(index_pert+1)
		  		
		  		
		  		do k = 1, New_Nc+1
		  			if (k.lt.index_pert) then
		  				New_loglikelihood_layers(k) = Acc_loglikelihood_layers(k)
		  				New_model_beta(k) = model_beta(k)
		  				New_model_mu(k) = model_mu(k)
		  				New_model_sigma(k) = model_sigma(k)
		  			else if (k.gt.index_pert) then
		  				New_loglikelihood_layers(k) = Acc_loglikelihood_layers(k+1)
		  				New_model_beta(k) = model_beta(k+1)
		  				New_model_mu(k) = model_mu(k+1)
		  				New_model_sigma(k) = model_sigma(k+1)
		  			else 
		  				New_loglikelihood_layers(k) = 0
		  				New_model_beta(k) = New_Beta
		  				New_model_mu(k) = New_Mu
		  				New_model_sigma(k) = New_Sigma
		  			end if
		  		end do
		  		New_model_beta(New_Nc+2) = 99 			! Essayer de supprimer
		  		New_model_mu(New_Nc+2) = 99
		  		New_model_sigma(New_Nc+2) = 99
		   	else 
		  		!write(*,*) "... Death rejected : already no layers "
		  		flag_unaccepted = 1
		  		New_Nc=0
		  	end if 
		  	
		  	if (mod(it,10000000)==0) then ! talk		
			  	write(*,*) "........................................."
			  	write(*,*) "... Current model ibins Death (new Timegrid)", model_ibins(1:Nc)
				write(*,*) "... Current model  idata Death (new Timegrid)", model_idata(1:Nc)
				write(*,*) "... Current model times Death (new Timegrid)", model_times(1:Nc)
				write(*,*) "... Current model beta ",model_beta(1:Nc+1)
				write(*,*) "... Current model mu ",model_mu(1:Nc+1)
				write(*,*) "... Current model sigma ",model_sigma(1:Nc+1)
				write(*,*) "... Current model loglike ",Acc_loglikelihood_layers(1:Nc+1) ! Attention changer 
				write(*,*)
				write(*,*) "........................................."
			  	write(*,*) "... Proposed model ibins Death (new Timegrid)", New_model_ibins(1:New_Nc)
				write(*,*) "... Proposed model  idata Death (new Timegrid)", New_model_idata(1:New_Nc)
				write(*,*) "... Proposed model times Death (new Timegrid)", New_model_times(1:New_Nc)
				write(*,*) "... Proposed model beta ",New_model_beta(1:New_Nc+1)
				write(*,*) "... Proposed model mu ",New_model_mu(1:New_Nc+1)
				write(*,*) "... Proposed model sigma ",New_model_sigma(1:New_Nc+1)
				write(*,*) "... Proposed model loglike ",New_loglikelihood_layers(1:New_Nc+1)
				write(*,*) "........................................."
				write(*,*) "........................................."
			end if
		  	
		else 													! Case birth
		  	!write(*,*) "... Case Birth : selecting a new discontinuity... "
		  	N_birth = N_birth + 1
		  	seg_N_birth = seg_N_birth+1
			index_pert = Nc + 1 
			!write(*,*) "Index pert", index_pert
			if (index_pert .le. max_layers) then 
				New_Nc = Nc +1
		  		t_new = time_lim + (ran3(ra)*(time_max-(2*time_lim)))
		  		New_model_times(index_pert) = t_new
		  		
				index_ibins = ceiling(T_bins * t_new/time_max)	! findig the index associated to this time in the time_bins = prop
		  		New_model_T(index_ibins) = 1
		  		New_model_ibins(index_pert) = index_ibins
		  		index_data = 0
		    		do i = 1, N_data_mag 
		    			if (times(i).le.t_new) then ! findig the index associated to this time in the dataset
		    				index_data = i
		    			end if 
		    		end do	
		    		New_model_idata(index_pert) = index_data
		  		! sorting the vectors of the temporal model
		  		call sort_shell_real(New_model_times)
		  		call sort_shell(New_model_ibins)
		  		call sort_shell(New_model_idata)
		  		
		  		
		  		! Searching for the index of the new discontinuity 
				do k=1,New_Nc
					if (New_model_times(k).eq.model_times(k)) then
						continue
					else
						index_pert = k
						!write(*,*) "The layer of the birth is layer", index_pert
						exit
					end if
				end do 
				
				!write(*,*) "Index pert", index_pert
				
				! Beta current vs proposed
				Acc_Beta = model_beta(index_pert)
				New_Beta = Acc_Beta + gasdev(ra) * beta_std_birth			! tirage aleatoire selon gaussienne (attention term_added)
				
				!New_Beta = beta_min + (ran3(ra)*(beta_max-beta_min))			! tirage aleatoire selon le prior (term_added=0)
				if ((New_Beta.lt.beta_min) .or. (New_Beta.gt.beta_max)) then 
		  			New_Beta = New_Beta
		  			flag_prior_beta = 1
		  		end if
		  		
		  		
		  		! Mu current vs proposed 
		  		Acc_Mu = model_mu(index_pert)
				New_Mu = Acc_Mu + gasdev(ra) * mu_std_birth				! tirage aleatoire selon gaussienne (attention term_added)
				if ((New_Mu.lt.mu_min) .or. (New_Mu.gt.mu_max)) then 
		  			New_Mu =  New_Mu
		  			flag_prior_mu = 1
		  		end if
		  		
		  		
		  		Acc_Sigma = model_sigma(index_pert)
				New_Sigma = Acc_Sigma + gasdev(ra) * sigma_std_birth			! tirage aleatoire selon gaussienne (attention term_added)
				!New_Sigma = sigma_min + (ran3(ra)*(sigma_max-sigma_min))		! tirage aleatoire selon le prior (term_added=0)
				if ((New_Sigma.lt.sigma_min) .or. (New_Sigma.gt.sigma_max)) then 
		  			New_Sigma =  New_Sigma
		  			flag_prior_sigma = 1
		  		end if
		  		
		  		if ((flag_prior_beta.eq.1).or.(flag_prior_mu.eq.1).or.(flag_prior_sigma.eq.1)) then
		  			flag_unaccepted = 1
		  			!write(*,*) New_Beta, New_Mu, New_Sigma
		  			New_Beta=Acc_Beta
		  			New_Mu = Acc_Mu
		  			New_Sigma = Acc_Sigma
		  			!write(*,*) "... Pas de birth : erreur sur priors GR"
		  		else
				
				
		  			!write(*,*) "... Birth accepted"
		  			
		  			do k = 1, New_Nc+1
			  			if (k.lt.index_pert) then
			  				New_loglikelihood_layers(k) = Acc_loglikelihood_layers(k)
			  				New_model_beta(k) = model_beta(k)
			  				New_model_mu(k) = model_mu(k)
			  				New_model_sigma(k) = model_sigma(k)
			  			else if (k.gt.index_pert+1) then
			  				New_loglikelihood_layers(k) = Acc_loglikelihood_layers(k-1)
			  				New_model_beta(k) = model_beta(k-1)
			  				New_model_mu(k) = model_mu(k-1)
			  				New_model_sigma(k) = model_sigma(k-1)
			  			else if (k.eq.index_pert) then
			  				New_loglikelihood_layers(k) = 0
			  				New_model_beta(k) = New_Beta
			  				New_model_mu(k) = New_Mu
			  				New_model_sigma(k) = New_Sigma
			  			else if (k.eq.index_pert+1) then
			  				New_loglikelihood_layers(k) = 0
			  				New_model_beta(k) = model_beta(k-1)
			  				New_model_mu(k) = model_mu(k-1)
			  				New_model_sigma(k) = model_sigma(k-1)
			  			end if
		  			end do
		  		end if
		  	else 
		  		!write(*,*) "... Pas de birth : conditions aux limites nb de couches"
		  		flag_unaccepted = 1
		  	end if 
		  	if (mod(it,10000000)==0) then ! talk
			  	write(*,*) "........................................."
			  	write(*,*) "... Current model ibins Birth (new Timegrid)", model_ibins(1:Nc)
				write(*,*) "... Current model  idata Birth (new Timegrid)", model_idata(1:Nc)
				write(*,*) "... Current model times Birth (new Timegrid)", model_times(1:Nc)
				write(*,*) "... Current model beta ",model_beta(1:Nc+1)
				write(*,*) "... Current model mu ",model_mu(1:Nc+1)
				write(*,*) "... Current model sigma ",model_sigma(1:Nc+1)
				write(*,*) "... Current model loglike ",Acc_loglikelihood_layers(1:Nc+1) 
				write(*,*)
				write(*,*) "........................................."
			  	write(*,*) "... Proposed model ibins Birth (new Timegrid)", New_model_ibins(1:New_Nc)
				write(*,*) "... Proposed model  idata Birth (new Timegrid)", New_model_idata(1:New_Nc)
				write(*,*) "... Proposed model times Birth (new Timegrid)", New_model_times(1:New_Nc)
				write(*,*) "... Proposed model beta ",New_model_beta(1:New_Nc+1)
				write(*,*) "... Proposed model mu ",New_model_mu(1:New_Nc+1)
				write(*,*) "... Proposed model sigma ",New_model_sigma(1:New_Nc+1)
				write(*,*) "... Proposed model loglike ",New_loglikelihood_layers(1:New_Nc+1)
				write(*,*) "........................................."
			end if
		end if



	else ! Case we perturb the GR parameters

		N_GR = N_GR + 1
		seg_N_GR = seg_N_GR +1	
	  	New_Nc = Nc		! the number of discontinuities does not change	
	  	
	  	alpha_GR = ran3(ra)
	  	index_pert = min(New_Nc+1, 1 + int(ran3(ra)*(New_Nc+1)))		! we randomly draw one of the Nc+1 segments
		Acc_Beta = model_beta(index_pert)			! we retrieve the current value of Beta in the concerned segment
		Acc_Mu = model_mu(index_pert)				! we retrieve the current value of Mu in the concerned segment
		Acc_Sigma = model_sigma(index_pert)			! we retrieve the current value of Sigma in the concerned segment
		!New_loglikelihood_layers = Acc_loglikelihood_layers	! temporary
	  	!write(*,*) ".. Case Perturbation Model Parameters on layer :", index_pert
		if (alpha_GR.lt.1./3) then 							! GR parameters > Beta move
			!write(*,*) "... Case beta : changing beta in the layer"
			N_beta = N_beta + 1
			seg_N_beta = seg_N_beta+1
		  	! update proposal Beta/Mu/sigma
		  	New_Beta = Acc_Beta + gasdev(ra) * beta_std
		  	New_Mu = Acc_Mu
		  	New_Sigma = Acc_Sigma
		  	
		  	if ((New_Beta.lt.beta_min) .or. (New_Beta.gt.beta_max)) then 
		  		New_Beta =  Acc_Beta
		  		flag_unaccepted = 1
		  	end if
		  	New_model_beta(index_pert) = New_Beta
		  	!write(*,*) "... Proposed model beta ",New_model_beta(1:New_Nc+1)
		
		else if ((alpha_GR.ge.1./3) .and. (alpha_GR.lt.2./3))then 			! GR parameters > Mu move
			!write(*,*) "... Case mu : changing mu in the layer"
			N_mu = N_mu + 1
			seg_N_mu = seg_N_mu +1
		  	! update proposal Beta/Mu/sigma
		  	New_Beta = Acc_Beta
		  	New_Mu = Acc_Mu + gasdev(ra) * mu_std
		  	New_Sigma = Acc_Sigma
		  	if ((New_Mu.lt.mu_min) .or. (New_Mu.gt.mu_max)) then 
		  		New_Mu = Acc_Mu
		  		flag_unaccepted = 1
		  	end if
		  	New_model_mu(index_pert) = New_Mu
		  	!write(*,*) "... Proposed model mu ",New_model_mu(1:New_Nc+1)
		 
		  
		  else 										! GR parameters > Sigma move
		  	!write(*,*) "... Case sigma : changing sigma in the layer"
			N_sigma = N_sigma + 1
			seg_N_sigma = seg_N_sigma+1
		  	
		  	New_Beta = Acc_Beta
		  	New_Mu = Acc_Mu 
		  	New_Sigma = Acc_Sigma + gasdev(ra) * sigma_std
		  	if ((New_Sigma.lt.sigma_min) .or. (New_Sigma.gt.sigma_max)) then 
		  		New_Sigma = Acc_Sigma
		  		flag_unaccepted = 1
		  	end if
		  	New_model_sigma(index_pert) = New_Sigma
		  	!write(*,*) "... Proposed model sigma ",New_model_sigma(1:New_Nc+1)
		end if
		if (mod(it,10000000)==0) then ! choose when to talk for debug
			write(*,*) "........................................."
		  	write(*,*) "... Current model ibins Move (new Timegrid)", model_ibins(1:Nc)
			write(*,*) "... Current model  idata Move (new Timegrid)", model_idata(1:Nc)
			write(*,*) "... Current model times Move (new Timegrid)", model_times(1:Nc)
			write(*,*) "... Current model beta ",model_beta(1:Nc+1)
			write(*,*) "... Current model mu ",model_mu(1:Nc+1)
			write(*,*) "... Current model sigma ",model_sigma(1:Nc+1)
			write(*,*) "... Current model loglike ",Acc_loglikelihood_layers(1:Nc+1) 
			write(*,*)
			write(*,*) "........................................."
		  	write(*,*) "... Proposed model ibins Move (new Timegrid)", New_model_ibins(1:New_Nc)
			write(*,*) "... Proposed model  idata Move (new Timegrid)", New_model_idata(1:New_Nc)
			write(*,*) "... Proposed model times Move (new Timegrid)", New_model_times(1:New_Nc)
			write(*,*) "... Proposed model beta ",New_model_beta(1:New_Nc+1)
			write(*,*) "... Proposed model mu ",New_model_mu(1:New_Nc+1)
			write(*,*) "... Proposed model sigma ",New_model_sigma(1:New_Nc+1)
			write(*,*) "... Proposed model loglike ",New_loglikelihood_layers(1:New_Nc+1)
			write(*,*) "........................................."
		end if


	end if  
	  
	! ----------------------------------------------------------------------------------------------------------------------------------------------------------
	  
	  
	if (flag_unaccepted == 0) then
	   New_model_T = 0
	   do k = 1, New_Nc
	      New_model_T(New_model_ibins(k)) = New_model_T(New_model_ibins(k)) + 1
	   end do
	end if
	  
	if (flag_unaccepted.eq.0) then 
	  
	  
	  !! 2.1 Compute loglikelihood of the new model

		do k = 1, New_Nc+1   ! Attention verifier si k=1 et New_Nc=0
			 if ((k==1).and.(k.ne.New_Nc+1)) then
			 	
		    		flag_ndata(1) = New_model_idata(1) 
		    		a_data = 1 
		      		b_data = a_data + flag_ndata(1) -1
		      		
		      		a_ibins = 1
		      		b_ibins = New_model_ibins(1)	
		    	
		    	else if ((k.ne.1).and.(k == New_Nc+1)) then
		    		flag_ndata(k) = N_data_mag - New_model_idata(k-1)  
		    		a_data = a_data + flag_ndata(k-1)
			    	b_data = a_data + flag_ndata(k) - 1
		    		
		    		a_ibins = New_model_ibins(k-1)
		      		b_ibins = T_bins
		      		
		      	else if ((k.eq.1).and.(k == New_Nc+1)) then
		    		flag_ndata(k) = N_data_mag  
		    		a_data = 1
			    	b_data = N_data_mag
		    		
		    		a_ibins = 1
		      		b_ibins = T_bins
		      		
		    	else
		    		flag_ndata(k) = New_model_idata(k) - New_model_idata(k-1) 
		    		a_data = a_data + flag_ndata(k-1)
			    	b_data = a_data + flag_ndata(k) - 1
			    	
			      	a_ibins = New_model_ibins(k-1)
			      	b_ibins = New_model_ibins(k)
		    	end if

		    	
		    	if (flag_ndata(k).eq.0) then 
			    		loglikelihood = 0
			else
			    		
					N_data_subset = size(mag_obs(a_data:b_data))
					
					call LogLikelihood_function_nobin_detect(mag_obs(a_data:b_data),N_data_subset, min_mags,&
									 New_model_beta(k),New_model_mu(k), New_model_sigma(k), loglikelihood)
					
				
			end if	
				
			    	New_loglikelihood_layers(k) = loglikelihood
			    	if (mod(it,10000000)==0) then ! talk
		    			write(*,*) "... Computing loglikelihood on new segment :", k, a_ibins, b_ibins, loglikelihood
		    		end if
		    	

		  end do
		  
		  
		  !if (sum(flag_ndata(1:New_Nc+1)).ne.N_data_mag) stop "HEIN" 
		  !write(*,*) ". Initial Log-likelihood :", Acc_loglikelihood_layers(1:Nc+1)
		  Acc_Sum_Loglikelihood =  sum(Acc_loglikelihood_layers(1:Nc+1))  ! Normal on compare le proposed index_pert+1 avec le current index_pert
		  New_Sum_Loglikelihood  = sum(New_loglikelihood_layers(1:New_Nc+1))
		  diff_log = New_Sum_Loglikelihood  - Acc_Sum_Loglikelihood
		  if (prior_only) diff_log = 0.0d0
		  
		  if (.not.(New_Sum_Loglikelihood > -huge(1.0d0))) flag_unaccepted = 1

		 
		  if (alpha_A.lt.1./2) then
	    		if (alpha_T.lt.1./3) then
	    			term_added = 0 
	    			
			else if ((alpha_T.ge.1./3) .and. (alpha_T.lt.2./3)) then ! Case Death

				terme_beta_1  = beta_delta /(beta_std_birth *sqrt(2.0d0*pi))
				terme_beta_2  = (New_Beta  - Acc_Beta )**2/(2.0d0*beta_std_birth**2)
				terme_mu_1    = mu_delta   /(mu_std_birth   *sqrt(2.0d0*pi))
				terme_mu_2    = (New_Mu    - Acc_Mu   )**2/(2.0d0*mu_std_birth**2)
				terme_sigma_1 = sigma_delta/(sigma_std_birth*sqrt(2.0d0*pi))
				terme_sigma_2 = (New_Sigma - Acc_Sigma)**2/(2.0d0*sigma_std_birth**2)
				term_added = log(terme_beta_1) + log(terme_mu_1) + log(terme_sigma_1) &
				       - terme_beta_2 - terme_mu_2 - terme_sigma_2
				
			else							! Case Birth
			
				terme_beta_1  = beta_std_birth *sqrt(2.0d0*pi)/beta_delta
				terme_beta_2  = (New_Beta  - Acc_Beta )**2/(2.0d0*beta_std_birth**2)
				terme_mu_1    = mu_std_birth   *sqrt(2.0d0*pi)/mu_delta
				terme_mu_2    = (New_Mu    - Acc_Mu   )**2/(2.0d0*mu_std_birth**2)
				terme_sigma_1 = sigma_std_birth*sqrt(2.0d0*pi)/sigma_delta
				terme_sigma_2 = (New_Sigma - Acc_Sigma)**2/(2.0d0*sigma_std_birth**2)
				
			    	term_added = log(terme_beta_1) + log(terme_mu_1) + log(terme_sigma_1) &
				       + terme_beta_2 + terme_mu_2 + terme_sigma_2                                                       
	
			end if
	    	  else 
	    		term_added = 0
	  
	  	  end if	  
		  condition = term_added + diff_log	
	end if
	  
	  
	  
	  
	  
	  !! 2.2 TEST CONDITION MCMC - Metropolis Hastings
	  !-------------------------------------------------------------------------
	  
	  u = log(ran3(ra))
	  if (mod(it,10000)==0) then ! talk
		  write(*,*) "-------------------------------------------------------------------------"
		  write(*,*) ". Metropolis-Hastings test : new / old max loglikelihood"
		  write(*,*) ". Metropolis-Hastings test : ", New_Sum_Loglikelihood , Acc_Sum_Loglikelihood
		  write(*,*) ". Term added: ", term_added
		  write(*,*) " Acceptance probability :", condition
		  write(*,*) " Random value for MH :", u
	  end if
		  
	  
	  

	  if ((u.gt.condition) .or. (flag_unaccepted ==1) ) then  ! rejected
	     	!write(*,*) ".. Model discarded (old>new) ! "  
    		flag_unaccepted = 0
		flag_std = ceiling(real(it )/real(it_std))
		flag_thin = ceiling(real(it )/real(it_thin))
		
		if (N_thin .lt. it_thin) then 
			N_thin = N_thin+1
			do k = 1, Nc+1
	
				if ((k==1).and.(k.ne.Nc+1)) then
					flag_ndata(1) = model_idata(1) 
    					a_data = 1 
      					b_data = a_data + flag_ndata(1) - 1
      		
      					a_ibins = 1
      					b_ibins = model_ibins(1)	
    	
    				else if ((k.eq.Nc+1).and.(k.ne.1)) then
    					flag_ndata(k) = N_data_mag - model_idata(k-1)  
    					a_data = a_data + flag_ndata(k-1)
	    				b_data = a_data + flag_ndata(k) - 1
    		
    					a_ibins = model_ibins(k-1)
      					b_ibins = T_bins
      		
      		
      				else if ((k.eq.1).and.(k.eq.Nc+1)) then
    					flag_ndata(k) = N_data_mag
    					a_data = 1
	    				b_data = N_data_mag
    		
    					a_ibins = 1
      					b_ibins = T_bins
      		
    				else
    					flag_ndata(k) = model_idata(k) - model_idata(k-1) 
    					a_data = a_data + flag_ndata(k-1)
	    				b_data = a_data + flag_ndata(k) - 1
	    	
	      				a_ibins = model_ibins(k-1)
	      				b_ibins = model_ibins(k)
    				end if
    				Acc_Beta = model_beta(k)
				Acc_Mu = model_mu(k)
				Acc_Sigma = model_sigma(k)
				if (mod(it,it_std)==0) then
					Store_Beta_proposed(k, flag_std+1) = Acc_Beta
					Store_Mu_proposed(k, flag_std+1) = Acc_Mu
					Store_Sigma_proposed(k, flag_std+1) = Acc_Sigma
					Store_times_proposed(k, flag_std+1) = model_times(k)
				end if 
			end do

		else
			N_thin = 1  
			Nc_total(flag_thin) = Nc
			
			
			if (it .gt. it_burnin) then
				do j = 1, T_bins	
					Sum_accepted_T(j) = Sum_accepted_T(j) + model_T(j)
		
				end do
			end if 
			do k = 1, Nc+1
	
				if ((k==1).and.(k.ne.Nc+1)) then
					flag_ndata(1) = model_idata(1) 
    					a_data = 1 
      					b_data = a_data + flag_ndata(1) - 1
      		
      					a_ibins = 1
      					b_ibins = model_ibins(1)	
    	
    				else if ((k.eq.Nc+1).and.(k.ne.1)) then
    					flag_ndata(k) = N_data_mag - model_idata(k-1)  
    					a_data = a_data + flag_ndata(k-1)
	    				b_data = a_data + flag_ndata(k) - 1
    		
    					a_ibins = model_ibins(k-1)
      					b_ibins = T_bins
      		
      		
      				else if ((k.eq.1).and.(k.eq.Nc+1)) then
    					flag_ndata(k) = N_data_mag
    					a_data = 1
	    				b_data = N_data_mag
    		
    					a_ibins = 1
      					b_ibins = T_bins
      		
    				else
    					flag_ndata(k) = model_idata(k) - model_idata(k-1) 
    					a_data = a_data + flag_ndata(k-1)
	    				b_data = a_data + flag_ndata(k) - 1
	    	
	      				a_ibins = model_ibins(k-1)
	      				b_ibins = model_ibins(k)
    				end if
				Acc_Beta = model_beta(k)
				Acc_Mu = model_mu(k)
				Acc_Sigma = model_sigma(k)
				
				bin_b     = min(beta_bins,  max(1, int((Acc_Beta  - beta_min )/beta_store )+1))
				bin_mu    = min(mu_bins,    max(1, int((Acc_Mu    - mu_min   )/mu_store   )+1))
				bin_sigma = min(sigma_bins, max(1, int((Acc_Sigma - sigma_min)/sigma_store)+1))
				
				if (mod(it,it_std)==0) then
					Store_Beta_proposed(k, flag_std+1) = Acc_Beta
					Store_Mu_proposed(k, flag_std+1) = Acc_Mu
					Store_Sigma_proposed(k, flag_std+1) = Acc_Sigma
					Store_times_proposed(k, flag_std+1) = model_times(k)
				end if 
				
				if (it .gt. it_burnin) then
					j_start = a_ibins
					if (k > 1) j_start = a_ibins + 1
					do j = j_start, b_ibins
						Beta_T_posterior( bin_b,j) = Beta_T_posterior( bin_b,j)+1 
						Mu_T_posterior( bin_mu,j) = Mu_T_posterior( bin_mu,j)+1 
						Sigma_T_posterior( bin_sigma,j) = Sigma_T_posterior( bin_sigma,j)+1 
					end do
				end if 
			end do
		end if
				
		    	

			
		
	  
	  else 
	  	Nc = New_Nc
		
		model_T = New_model_T		
    		model_times = New_model_times
    		model_idata = New_model_idata
    		model_ibins = New_model_ibins
    		model_beta = New_model_beta
    		model_mu = New_model_mu
    		model_sigma = New_model_sigma
		Acc_loglikelihood_layers = New_loglikelihood_layers
		flag_thin = ceiling(real(it )/real(it_thin))
		flag_std = ceiling(real(it )/real(it_std))
		

		if (N_thin .lt. it_thin) then 
			N_thin = N_thin+1
			do k = 1, Nc+1
				
		 		if ((k==1).and.(k.ne.Nc+1)) then
	    				flag_ndata(1) = New_model_idata(1) 
			    		a_data = 1 
			      		b_data = a_data + flag_ndata(1) - 1
			      		
			      		a_ibins = 1
			      		b_ibins = New_model_ibins(1)	
			    	
			    	else if ((k.eq.Nc+1).and.(k.ne.1)) then
			    		flag_ndata(k) = N_data_mag - New_model_idata(k-1)  
			    		a_data = a_data + flag_ndata(k-1)
				    	b_data = a_data + flag_ndata(k) - 1
			    		
			    		a_ibins = New_model_ibins(k-1)
			      		b_ibins = T_bins
			      		
			      		
			      	else if ((k.eq.1).and.(k.eq.Nc+1)) then
			    		flag_ndata(k) = N_data_mag
			    		a_data = 1
				    	b_data = N_data_mag
			    		
			    		a_ibins = 1
			      		b_ibins = T_bins
			      		
			    	else
			    		flag_ndata(k) = New_model_idata(k) - New_model_idata(k-1) 
			    		a_data = a_data + flag_ndata(k-1)
				    	b_data = a_data + flag_ndata(k) - 1
				    	
				      	a_ibins = New_model_ibins(k-1)
				      	b_ibins = New_model_ibins(k)
			    	end if

			    				
				
				Acc_Beta = New_model_beta(k)
				Acc_Mu = New_model_mu(k)
				Acc_Sigma = New_model_sigma(k)
				
				if (mod(it,it_std)==0) then
					Store_Beta_proposed(k, flag_std+1) = Acc_Beta
					Store_Mu_proposed(k, flag_std+1) = Acc_Mu
					Store_Sigma_proposed(k, flag_std+1) = Acc_Sigma
					Store_times_proposed(k, flag_std+1) = model_times(k)
				end if 
			end do

		else 
			N_thin = 1
			Nc_total(flag_thin) = New_Nc	
			if (it .gt. it_burnin) then
				do j = 1, T_bins	
    					Sum_accepted_T(j) = Sum_accepted_T(j) + New_model_T(j)
   				end do
   			end if
   			do k = 1, Nc+1
				
		 		if ((k==1).and.(k.ne.Nc+1)) then
	    				flag_ndata(1) = New_model_idata(1) 
			    		a_data = 1 
			      		b_data = a_data + flag_ndata(1) - 1
			      		
			      		a_ibins = 1
			      		b_ibins = New_model_ibins(1)	
			    	
			    	else if ((k.eq.Nc+1).and.(k.ne.1)) then
			    		flag_ndata(k) = N_data_mag - New_model_idata(k-1)  
			    		a_data = a_data + flag_ndata(k-1)
				    	b_data = a_data + flag_ndata(k) - 1
			    		
			    		a_ibins = New_model_ibins(k-1)
			      		b_ibins = T_bins
			      		
			      		
			      	else if ((k.eq.1).and.(k.eq.Nc+1)) then
			    		flag_ndata(k) = N_data_mag
			    		a_data = 1
				    	b_data = N_data_mag
			    		
			    		a_ibins = 1
			      		b_ibins = T_bins
			      		
			    	else
			    		flag_ndata(k) = New_model_idata(k) - New_model_idata(k-1) 
			    		a_data = a_data + flag_ndata(k-1)
				    	b_data = a_data + flag_ndata(k) - 1
				    	
				      	a_ibins = New_model_ibins(k-1)
				      	b_ibins = New_model_ibins(k)
			    	end if

			    				
				
				Acc_Beta = New_model_beta(k)
				Acc_Mu = New_model_mu(k)
				Acc_Sigma = New_model_sigma(k)
				bin_b     = min(beta_bins,  max(1, int((Acc_Beta  - beta_min )/beta_store )+1))
				bin_mu    = min(mu_bins,    max(1, int((Acc_Mu    - mu_min   )/mu_store   )+1))
				bin_sigma = min(sigma_bins, max(1, int((Acc_Sigma - sigma_min)/sigma_store)+1))
				if (mod(it,it_std)==0) then
					Store_Beta_proposed(k, flag_std+1) = Acc_Beta
					Store_Mu_proposed(k, flag_std+1) = Acc_Mu
					Store_Sigma_proposed(k, flag_std+1) = Acc_Sigma
					Store_times_proposed(k, flag_std+1) = model_times(k)
				end if 
				if (it .gt. it_burnin) then
					j_start = a_ibins
					if (k > 1) j_start = a_ibins + 1
					do j = j_start, b_ibins
						Beta_T_posterior( bin_b,j) = Beta_T_posterior(bin_b, j)+1
						Mu_T_posterior( bin_mu,j) = Mu_T_posterior( bin_mu,j)+1 
						Sigma_T_posterior( bin_sigma,j) = Sigma_T_posterior( bin_sigma,j)+1
					end do
   				end if
			end do
			
    		end if
 
    		if (mod(it,10000000)==0) then ! talk
			write(*,*) "------------------------------------------------------------"
	    		write(*,*) ".. Accepted model times (old<new) : ", model_times(1:Nc)
			write(*,*) ".. Accepted model beta ",model_beta(1:Nc+1)
			write(*,*) ".. Accepted model mu ",model_mu(1:Nc+1)
			write(*,*) ".. Accepted model sigma ",model_sigma(1:Nc+1)
		  	write(*,*) "------------------------------------------------------------"
		end if
	  	
	  	
	  	! Update acceptance rates
		    if (alpha_A .lt. 1./2) then
			N_GEO_accept     = N_GEO_accept + 1
			seg_N_GEO_accept = seg_N_GEO_accept + 1
			if (alpha_T .lt. 1./3) then
			    N_move_accept     = N_move_accept + 1
			    seg_N_move_accept = seg_N_move_accept + 1
			else if ((alpha_T .ge. 1./3) .and. (alpha_T .lt. 2./3)) then
			    N_death_accept     = N_death_accept + 1
			    seg_N_death_accept = seg_N_death_accept + 1
			else
			    N_birth_accept     = N_birth_accept + 1
			    seg_N_birth_accept = seg_N_birth_accept + 1
			end if
		    else
			N_GR_accept     = N_GR_accept + 1
			seg_N_GR_accept = seg_N_GR_accept + 1
			if (alpha_GR .lt. 1./3) then
			    N_beta_accept     = N_beta_accept + 1
			    seg_N_beta_accept = seg_N_beta_accept + 1
			else if ((alpha_GR .ge. 1./3) .and. (alpha_GR .lt. 2./3)) then
			    N_mu_accept     = N_mu_accept + 1
			    seg_N_mu_accept = seg_N_mu_accept + 1
			else
			    N_sigma_accept     = N_sigma_accept + 1
			    seg_N_sigma_accept = seg_N_sigma_accept + 1
			end if
		    end if
	  end if
	  
	  
	  
	 !! 2.3 ON TIRE UN NOUVEAU MODELE DE DISCONTINUITES
    	!-------------------------------------------------------------------------
	CALL cpu_time(t3)  !Toc. end counting time
	progress = (t3-t1)
	t_it(it) = (t3-t4)
	mean_t_it = (progress/it) * it_max
	
	if (mod(it,10000)==0) then ! talk
		write(*,*) "-------------------------------------------------------------------------"
		write(*,*) "Time (secondes) : ", progress ,'/', mean_t_it
		write(*,*) "Time (minuts) : ", (progress/60) ,'/', (mean_t_it/60)
		write(*,*) "< End of iteration >"
	  	write(*,*)
	    	write(*,*)
		if (it.EQ.it_max) then
			write(*,*) "> Finalisation :", it
		else
			write(*,*)
			write(*,*)
	    		write(*,*) "> New iteration : ", it
		end if
	end if
    	! reinitialisation des variables temporaires
    	
    	 
  end do  
!----------------------------------------------------------------------------------------------------
! 2- FIN TRANSDIMENSIONAL ITERATIONS
!------------------------------------------------------------------------------------------------------ 
  CALL cpu_time(t2)  !Toc. end counting time
  comp_time = nint((t2-t1)*100)
  write(*,*) "< End of reversible-jump >"
  write(*,*) " -------------------------"
  write(*,*) "> End of inversion (resume)" 
  write(*,*) ". Total time of inversion in minutes : " , (comp_time/100.0/60.0)
  write(*,*)
  write(*,*) size_thin, flag_thin
  ! Petit bilan des taux d'acceptance 
  AR_GEO = real(N_GEO_accept)/real(N_GEO)*100
  AR_move = real(N_move_accept)/real(N_move)*100
  AR_death = real(N_death_accept)/real(N_death)*100
  AR_birth = real(N_birth_accept)/real(N_birth)*100
  AR_GR = real(N_GR_accept)/real(N_GR)*100
  AR_beta = real(N_beta_accept)/real(N_beta)*100
  AR_mu = real(N_mu_accept)/real(N_mu)*100
  AR_sigma = real(N_sigma_accept)/real(N_sigma)*100
  write(*,*) " Final acceptance rates : "
  write(*,*) " Name : number_accepted  number_proposed    rate(pct)"
  write(*,*) " # GEOMETRY :", N_GEO_accept, N_GEO, AR_GEO
  write(*,*) "Births : ", N_birth_accept, N_birth, AR_birth
  write(*,*) "Deaths : ", N_death_accept, N_death, AR_death
  write(*,*) "Moves  : ", N_move_accept, N_move, AR_move
  write(*,*) " # GR parameters :", N_GR_accept, N_GR, AR_GR
  write(*,*) "Beta  : ", N_beta_accept, N_beta, AR_beta
  write(*,*) "Mu  : ", N_mu_accept, N_mu, AR_mu
  write(*,*) "Sigma  : ", N_sigma_accept, N_sigma, AR_sigma
  
!----------------------------------------------------------------------------------------------------
! #3- END OF PARALLELISATION OPENMPI
!------------------------------------------------------------------------------------------------------
  
  
!----------------------------------------------------------------------------------------------------
! #4- WRITING RESULTS - 
!------------------------------------------------------------------------------------------------------
  write(*,*) "..Writting outputs files... "
  close(75)
  open(unit=62, file="../outputs/Acceptances.txt",  status="replace", action="write")
	  write(62,*) " Acceptance rates at iteration : ", it_max
	  write(62,*) " Name : number_accepted  number_proposed    rate(pct)"
	  write(62,*) " # GEOMETRY :", N_GEO_accept, N_GEO, AR_GEO
	  write(62,*) "Births : ", N_birth_accept, N_birth, AR_birth
	  write(62,*) "Deaths : ", N_death_accept, N_death, AR_death
	  write(62,*) "Moves  : ", N_move_accept, N_move, AR_move
	  write(62,*) " # GR parameters :", N_GR_accept, N_GR, AR_GR
	  write(62,*) "Beta  : ", N_beta_accept, N_beta, AR_beta
	  write(62,*) "Mu  : ", N_mu_accept, N_mu, AR_mu
	  write(62,*) "Sigma  : ", N_sigma_accept, N_sigma, AR_sigma
  close(62)
  
  
  open(unit=63, file="../outputs/NC_tot.txt",  status="replace", action="write")
  do line = 1, size_thin
  	write(63, *) Nc_total(line)
  end do
  close(63)


  open(unit=69, file="../outputs/Beta_tot.txt",  status="replace", action="write")
  open(unit=70, file="../outputs/Mu_tot.txt",    status="replace", action="write")
  open(unit=71, file="../outputs/Sigma_tot.txt", status="replace", action="write")
  do line = 1, max_layers+2
   	do column = 1, size_std+1
	      write(69,'(F15.5,2X)',advance='no') Store_Beta_proposed(line, column)
	      write(70,'(F15.5,2X)',advance='no') Store_Mu_proposed(line, column)
	      write(71,'(F15.5,2X)',advance='no') Store_Sigma_proposed(line, column)
   	end do
   	write(69,*); write(70,*); write(71,*)
  end do
  close(69); close(70); close(71)

  open(unit=72, file="../outputs/Times_tot.txt", status="replace", action="write")
  do line = 1, max_layers+1
   	do column = 1, size_std+1
      		write(72,'(F15.5,2X)',advance='no') Store_times_proposed(line, column)
   	end do
   	write(72,*)
  end do
  close(72)  

  
  
  open(unit=64, file="../outputs/prop.txt",  status="replace", action="write")
  do line = 1, T_bins
  
  	write(64, *) Sum_accepted_T(line)
  end do
  close(64)
  

   open(unit=60, file="../outputs/hypercube.txt",  status="replace", action="write")
	write(60,*) mu_min, mu_max, mu_bins, mu_init,  mu_stdi
	write(60,*) beta_min, beta_max, beta_bins, beta_init, beta_stdi
	write(60,*) sigma_min, sigma_max, sigma_bins, sigma_init, sigma_stdi
	write(60,*) it_max, T_bins, comp_time/100.0/60.0
	write(60,*) t_offset, time_max
	write(60,*) beta_store, mu_store, sigma_store
 close(60)


  
  open(unit=55, file="../outputs/Likelihood_Mu.txt",  status="replace", action="write")
  do line = 1, mu_bins
  	do column = 1, T_bins
  		write(55,'(F50.5,2X)', advance='no') Mu_T_posterior(line, column)
  	end do
  	write(55,*)
  end do
  close(55)
  
  
  open(unit=56, file="../outputs/Likelihood_b.txt",  status="replace", action="write")
  do line = 1, beta_bins
  	do column = 1, T_bins
  		write(56, '(F50.5,2X)', advance='no') Beta_T_posterior(line, column)
  	end do
  	write(56,*)
  end do
  close(56)
  
  open(unit=59, file="../outputs/Likelihood_Sigma.txt",  status="replace", action="write")
  do line = 1, sigma_bins
  	do column = 1, T_bins
  		write(59, '(F50.5,2X)', advance='no') Sigma_T_posterior(line, column)
  	end do
  	write(59,*)
  end do
  close(59)

  write(*,*) ".. Completed"
  write(*,*) "< End of b-Bayesian program >"

!------------------------------------------------------------------------------------------------------
end program b_Bayesian_binned_McMC_detect










!-----------------------------------------
! UTILS
!-----------------------------------------


 
FUNCTION ran3(idum)
  INTEGER idum
  REAL ran3                 ! garder le type d'origine
  CALL random_number(ran3)
  ran3 = 1.0 - ran3
END FUNCTION ran3


FUNCTION GASDEV(idum)
	!     ..Arguments..
	integer          idum
	real GASDEV

	!     ..Local..
	real v1,v2,r,fac
	real ran3

	
	10 v1=2*ran3(idum)-1
	v2=2*ran3(idum)-1
	r=v1**2+v2**2
	if(r.ge.1.or.r.eq.0) GOTO 10
	fac=sqrt(-2*log(r)/r)
	GASDEV=v2*fac
RETURN
END

