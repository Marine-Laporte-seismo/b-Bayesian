program b_Bayesian

  implicit none

  !            ----------------------------------
  !              # BEGINING OF THE USER-DEFINED PARAMETERS
  !            ----------------------------------

  real , EXTERNAL    ::   ran3, gasdev

!-----------------------------------------
! # Data path
!----------------------------------------- 
  character(LEN=*), parameter :: Magnitude_file = '../data/synthetic_tests/ETAS_STAI_10000evts.in'
!-----------------------------------------
! # Parameters of the Markov chain
!-----------------------------------------
  integer, parameter :: itmax = 30000 		! Nombre d'iterations totales
  integer, parameter :: burnin = 5000		! Nombre d'iterations burned
  integer, parameter :: thinning = 1		! Thinning pour pour independance des modeles
  integer, parameter :: max_layers = 300 	! Nombre de couches max autorisees
  
  integer, parameter :: Nc_start_min = 10 	! Nombre de couches initiales
  integer, parameter :: Nc_start_max = 20 	! Nombre de couches initiales
  integer, parameter :: std = 10		! Pour la perturbation INITIALE de la discontinuite 
  						! = standard deviation autour de l'ancienne discontinuite (en pourcentage du tps max des data)
  integer, parameter :: it_std = 5000		! Nombre d'iterations entre deux changements de std
  integer, parameter :: T_bins = 200		! Taille du vecteur temps discretise
  integer, parameter :: ndata_max = 800000	! Nombre maximum de donnees lues dans le fichier donnees

!-----------------------------------------
! # Definition de l'hypercube a explorer
!-----------------------------------------
  integer, parameter :: n_tirs = 950000
  integer, parameter :: Sigma_bins = 100	! Taille de l'exploration du prior du parametre Sigma de la loi de detection (Gumbel)
  real*8, parameter :: Sigma_min = 0.01		! Valeur minimale de l'exploration du parametre Sigma de la loi de detection (Gumbel)
  real*8, parameter :: Sigma_max = 0.6		! Valeur maximale de l'exploration du parametre Sigma de la loi de detection (Gumbel)
  
  integer, parameter :: Mu_bins = 100		! Taille de l'exploration du prior du parametre Mu de la loi de detection (Gumbel)
  real*8, parameter :: Mu_min =-1		! Valeur minimale de l'exploration du parametre Mu de la loi de detection (Gumbel)
  real*8, parameter :: Mu_max = 5		! Valeur maximale de l'exploration du parametre Mu de la loi de detection (Gumbel)
  
  integer, parameter :: B_bins = 100		! Taille de l'exploration du prior du parametre Beta de la loi de Gutenberg-Richter (b_value*ln(10))
  real*8, parameter ::  B_min = 0.4*2.3	! Valeur minimale de l'exploration du parametre Beta de la loi de Gutenberg-Richter (b_value*ln(10))
  real*8, parameter ::  B_max = 1.6*2.3	! Valeur maximale de l'exploration du parametre Beta de la loi de Gutenberg-Richter (b_value*ln(10))
  
  
  
  
  
  
!-----------------------------------------
! # Interface subroutines : declaration des fonctions
!-----------------------------------------		
interface


  subroutine bMuSigma_rand(mag_obs, Mmin, n_tirs, bvec, Muvec, Sigmavec, &
				mat_out, maxval2, MLike_b, MLike_mu, MLike_sigma)
	  real*8, intent(in) :: mag_obs(:)
	  real*8, intent(in) :: Mmin
	  integer, intent(in):: n_tirs
	  real*8, intent(in) :: bvec(:), Muvec(:), Sigmavec(:)
	  real*8, intent(out) :: mat_out(n_tirs)
	  real*8, intent(out) :: maxval2
	  real*8, intent(out) :: MLike_b(size(bvec)), MLike_mu(size(Muvec)), MLike_sigma(size(Sigmavec))
  end subroutine bMuSigma_rand

  subroutine bMuSigma_full(mag_obs, Mmin, Svec, Muvec, bvec, mat1, maxval2)
  	real*8, intent(in) :: mag_obs(:)
	real*8, intent(in) :: Mmin
  	real*8, intent(in) :: Muvec(:)
  	real*8, intent(in) :: bvec(:)  
  	real*8, intent(in) :: Svec(:)
  	real*8, intent(out) :: mat1(size(Muvec), size(bvec),size(Svec))
  	real*8, intent(out) :: maxval2
  end subroutine bMuSigma_full
  
  
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
   
   
   subroutine likelihood_death_rand(mag_obs, Mmin, n_tirs, bvec, Muvec, Sigma_vec, cpert, tmp_Nc,T_bins, &
				      tmp_model_idata, tmp_model_iprop, & 
				      New_Like_Mu, New_Like_b, New_Like_Sigma, &
				      acc_Like_Mu, acc_Like_b, acc_Like_Sigma, &
				      Likelihood, MaxLikelihood, New_Like, New_MaxLike)
	real*8, intent(in) :: mag_obs(:)
	real*8, intent(in) :: Mmin
	integer, intent(in) :: n_tirs
	real*8, intent(in) :: bvec(:), Muvec(:), Sigma_vec(:)
	integer,intent(in) :: cpert, tmp_Nc, T_bins
	integer ,intent(in):: tmp_model_idata(:), tmp_model_iprop(:)
	
	real*8 ,intent(in):: Likelihood(:), MaxLikelihood(:)
	real*8 ,intent(in):: acc_Like_Mu(:,:), acc_Like_b(:,:), acc_Like_Sigma(:,:)
	real*8, intent(out):: New_Like(:), New_MaxLike(:)
	real*8, intent(out) :: New_Like_Mu(:,:), New_Like_b(:,:),New_Like_Sigma(:,:)
   end subroutine likelihood_death_rand
   subroutine likelihood_birth_rand(mag_obs,Mmin, n_tirs, bvec, Muvec, Sigma_vec, cpert, tmp_Nc,T_bins, &
				      tmp_model_idata, tmp_model_iprop, & 
				      New_Like_Mu, New_Like_b, New_Like_Sigma, &
				      acc_Like_Mu, acc_Like_b, acc_Like_Sigma, &
				      Likelihood, MaxLikelihood, New_Like, New_MaxLike)
	real*8, intent(in) :: mag_obs(:)
	real*8, intent(in) :: Mmin
	integer, intent(in) :: n_tirs
	real*8, intent(in) :: bvec(:), Muvec(:), Sigma_vec(:)
	integer,intent(in) :: cpert, tmp_Nc, T_bins
	integer ,intent(in):: tmp_model_idata(:), tmp_model_iprop(:)
	
	real*8 ,intent(in):: Likelihood(:), MaxLikelihood(:)
	real*8 ,intent(in):: acc_Like_Mu(:,:), acc_Like_b(:,:), acc_Like_Sigma(:,:)
	real*8, intent(out):: New_Like(:), New_MaxLike(:)
	real*8, intent(out) :: New_Like_Mu(:,:), New_Like_b(:,:),New_Like_Sigma(:,:)
   end subroutine likelihood_birth_rand
end interface  
   
!-----------------------------------------
! # DECLARATIONS OF VARIABLES
!-----------------------------------------
integer :: ra							! seed of our run (change if we want a different run)

real*8 :: mag_obs0(ndata_max), times0(ndata_max)   
real*8, allocatable :: mag_obs(:), times(:) 			! DONNEES
integer :: ndata_mag 						! nombre de donnees lues dans le fichier
real*8 :: tlim, tmax, tmin					! tps min et tps max lu dans le fichier
real*8 :: Mmin

real*8 ::  B_int, Mu_int, Sigma_int	
real*8 :: Muvec(Mu_bins), bvec(B_bins), Sigma_vec(Sigma_bins)		! Espace exploration de la Mu et de la b
real*8 :: mat1(size(Muvec), size(bvec),size(Sigma_vec))			! sorties du code bMuSigma_full
real*8 :: mat_out(n_tirs)						! sorties du code bMuSigma_rand
real*8 :: mat_MB(size(Muvec), size(bvec)) , mat_MS(size(Muvec), size(Sigma_vec)), mat_BS(size(bvec), size(Sigma_vec)) ! Likelihood 2D : full explo
real*8 :: maxval2							! sorties du code bMuSigma_full


integer :: Nc, tmp_Nc				! nombre de couches 
integer :: Nc_total(itmax)			! store nb de couches a chaque iteration


real*8 :: Likelihood(max_layers+1),MaxLikelihood(max_layers+1) 	! likelihood stackee
real*8 :: acc_like(max_layers+1),acc_maxlike(max_layers+1) 	! likelihood acceptee
real*8 :: New_Like(max_layers+1) , New_MaxLike(max_layers+1)	! nouvelle likelihood
real*8 :: tmp_Like(max_layers+1) , tmp_MaxLike(max_layers+1)	! likelihood temporaire : acceptee ou non
real*8 :: Likelihood_all, MaxLikelihood_all, log_all		! likelihood du jeu de donnees total
real*8 :: logTemp, logP(itmax), logPb(itmax), dMdB 		! LogLike pour calculer condition

real*8 :: model_times(max_layers+1)		! vecteur temps des perturbations
integer :: model_iprop(max_layers+1)		! vecteur indices des perturbations dans le vecteur prop
integer :: model_idata(max_layers+1)		! vecteur indices des perturbations dans le vecteur data 

real*8 :: tmp_model_times(max_layers+1)		! vecteur temps des perturbations temporaire
integer :: tmp_model_iprop(max_layers+1)	! vecteur indices des perturbations dans le vecteur prop temporaire
integer :: tmp_model_idata(max_layers+1)	! vecteur indices des perturbations dans le vecteur data temporaire 
integer :: tmp_model_ddata(max_layers+1)
integer :: tmp_model_dprop(max_layers+1)

integer :: a_segt, b_segt, a_prop, b_prop	! Indices pour decoupe de catalogues entre deux discontinuites temporelles

integer ::  index_data, index_prop, index_prop_mu, index_prop_s	! Indices des vecteurs donnees/prop


integer :: cpert					! Indice de la discontinuite a perturber
real*8 :: tp, s						! Temps de la discontinuite a perturber
real*8 :: tc, tc_init, tstd				! Temps de discontinuite uen fois perturbe selon tstd, temps minimal des perturbations, standard deviation temporelle

integer :: init_full
integer :: Npert, Nbirth, Ndeath , store_death		! Nombre de perturbations / birth / death 
integer :: Nperta, Nbirtha, Ndeatha			! Nombre de perturbations / birth / death acceptees
real*8 ::  Acc_pert, Acc_birth, Acc_death		! Enregistre le taux d'acceptation

integer ::  c, flag

integer :: prop(T_bins), Sumprop(T_bins)		! Matrice des discontinuites acc/ Stack des discontinuites acc
integer :: tmp_prop(T_bins)				! Matrice temporaire des discontinuites proposees

real*8 :: MLike_mu( Mu_bins), MLike_b(B_bins), MLike_sigma(Sigma_bins)			! Marginaux estimes par exploration random : sortie de BMuSigma_rand
real*8 :: Like_Mu(T_bins, Mu_bins), tmp_Like_Mu(T_bins, Mu_bins),Like_Mu_all(T_bins,Mu_bins) ! Marginaux Mu : stack/temporaire/all
real*8 :: New_Like_Mu(T_bins, Mu_bins), acc_Like_Mu(T_bins, Mu_bins)			! Marginaux Mu : nouveau/accepte
real*8 :: Like_b(T_bins, B_bins),  tmp_Like_b(T_bins, B_bins),Like_b_all(T_bins,B_bins) ! Marginaux B : stack/temporaire/all
real*8 :: New_Like_b(T_bins, B_bins),  acc_Like_b(T_bins, B_bins)			! Marginaux B : nouveau/accepte
real*8 :: Like_Sigma(T_bins, Sigma_bins),  tmp_Like_Sigma(T_bins, Sigma_bins),Like_Sigma_all(T_bins,Sigma_bins) ! Marginaux Sigma : stack/temporaire/all
real*8 :: New_Like_Sigma(T_bins, Sigma_bins),  acc_Like_Sigma(T_bins, Sigma_bins)	! Marginaux Sigma : stack/temporaire/all

! iteratives
integer :: j, i, ij, a1, a2, a3, k, kj	, kk , kp, iw, idata, tci, i_mag, k_mag, k_prop	
integer :: kb, kMu
integer :: io, ioo


integer ::  flag_unaccepted, flag_ndata(max_layers)			! some models are ill constructed - reject directly
real*8 :: alpha 							! condition birth/death/pert reversible jump
real*8 :: u, condition							! condition acceptation Metropolis Hastings
real*8:: New_std						! Methode T.Santos : adaptation du std en fonction de Acc_pert
integer :: N_thin						! Pour independance des parametres : garder les modeles acceptes tous les dix iterations


real*8 :: t1, t2, comp_time					! Pour estimer le temps de calcul
integer :: line_data, line61,  line52, line53, line60,line55,column55,line56,column56, line59, column59, line64, line63	! Ecriture des fichiers de sortie

!-----------------------------------------
! Variables initialisation
!-----------------------------------------

init_full = 0 ! For test if full = random during 1st iteration (init)
model_times =  9999999
model_iprop =  9999999
model_idata =  9999999


Likelihood=0
Sumprop=0
tmp_prop = 0 
tmp_model_times =  9999999
tmp_model_iprop =  9999999
tmp_model_idata =  9999999

Nbirth=0
Npert=0
Ndeath=0


idata=1
index_data =0
flag_ndata = 0


logP = 0.00
logPb = 0.00
logTemp = 0.00


Like_b=0
Like_Mu=0
Like_Sigma=0

N_thin = 0 

!-------------------------------------------------
! Etape 0.A : Affectation des vecteurs Muvec , Bvec, Sigma_vec
! Equivalent as np.linspace 
!-------------------------------------------------
do a1 = 1, Mu_bins
  Muvec(a1) = (a1-1) * (Mu_max-Mu_min)/(Mu_bins-1) + Mu_min
end do
do a2 = 1, B_bins
  bvec(a2) =  (a2-1) * (B_max-B_min)/(B_bins-1) + B_min
end do

do a3 = 1, Sigma_bins
  Sigma_vec(a3) =  (a3-1) * (Sigma_max-Sigma_min)/(Sigma_bins-1) + Sigma_min
end do

dMdB = size(Muvec)*size(bvec)*size(Sigma_vec)	! Poids Full explo
B_int = (B_max-B_min)/B_bins
Mu_int = (Mu_max-Mu_min)/Mu_bins
Sigma_int = (Sigma_max-Sigma_min)/Sigma_bins

!-----------------------------------------
! Etape 0.B : Read the magnitude data
! The magnitude catalog must only have two colums
! C1 : Time between 0 and TMAX in minuts or hour or day 
! C2 : Event magnitude
!-----------------------------------------

open(50,file= Magnitude_file ,status='old')
    ndata_mag = 0
    do line_data =1, ndata_max
	      read(50,*,IOSTAT=io)times0(line_data),mag_obs0(line_data)	
	      if (io > 0) then				
			    stop "Check input.  Something was wrong"
	      elseif (io < 0) then
			    exit
	      else
	      	ndata_mag = ndata_mag+1
	      endif
	      if (line_data == ndata_max) stop "number of magnitudes data >= ndata_magmax"
    enddo
 close(50)! close the file

!-----------------------------------------
! allocate de tous les tableaux qui dependent de ndata 
allocate(mag_obs(ndata_mag))
allocate(times(ndata_mag))
mag_obs = mag_obs0(1:ndata_mag)
times =  times0(1:ndata_mag)
Mmin = minval(mag_obs)
! initialisation of time parameters
tmax = maxval(times)
tmin = minval(times)
if (tmax .EQ. tmin ) then 
	if (tmin .NE. 0) then
		tmax = 2 * tmin 
	else 
		tmax = 1
	end if
else if (tmax == 0) then
	tmax = 1
end if
tlim =  tmax/100 * 1			! on accepte pas de discontinuite sur les bords std% du catalogue 
tstd = tmax/100  * std			! on perturbe la discontinuite de (std%) du catalogue


!------------------------------------------------------------
! Etape 0.C - WRITE PARAMETRISATION AND DATA IN ORDER TO KEEP IT IN MIND
!------------------------------------------------------------


 open(unit=61, file="../out/input.txt",  status="replace", action="write")
 do line61 = 1, ndata_mag
	write(61, *) times(line61),mag_obs(line61)
 end do
 close(61)
  
  
  
!------------------------------------------------------------------------
! 1- INITIALIZATION OF THE TRANSDIMENSIONAL INVERSION
! THE PURPOSE OF PARALLELISATION IS TO BEGIN AT DISTINCT INITIALIZATION 
!-----------------------------------------------------------------------
!-----------------------------------------
! #Start MPI
!-----------------------------------------


!call MPI_INIT(ierror)
!call MPI_COMM_SIZE(MPI_COMM_WORLD, nbproc, ierror)
!call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
!call MPI_Comm_group(MPI_COMM_WORLD,group_world,ierror)



  CALL cpu_time(t1)  !Tic. start counting time 
	
  !!1.1 Initialization !!
  !-----------------------------------------
  write(*,*) "> Begin __init__... "
  write(*,*) "------------------- "
  Nc = Nc_start_min + (ran3(ra)*(Nc_start_max-Nc_start_min))
  write(*,*) ". Initial number of discontinuities ", Nc


  do i=1, Nc
	
    	tc_init = tlim + (ran3(ra)*(tmax-(2*tlim)))		! temps d'observations random entre tpert_lim et tmax-tpert_lim (REEL)
	!tc_init = (tmin+tlim) + (ran3(ra)*((tmax-tlim)-(tmin+tlim)))
    	index_prop = ceiling(T_bins * tc_init/tmax)		! on retrouve l'index associe a ce temps d'observation random dans le vecteur prop  	
    	
    	! on retrouve l'index associe a ce temps d'observation random dans les donnees
    	idata = 1
    	index_data =1
    	do while (((idata) .lt. ndata_mag) .and. (times(idata) .le. tc_init) ) ! on retrouve l'index associe a ce temps d'observation random dans les donnees
    		index_data = idata
    		idata = idata + 1
    	end do							
    	model_times(i) = tc_init				! Temps des discontinuites tirees (ne correspond pas forcement a un temps du catalogue) = reels					
    	model_iprop(i) = index_prop				! on garde en memoire les deux index 
    	model_idata(i) = index_data
    	prop(index_prop) = 1 					! assigne 1 à l'endroit où l'on place une discontinuite [0,0,0,1,0,...] dans le vecteur prop(T_bins)
    	Nc_total(1) = Nc						! on garde en memoire le nombre de couches acceptees Nc a chaque iteration 
   end do
   call sort_shell_real(model_times)			! et on replace tout dans le bon ordre pour creer le modele de discontinuites !! 
   call sort_shell(model_iprop)
   call sort_shell(model_idata)
	
   write(*,*) ". Initial times of discontinuities               ", model_times(1:Nc)
   write(*,*) ". Initial position of discontinuities in Timegrid", model_iprop(1:Nc)
   write(*,*) ". Initial position of discontinuities in data    ", model_idata(1:Nc)
   write(51, *) model_times(1:Nc)


 
  !! 1.2 On calcule la likelihood pour chaque segment !! 
  !------------------------------------------------------
  ! 1.2.1 On cherche les indices des segments
  write(*,*) ". Now, computing first likelihood..."
  write(*,*)".. Likelihood :       id        begin       end       n_data	    Like_b(sample)     MaxLogLike       LogLike(1)"
  do ij = 1, Nc+1
      	if (ij==1) then
    		flag_ndata(1) = model_idata(1) 
    		a_segt = 1 
      		b_segt = a_segt + flag_ndata(1) - 1
      		
      		a_prop = 1
      		b_prop = model_iprop(ij)	
    	
    	else if (ij == Nc+1) then
    		flag_ndata(ij) = ndata_mag - model_idata(ij-1)  
    		a_segt = a_segt + flag_ndata(ij-1)
	    	b_segt = a_segt + flag_ndata(ij) - 1
    		
    		a_prop = model_iprop(ij-1)
      		b_prop = T_bins
    	else
    		flag_ndata(ij) = model_idata(ij) - model_idata(ij-1) 
    	
    		a_segt = a_segt + flag_ndata(ij-1)
	    	b_segt = a_segt + flag_ndata(ij) - 1
	    	
	      	a_prop = model_iprop(ij-1)
	      	b_prop = model_iprop(ij)
    	end if
    	
    	!! ON COMPARE FULL VS RANDOM POUR LA PREMIERE EXPLORATION 
    	! 1.2.2 On calcule la likelihood sur chaque segment
    	! A - MANIERE FULL EXPLORATION
    	if (init_full == 1) then 
			write(*,*) "... Computing full likelihood for comparison..."
	    	if (flag_ndata(ij)==0) then 								! Pas d evenement dans le segment
	    		write(*,*) '...No events in partition'
	    	    	Likelihood(ij) = 1
	    		MaxLikelihood(ij) = 0	
	    		do i_mag = a_prop, b_prop
				Like_Mu(i_mag,  :) = 1/Mu_bins
				Like_b( i_mag , :) = 1/B_bins
				Like_Sigma( i_mag , :) = 1/Sigma_bins
			end do	
	    	else 											! Au moins un evenement dans le segment								
			call bMuSigma_full(mag_obs(a_segt:b_segt), Mmin,Sigma_vec, Muvec, bvec, mat1, maxval2)
			
		    Likelihood(ij) = sum(mat1)
		    MaxLikelihood(ij) = maxval2
			!!!!! ATTENTION NE PAS DIVISER PAR LE MAX MAIS PAR LA SOMME
		    mat_MB = sum(mat1, dim=3) 
			mat_MS = sum(mat1, dim=2) 
			mat_BS = sum(mat1, dim=1) 
			do i_mag = a_prop, b_prop				
				Like_b( i_mag , :) =  sum(mat_MB, dim=1)/(B_int * sum(mat1))
				Like_Mu(i_mag,  :) =  sum(mat_MB, dim=2)/(Mu_int * sum(mat1))
				Like_Sigma( i_mag , :) =  sum(mat_BS, dim=1)/(Sigma_int * sum(mat1))
			end do
		end if
		if (Likelihood(ij) == 1) then 
	    		logPb(1) =logPb(1) 
	    	else 
	    		logPb(1) = logPb(1) + log(Likelihood(ij))+ MaxLikelihood(ij)- log(dMdB)
	    	end if
		write(*,*) ".. F",ij,a_segt,b_segt,flag_ndata(ij),Likelihood(ij)/(B_bins*Sigma_bins*Mu_bins),Like_b(2,2),MaxLikelihood(ij)
	
	end if  
	! B - MANIERE RANDOM EXPLORATION
	
	if (flag_ndata(ij)==0) then 
	
	    	write(*,*) '...No events in partition'
    	    	Likelihood(ij) = 1
    		MaxLikelihood(ij) = 0	
    		do i_mag = a_prop, b_prop
			Like_Mu(i_mag,  :) = 1/Mu_bins
			Like_b( i_mag , :) = 1/B_bins
			Like_Sigma( i_mag , :) = 1/Sigma_bins
		end do	
	else
		call bMuSigma_rand(mag_obs(a_segt:b_segt),Mmin, n_tirs, bvec, Muvec, Sigma_vec, &
					mat_out, maxval2, MLike_b, MLike_mu, MLike_sigma)
		Likelihood(ij) = sum(mat_out)
	    	MaxLikelihood(ij) = maxval2
			do i_mag = a_prop, b_prop				
				Like_b( i_mag , :) =  MLike_b
				Like_Mu( i_mag , :) =  MLike_mu
				Like_Sigma( i_mag , :) =  MLike_sigma
			end do
	end if
	
	if (Likelihood(ij) == 1) then 
    		logP(1) =logP(1)
    	else 
    		logP(1) = logP(1) + log(Likelihood(ij))+ MaxLikelihood(ij)- log(dMdB)
    end if
	write(*,*)".. Likelihood :",ij,a_segt,b_segt,flag_ndata(ij),Like_b(2,2),maxval2,logP(1)

  end do

  
  !! 1.3 Premier jet de PERT/BIRTH/DEATH !!
  !-----------------------------------------
  write(*,*) ". First random draw Move/Birth/Death"
  ! reinitialisation des vecteurs temp avec le dernier modele
  flag_unaccepted = 0 
  index_prop = 0
  index_prop_mu = 0
  
  tmp_Nc = Nc
  tmp_prop = prop
  tmp_model_times =  model_times
  tmp_model_iprop =  model_iprop
  tmp_model_idata =  model_idata
  
  tmp_model_ddata =  model_idata
  tmp_model_dprop =  model_iprop
  
  tmp_Like = Likelihood
  tmp_MaxLike = MaxLikelihood

  
  
  ! tir du type de perturbation alpha
  alpha = ran3(ra) 
  if (alpha.lt.1./3) then ! Case we perturb one layer
  	Npert=Npert+1
  	write(*,*) ".. Case Move : selecting discontinuity... "
  	cpert = ceiling(ran3(ra)*(Nc))			! on tire aleatoire une des Nc perturbations
  	tp = model_times(cpert)				! on recupere le temps associe a cette perturbation
  	tc = tp + gasdev(ra) * tstd			! on perturbe selon une gaussienne
	
  	
  	if ((tc.ge.tlim) .and. (tc.le.(tmax-tlim))) then ! sinon do nothing
  		tmp_model_times(cpert) = tc		
  		index_prop_mu = ceiling(T_bins * tp/tmax)		! on retrouve l'index associe a ce temps d'observation random dans le vecteur prop
  		index_prop = ceiling(T_bins * tc/tmax)			! on retrouve l'index associe a ce temps d'observation random dans le vecteur prop
  		tmp_prop(index_prop_mu) = 0
  		tmp_prop(index_prop) = 1
  		tmp_model_iprop(cpert) = index_prop
  		
  		tmp_model_ddata(cpert) = 9999999
  		tmp_model_dprop(cpert) = 9999999
  		idata = 1
    		do while (times(idata) .lt. tc) 
    			index_data = idata
    			idata = idata + 1
    		end do	
    		
    		tmp_model_idata(cpert) = index_data
    		tmp_Nc = Nc 
  		!write(*,*) ".. Perturbation acceptee"
  	else 
  		!write(*,*) ".. Pas de perturbation : conditions aux limites"
  		flag_unaccepted = 0
  	end if 
  	call sort_shell_real(tmp_model_times)
  	call sort_shell(tmp_model_iprop)
  	call sort_shell(tmp_model_idata)
  	
	call sort_shell(tmp_model_dprop)
	call sort_shell(tmp_model_ddata)
	write(*,*) ".. Proposed model Perturbation (new Timegrid) ",tmp_model_iprop(1:tmp_Nc)
	
		
  else if ((alpha.ge.1./3) .and. (alpha.lt.2./3))  then   !Case death 
  	write(*,*) ".. Case Death : selecting discontinuity..."
  	Ndeath = Ndeath +1
	if (Nc/=0) then 
  		cpert = ceiling(ran3(ra)*(Nc))		
  		tc = model_times(cpert)
  		tmp_model_times(cpert) = 9999999
  		tmp_model_idata(cpert) = 9999999
  		tmp_model_iprop(cpert) = 9999999
  		
  		index_prop = ceiling(T_bins * tc/tmax)
  		tmp_prop(index_prop) = 0	! on supprime la discontinuite
  		call sort_shell_real(tmp_model_times)
  		call sort_shell(tmp_model_iprop)
  		call sort_shell(tmp_model_idata)
  		
  		tmp_Nc = Nc-1
  		!write(*,*) ".. Death acceptee"
  	else 
  		!write(*,*) ".. Pas de death : 0 couches "
  		flag_unaccepted = 0
  		tmp_Nc=0
  	end if   		
  	write(*,*) ".. Proposed model Death (new Timegrid)", tmp_model_iprop(1:tmp_Nc)	

  	
  else ! Case birth
  	write(*,*) ".. Case Birth : selecting a new discontinuity... "
  	Nbirth = Nbirth + 1
	cpert = Nc + 1 
	if (cpert .le. max_layers) then 
  		tc = tlim+(ran3(ra)*(tmax-(2*tlim)))
  		!write(*,*) cpert, tc
  		tmp_model_times(cpert) = tc
		index_prop = ceiling(T_bins * tc/tmax)			! on retrouve l'index associe a ce temps d'observation random dans le vecteur prop
  		tmp_prop(index_prop) = 1
  		tmp_model_iprop(cpert) = index_prop
  		
  		idata=1
    		do while (times(idata) .lt. tc) 
    			index_data = idata
    			idata = idata + 1
    		end do	
    		tmp_model_idata(cpert) = index_data

  		call sort_shell_real(tmp_model_times)
  		call sort_shell(tmp_model_iprop)
  		call sort_shell(tmp_model_idata)
		tmp_Nc = Nc +1
  		!write(*,*) ".. Birth acceptee"
  	else 
  		!write(*,*) ".. Pas de birth : conditions aux limites nb de couches"
  		flag_unaccepted = 0
  	end if 
  	write(*,*) ".. Proposed model Birth (New TimeGrid)",tmp_model_iprop(1:tmp_Nc)
  end if
  
  
  !! 1.4 CALCULE LIKELIHOOD DATASET TOTAL ie Nc=0 !!
  !-------------------------------------------------
  write(*,*)". Computing once the total Likelihood (case of no discontinuities)..."
  write(*,*)". may take several seconds..."
  !call bMuSigma_full(mag_obs, Mmin,Sigma_vec, Muvec, bvec, mat1, maxval2)
  
  !write(*,*) ".. Total Likelihood (Nc=0) full gridsearch (long) :						",  Likelihood_all,  MaxLikelihood_all
  call bMuSigma_rand(mag_obs,Mmin, n_tirs, bvec, Muvec, Sigma_vec, &
				mat_out, maxval2, MLike_b, MLike_mu, MLike_sigma)
  write(*,*) ".. Total Likelihood (Nc=0) without full gridsearch (random exploration) :  ", sum(mat1), maxval2
  Likelihood_all = sum(mat_out)
  MaxLikelihood_all = maxval2
  dMdB = n_tirs
  log_all = dlog(Likelihood_all)  + MaxLikelihood_all - log(dMdB)
  do k_prop = 1, T_bins				
	Like_b_all( k_prop , :) =  MLike_b
	Like_Mu_all(k_prop,  :) = MLike_mu
	Like_Sigma_all( k_prop , :) =  MLike_sigma
  end do
   write(*,*) ".. If Likelihood full gridsearch different of likelihood random exploration !! Change exploration parameters!!"
   write(*,*) "< End of Initialization >"
   write(*,*)
   write(*,*)
   
   
   
   
   
   
   
   
   
 
   
   
   
   
!----------------------------------------------------------------------------------------------------
! 2- BEGIN TRANSDIMENSIONAL ITERATIONS
!------------------------------------------------------------------------------------------------------
  write(*,*) "> Begin reversible-jump..."
  write(*,*) "--------------------------"
  write(*,*) ". Total number of iterations : ", itmax
  write(*,*) ". Thinning : ", it_std
  write(*,*)
  write(*,*) "> First iteration : 1"
  write(*,*) ". Computing the likelihood of the new model (for comparison)..."
  c = 0
  flag=0
  
  Nperta=0
  Ndeatha=0
  Nbirtha=0

  dMdB = n_tirs ! case random exploration pour condition MCMC
  New_std=std
  do k = 2, itmax
	logTemp=0
  	if (mod(k,it_std)==0) then ! talk
  		
 		  Acc_pert = real(Nperta)/real(Npert)
 		  write(*,*) ". Old perturbation size",New_std
		  Acc_death = real(Ndeatha)/real(Ndeath)
		  Acc_birth = real(Nbirtha)/real(Nbirth)
  		  write(*,*) ". We change the perturbation size :", k, Acc_pert, Acc_death, Acc_birth
  		  New_std = max(0.5, New_std *  Acc_pert/0.3)	!STD in purcentage
  		  write(*,*) ". New perturbation size", New_std
  		  tstd = tmax/100  * New_std
  	end if
 
  	!! 2.1 On calcule likelihood methode random 
  	!-------------------------------------------------------------------------
  	logTemp = 0
  	tmp_Like =0
  	tmp_MaxLike =0
  	
  	if (tmp_Nc==0) then
  		Likelihood(1) = Likelihood_all
		MaxLikelihood(1) = MaxLikelihood_all
		logTemp = log_all
  	else ! On ne calcule la likelihood que la ou elle a change
		if (alpha.lt.1./3) then
			!write(*,*) 'Step 1: Calcul like death'
			tmp_Nc = tmp_Nc-1
			 if (tmp_Nc==0) then
		  		Likelihood(1) = Likelihood_all
				MaxLikelihood(1) = MaxLikelihood_all
			else
			
				call likelihood_death_rand(mag_obs,Mmin, n_tirs, bvec, Muvec, Sigma_vec, &
								cpert, tmp_Nc, T_bins,  &
								tmp_model_ddata, tmp_model_dprop, &  
								New_Like_Mu, New_Like_b, New_Like_Sigma, &
								acc_Like_Mu, acc_Like_b, acc_Like_Sigma, &
								Likelihood, MaxLikelihood, New_Like, New_MaxLike)
				! Like/MaxLike pour condition MCMC				
				tmp_Like = New_Like
				tmp_MaxLike = New_MaxLike
				! Marginaux
				tmp_Like_b = New_Like_b
				tmp_Like_Mu = New_Like_Mu
				tmp_Like_Sigma = New_Like_Sigma
  			end if
			
			tmp_Nc = tmp_Nc+1
			!write(*,*)
			!write(*,*) 'Step 2: Calcul like birth', tmp_model_idata(1:tmp_Nc), tmp_model_iprop(1:tmp_Nc)

			call likelihood_birth_rand(mag_obs,Mmin, n_tirs, bvec, Muvec, Sigma_vec, &
							cpert, tmp_Nc, T_bins,  &
							tmp_model_idata, tmp_model_iprop, & 
							New_Like_Mu, New_Like_b, New_Like_Sigma, &
			 				acc_Like_Mu, acc_Like_b, acc_Like_Sigma, &
			 				tmp_Like, tmp_MaxLike, New_Like, New_MaxLike)
			! Like/MaxLike pour condition MCMC
			tmp_Like = New_Like
			tmp_MaxLike = New_MaxLike
			! Marginaux
			tmp_Like_b = New_Like_b
			tmp_Like_Mu = New_Like_Mu
			tmp_Like_Sigma = New_Like_Sigma
			
		else if ((alpha.ge.1./3) .and. (alpha.lt.2./3)) then
			call likelihood_death_rand(mag_obs,Mmin, n_tirs, bvec, Muvec, Sigma_vec, &
							cpert, tmp_Nc, T_bins,  &
			 				tmp_model_idata, tmp_model_iprop, & 
			 				New_Like_Mu, New_Like_b, New_Like_Sigma, &
			 				acc_Like_Mu, acc_Like_b, acc_Like_Sigma, &
							Likelihood, MaxLikelihood, New_Like, New_MaxLike)
			! Like/MaxLike pour condition MCMC
			tmp_Like = New_Like
			tmp_MaxLike = New_MaxLike
			! Marginaux
			tmp_Like_b = New_Like_b
			tmp_Like_Mu = New_Like_Mu
			tmp_Like_Sigma = New_Like_Sigma
			
		else
			call likelihood_birth_rand(mag_obs, Mmin,n_tirs,bvec, Muvec, Sigma_vec, &
							cpert, tmp_Nc, T_bins,  &
			 				tmp_model_idata, tmp_model_iprop, & 
			 				New_Like_Mu, New_Like_b, New_Like_Sigma, &
							acc_Like_Mu, acc_Like_b, acc_Like_Sigma, &
			 				Likelihood, MaxLikelihood, New_Like, New_MaxLike)
			!  Like/MaxLike pour condition MCMC
			tmp_Like = New_Like
			tmp_MaxLike = New_MaxLike
			! Marginaux
			tmp_Like_b = New_Like_b
			tmp_Like_Mu = New_Like_Mu
			tmp_Like_Sigma = New_Like_Sigma
		end if 
		logTemp = 0 
		do kj = 1, (tmp_Nc +1)
			if (tmp_Like(kj)==1) then 
				logTemp = logTemp
			else 
				logTemp = logTemp + log(tmp_Like(kj)) + tmp_MaxLike(kj) - log(dMdB) ! New value LogP
			end if
		end do 
    	end if
    	






  	!! 2.2 TEST CONDITION MCMC - Metropolis Hastings
  	!-------------------------------------------------------------------------
    	write(*,*) "----------------------"
		write(*,*) ". Metropolis-Hastings test : new / old max loglikelihood"
		write(*,*) ". Metropolis-Hastings test : ", logTemp, logP(k-1)
    	u = log(ran3(ra))
    	condition = logTemp - logP(k-1)
    	
    	if ((u.gt.condition) .or. (flag_unaccepted ==1) ) then   	! 2.2.A TEST CONDITION MH : REJECTION
    	  								!--------------------------------------
    		write(*,*) ".. Model discarded (old>new) ! " 
    		logP(k) = logP(k-1) ! on garde le modele precedent
    		
    		flag_unaccepted =0
    		Nc_total(k) = Nc
    		do kp = 1, T_bins	
    			Sumprop(kp) = Sumprop(kp) + prop(kp)
   		end do
	
	
	
    	else 								! 2.2.B TEST CONDITION MH : ACCEPT
    									!-------------------------------------
		Nc = tmp_Nc
		logP(k) = logTemp
		prop = tmp_prop		
    		model_times = tmp_model_times
    		model_idata = tmp_model_idata
    		model_iprop = tmp_model_iprop
		
		Likelihood = tmp_Like
		MaxLikelihood = tmp_MaxLike
		
		acc_Like_b = tmp_Like_b
		acc_Like_Mu = tmp_Like_Mu
		acc_Like_Sigma = tmp_Like_Sigma

		if (k .gt. burnin) then					! APRES LE BURNIN ON STACKE LES LIKELIHOODS ACCEPTEES
			!if (N_thin .lt. thinning) then 
			!	N_thin = N_thin+1

			!else 
				N_thin=0
	    			Like_b = Like_b + acc_Like_b
	    			Like_Mu = Like_Mu + acc_Like_Mu
	    			Like_Sigma = Like_Sigma + acc_Like_Sigma
				do kp = 1, T_bins	
	    				Sumprop(kp) = Sumprop(kp) + tmp_prop(kp)
	   			end do
			!end if 
    		end if
    		
    		if (alpha.lt.1./3) then
    			Nperta = Nperta +1
    		else if ((alpha.ge.1./3) .and. (alpha.lt.2./3)) then
    			Ndeatha = Ndeatha+1
    		else
    			Nbirtha = Nbirtha + 1 
    		end if 
    		
    		Nc_total(k) = tmp_Nc
    		
    		write(*,*) ".. Model accepted (old<new) : ", model_times(1:Nc)
    		
    	end if  
    	
    	
    	
    	
    	
    	!! 2.3 ON TIRE UN NOUVEAU MODELE DE DISCONTINUITES
    	!-------------------------------------------------------------------------

  		write(*,*) "< End of iteration >"
    	write(*,*)
		write(*,*)
		if (k.eq.itmax) then
			write(*,*) "> Finalisation :", k
		else
    			write(*,*) "> New iteration : ", k 
		endif	
    	! reinitialisation des variables temporaires
    	flag_unaccepted = 0 
    	index_prop = 0
	index_prop_mu = 0
	
	tmp_Nc = Nc	
	tmp_prop = prop
	tmp_model_times =  model_times
	tmp_model_iprop =  model_iprop
	tmp_model_idata =  model_idata
	
	tmp_model_dprop =  model_iprop
	tmp_model_ddata =  model_idata

	! tir du type de perturbation alpha
	write(*,*) ". Random draw : Birth/Death/Move..."
	alpha = ran3(ra) 
	if (alpha.lt.1./3) then 								! 2.3.A PERTURBATION OF ONE LAYER
	  	write(*,*) ".. Case Move : selecting discontinuity..."						!--------------------------------
		Npert = Npert + 1
		if (Nc/=0) then 
			cpert = ceiling(ran3(ra)*(Nc))						! on tire aleatoire une des Nc perturbations : indice Nc
	  		tp = model_times(cpert)							! on recupere le temps de l'ancienne discontinuite
	  		tc = tp + gasdev(ra) * tstd						! on perturbe le temps autour d'une gaussienne
	  		if ((tc.ge.tlim) .and. (tc.le.(tmax-tlim))  ) then 			! SI le nouveau temps n'est pas trop pres des bords : CL acceptees
		  		tmp_model_times(cpert) = tc
		  		tmp_model_ddata(cpert) = 9999999
		  		tmp_model_dprop(cpert) = 9999999
		  		index_prop_mu = ceiling(T_bins * tp/tmax)			! on retrouve l'index associe a l'ancienne discontinuite dans le vecteur prop
		  		index_prop = ceiling(T_bins * tc/tmax)				! on retrouve l'index associe a la nouvelle discontinuite random dans le vecteur prop
		  		tmp_prop(index_prop_mu) = 0					! ancienne discontinuite supprimee du vecteur prop
		  		tmp_prop(index_prop) = 1					! ajout de la nouvelle discontinuite
		  		tmp_model_iprop(cpert) = index_prop				! update du nouveau modele
		    		idata = 1
		    		index_data =0
		    		do while ((times(idata) .le. tc) .and. ((idata) .lt. ndata_mag)) ! on retrouve l'index associe a ce temps d'observation random dans les donnees
		    			index_data = idata
		    			idata = idata + 1
		    		end do	
		    		tmp_model_idata(cpert) = index_data
		    		tmp_Nc = Nc 
		  		
		  		!write(*,*) ".. Perturbation acceptee"
	  		else  									! sinon do nothing
	  			!write(*,*) ".. Pas de perturbation : conditions aux limites"
	  			flag_unaccepted = 1
	  		end if 
	  	else 
	  		!write(*,*) ".. Pas de perturbation : 0 couches "
	  		flag_unaccepted = 1
	  		tmp_Nc=0
	  	end if
	  	call sort_shell_real(tmp_model_times)	! SORT modele en fonction du temps des discontinuites
	  	call sort_shell(tmp_model_iprop)
	  	call sort_shell(tmp_model_idata)
	  	call sort_shell(tmp_model_dprop)
	  	call sort_shell(tmp_model_ddata)
		write(*,*) ".. Proposed model Move :",tmp_model_iprop(1:tmp_Nc)
	
	
	else if ((alpha.ge.1./3) .and. (alpha.lt.2./3))  then   					! 2.3.B DEATH OF ONE LAYER
	  	write(*,*) ".. Case Death : selecting discontinuity... "								!---------------------------
	  	Ndeath = Ndeath +1
		if (Nc/=0) then 									! s'il y a au moins une discontinuite dans le modele precedent
	  		cpert = ceiling(ran3(ra)*(Nc))							! on tire au hazard une des dscontinuites de l ancien modele
	  		tmp_model_times(cpert) = 9999999
	  		tmp_model_idata(cpert) = 9999999
	  		tmp_model_iprop(cpert) = 9999999
	  		tc= model_times(cpert)								! temps de la discontinuite que l on va supprimer
	  		!write(*,*) cpert, tc
	  		index_prop = ceiling(T_bins *tc/tmax)						! on retrouve l'index associe a la discontinuite dans le vecteur prop	
	  		tmp_prop(index_prop) = 0							! on supprime la discontinuite
	  		call sort_shell_real(tmp_model_times)
	  		call sort_shell(tmp_model_iprop)
	  		call sort_shell(tmp_model_idata)
	  		tmp_Nc = Nc-1
	  		!write(*,*) ".. Death acceptee"
	  	else 
	  		!write(*,*) ".. Pas de death : 0 couches "
	  		flag_unaccepted = 1
	  		tmp_Nc=0
	  	end if   		
	  	write(*,*) ".. Proposed model Death", tmp_model_iprop(1:tmp_Nc)
	  	
	  		
	else 												! 2.3.B BIRTH OF ONE LAYER
	  	write(*,*) ".. Case Birth : selecting new discontinuity... "								! ----------------------------
	  	Nbirth = Nbirth + 1
		cpert = Nc + 1 
		if (cpert .le. max_layers) then 
	  		tc = tlim + (ran3(ra)*(tmax-(2*tlim)))						! on tire aleatoirement le temps de la discontinuite que l on va ajouter
	  		tmp_model_times(cpert) = tc							! ajout de la discontinuite
			index_prop =  ceiling(T_bins * tc/tmax)						! on retrouve l'index associe a ce temps d'observation random dans le vecteur prop
	  		tmp_prop(index_prop) = 1
	  		tmp_model_iprop(cpert) = index_prop
	  		
	    		idata = 1
	    		index_data =0
	    		do while ((times(idata) .le. tc) .and. ((idata) .lt. ndata_mag)) 		! on retrouve l'index associe a ce temps d'observation random dans les donnees
	    			index_data = idata
	    			idata = idata + 1
	    		end do	
	    		tmp_model_idata(cpert) = index_data						! update du modele
	  		call sort_shell_real(tmp_model_times)
	  		call sort_shell(tmp_model_iprop)
	  		call sort_shell(tmp_model_idata)
			tmp_Nc = Nc +1
	  		!write(*,*) ".. Birth acceptee"
	  		do tci = 1, Nc+1					! Recherche la nouvelle place de la perturbation dans le vecteur tmp_model_times
     				if (tmp_model_times(tci) == tc) then
        				cpert = tci
     				end if
     			end do
	  	else 
	  		!write(*,*) ".. Pas de birth : conditions aux limites nb de couches"
	  		flag_unaccepted = 1
	  	end if 
	  	write(*,*) ".. Proposed model Birth :",tmp_model_iprop(1:tmp_Nc)
	end if
	write(*,*) ".. New number of discontinuities :", tmp_Nc	
	write(*,*) ".. --------------------------"
	write(*,*) ".. Recomputing likelihood..."
  
  
  end do  
!----------------------------------------------------------------------------------------------------
! 2- FIN TRANSDIMENSIONAL ITERATIONS
!------------------------------------------------------------------------------------------------------ 
  CALL cpu_time(t2)  !Toc. end counting time
  comp_time = t2-t1
  write(*,*) "< End of reversible-jump >"
  write(*,*) " -------------------------"
  write(*,*) "> End of inversion (resume) "
  write(*,*) ". Total time of inversion in seconds : " , comp_time
  write(*,*)
  
  ! Petit bilan des taux d'acceptance 
  Acc_pert = real(Nperta)/real(Npert)*100
  Acc_death = real(Ndeatha)/real(Ndeath)*100
  Acc_birth = real(Nbirtha)/real(Nbirth)*100
  write(*,*) " Final acceptance rates : "
  write(*,*) " Name : number_accepted  number_proposed    rate(pct)"
  write(*,*) "Births : ", Nbirtha, Nbirth, Acc_birth
  write(*,*) "Deaths : ", Ndeatha, Ndeath, Acc_death
  write(*,*) "Moves  : ", Nperta, Npert, Acc_pert
  

  
  
!----------------------------------------------------------------------------------------------------
! #4- WRITING RESULTS - 
!------------------------------------------------------------------------------------------------------
  
  write(*,*) "..Writting outputs files... "
  !close(51)

   open(unit=60, file="../out/hypercube.txt",  status="replace", action="write")
	write(60,*) Mu_min, Mu_max, Mu_bins 
	write(60,*) B_min, B_max, B_bins
	write(60,*) Sigma_min, Sigma_max, Sigma_bins
	write(60,*) itmax, n_tirs, T_bins, comp_time
  close(60)

  open(unit=62, file="../out/Acceptances.txt",  status="replace", action="write")

  write(62,*) Nbirtha,Ndeatha, Nperta  
  write(62,*) Nbirth, Ndeath,Npert
  write(62,*) Acc_birth,  Acc_death, Acc_pert
  close(62)
  
  
  open(unit=63, file="../out/NC_tot.txt",  status="replace", action="write")
  do line63 = 1, itmax
  	write(63, *) Nc_total(line63)
  end do
  close(63)
  
  open(unit=64, file="../out/prop.txt",  status="replace", action="write")
  do line64 = 1, T_bins
  
  	write(64, *) Sumprop(line64)
  end do
  close(64)
  
  open(unit=55, file="../out/Likelihood_Mu.txt",  status="replace", action="write")
  do line55 = 1, Mu_bins
  	do column55 = 1, T_bins
  		write(55, '(F25.5,2X)', advance='no') Like_Mu(column55,line55)
  	end do
  	write(55,*)
  end do
  close(55)
  
  
  open(unit=56, file="../out/Likelihood_b.txt",  status="replace", action="write")
  do line56 = 1, B_bins
  	do column56 = 1, T_bins
  		write(56, '(F25.5,2X)', advance='no') Like_b(column56,line56)
  	end do
  	write(56,*)
  end do
  close(56)
  
  open(unit=59, file="../out/Likelihood_Sigma.txt",  status="replace", action="write")
  do line59 = 1, Sigma_bins
  	do column59 = 1, T_bins
  		write(59, '(F25.5,2X)', advance='no') Like_Sigma(column59,line59)
  	end do
  	write(59,*)
  end do
  close(59)
  write(*,*) ".. Completed"
  write(*,*) "< End of b-Bayesian program > "
  

  !------------------------------------------------------------------------------------------------------

  
end program b_Bayesian










!-----------------------------------------
! UTILS
!-----------------------------------------


 
 FUNCTION ran3(idum)
INTEGER idum
INTEGER MBIG,MSEED,MZin
!     REAL MBIG,MSEED,MZ
REAL ran3,FAC
PARAMETER (MBIG=800000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
!     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
INTEGER i,iff,ii,inext,inextp,k
INTEGER mj,mk,ma(55)
!     REAL mj,mk,ma(55)
SAVE iff,inext,inextp,ma
DATA iff /0/
! write(*,*)' idum ',idum
if(idum.lt.0.or.iff.eq.0)then
	iff=1
	mj=MSEED-iabs(idum)
	mj=mod(mj,MBIG)
	ma(55)=mj
	mk=1
	do 11 i=1,54
	ii=mod(21*i,55)
	ma(ii)=mk
	mk=mj-mk
	if(mk.lt.MZ)mk=mk+MBIG
	mj=ma(ii)
!  write(*,*)' idum av',idum
11      continue
	do 13 k=1,4
	do 12 i=1,55
	ma(i)=ma(i)-ma(1+mod(i+30,55))
	if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
! write(*,*)' idum ap',idum
	inext=0
	inextp=31
	idum=1
endif
! write(*,*)' idum app ',idum
inext=inext+1
if(inext.eq.56)inext=1
inextp=inextp+1
if(inextp.eq.56)inextp=1
mj=ma(inext)-ma(inextp)
if(mj.lt.MZ)mj=mj+MBIG
ma(inext)=mj
ran3=mj*FAC
!  write(*,*)' idum ',idum
	
return
END


FUNCTION GASDEV(idum)

!     ..Arguments..
integer          idum
real GASDEV

!     ..Local..
real v1,v2,r,fac
real ran3

if (idum.lt.0) iset=0
10 v1=2*ran3(idum)-1
v2=2*ran3(idum)-1
r=v1**2+v2**2
if(r.ge.1.or.r.eq.0) GOTO 10
fac=sqrt(-2*log(r)/r)
GASDEV=v2*fac

RETURN
END











































		
	
	

