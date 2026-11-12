subroutine bMuSigma_rand(mag_obs, Mmin, n_tirs, bvec, Muvec, Sigmavec, &
				mat_out, maxval2, MLike_b, MLike_mu, MLike_sigma)
! C ----------------------------------------------------------------------
! C   Compute the posterior of the fitted magnitude distribution 
! C   considering a detection function (Gumbel) and a Gutenberg Richter law
! C   Inputs : 
! C    . M_data : list of magnitudes {array or list}
! C    . priors : optional, to choose if we impose a prior on the Gumbel function 
! C    Returns :
! C    . prob : 3D posterior distribution for the 3 parameters of the detection function (mu, sigma) 
! C    and of the Gutenberg Richter law
! C    . maxlike : Maximum of the 3D posterior distribution with the best three parameters values"""
! C ----------------------------------------------------------------------
  implicit none
  real*8 , EXTERNAL  :: gammp, ran3, gasdev
  
  real*8, intent(in) :: mag_obs(:)
  real*8, intent(in) :: Mmin
  integer, intent(in):: n_tirs
  real*8, intent(in) :: bvec(:), Muvec(:), Sigmavec(:)
  real*8, intent(out) :: mat_out(n_tirs)
  real*8, intent(out) :: maxval2
  real*8, intent(out) :: MLike_b(size(bvec)), MLike_mu(size(Muvec)), MLike_sigma(size(Sigmavec))
  
  real*8:: B_min, B_max, Mu_min, Mu_max, Sigma_min, Sigma_max, B_int, Mu_int, Sigma_int
  real*8 :: mat_in(n_tirs, 3)
  real*8 :: sigma_g, b_g, mu_g
  !real*8 :: gamma_val		 ! GUMBEL FUNCTION
  !real*8 :: c1, c2, c3, c4, croch	 ! SIMPLE EXP FUNCTION
  real*8 :: Q_Mmin, Q_Mc	 ! ERF FUNCTION
  real*8 :: Phi_cste
  real*8 :: q(size(mag_obs))
  real*8 :: maxvalb, maxvalmu, maxvalsigma
  integer :: N, Nb, Nm, Ns
  integer :: ip, iq,  k, ib, im, is, ic
  integer :: indices_b(n_tirs), indices_mu(n_tirs), indices_sigma(n_tirs)
  integer :: counts_b(size(bvec)), counts_mu(size(Muvec)), counts_sigma(size(Sigmavec))
  real*8 :: Like_GG, LogLike_GG
  real*8 :: M_min
  
  ! Init params for random exploration
  N = size(mag_obs) 
  Nb = size(bvec)
  Nm = size(Muvec)
  Ns = size(Sigmavec)
  
  B_min = minval(bvec)
  B_max = maxval(bvec)
  Mu_min = minval(Muvec)
  Mu_max = maxval(Muvec)
  Sigma_min = minval(Sigmavec)
  Sigma_max = maxval(Sigmavec)
  B_int = (B_max-B_min)/n_tirs
  Mu_int = (Mu_max-Mu_min)/n_tirs
  Sigma_int = (Sigma_max-Sigma_min)/n_tirs
  
  M_min = Mmin
  ! Init -Marginals
  MLike_b = 0.0
  MLike_mu = 0.0
  MLike_sigma = 0.0
  ! Sample over the 3 directions
  call random_number(mat_in)
  !write(*,*) "OK : exemple val b", mat_in(3,1)
  do ip = 1, n_tirs
	b_g     = mat_in(ip,1)* (B_max-B_min) + B_min
	mu_g    = mat_in(ip,2)* (Mu_max-Mu_min) + Mu_min
	sigma_g = mat_in(ip,3)* (Sigma_max-Sigma_min) + Sigma_min
	
	LogLike_GG = 0
	Like_GG = 0  
	! ERROR FUNCTION
        !-----------------------------------------------------------------------------------------
	Q_Mmin = 0.5 + 0.5 * erf( (M_min-mu_g) / (sqrt(2.0)*sigma_g) )
	Q_Mc = 0.5 + 0.5 * erf( ((M_min+sigma_g**2 * b_g)- mu_g) / (sqrt(2.0)*sigma_g) )
	Phi_cste = Q_Mmin + exp( ((sigma_g**2) * (b_g**2) /2) - (b_g*(mu_g-M_min)))*(1-Q_Mc)
	do k = 1, N
		!%---- q error function		
	  	q(k) = 0.5 + 0.5 * erf( (mag_obs(k)-mu_g)/(sqrt(2.0)*sigma_g))		
	        Like_GG =  q(k)   * b_g * exp(-b_g * (mag_obs(k)-M_min))  /   Phi_cste			
    	        LogLike_GG = LogLike_GG + dlog( Like_GG )
	end do

	! GUMBEL FUNCTION
        !-----------------------------------------------------------------------------------------
        !gamma_val = gammp( b_g * sigma_g ,exp(-(Mmin-mu_g)/sigma_g)) 
	!----integral for q Gumbel function (see David's paper)
        !Phi_cste = b_g * sigma_g * gamma(b_g*sigma_g) * exp(b_g*(Mmin-mu_g)) * gamma_val

   	!do k = 1, N ! # Instead of product over different likelihood for each M(i), take the sum of the log
        !        !---- q Gumbel function
	!	q(k)= exp(- exp(-(mag_obs(k)-mu_g)/sigma_g))
	!	Like_GG = q(k)  * b_g * exp(- b_g * (mag_obs(k)-Mmin)) / Phi_cste
	!	LogLike_GG = LogLike_GG  + dlog(Like_GG)
   	!end do
	! EXPONENTIAL FUNCTION
        !-----------------------------------------------------------------------------------------
	!%----integral for qsimple exponential
	!if (0 < mu_g) then
	!	Phi_cste = (1/(1+b_g*sigma_g)**3) * exp(-b_g*(mu_g-Mmin))       		
	!else 
	!    	c1 = 1/(1+sigma_g*b_g) 
	!    	c2 = 1/(1+sigma_g*b_g)**2 
	!    	c3 = ((Mmin-mu_g)/sigma_g) *(1 + (1/(1+sigma_g*b_g)))
	!    	c4 = (Mmin - mu_g)**2/(2*sigma_g**2)
	!       croch = 1 + c1 +  c2  - c3 + c4
	!       Phi_cste = 1-( b_g*sigma_g*exp(mu_g/sigma_g) / (1+b_g*sigma_g) * croch)
	        
	 !end if
	 !do k = 1, N
	 !       !%---- q simple exponential
	 !       if (mag_obs(k) > mu_g) then
	 !           	q(k) = 1 - exp(-(mag_obs(k) - mu_g)/sigma_g) * (1 + (mag_obs(k)-mu_g)/sigma_g + (mag_obs(k)-mu_g)**2/(2*sigma_g**2))		
	 !           	Like_GG =  q(k)   * b_g * exp(-b_g * (mag_obs(k)-Mmin))  /   Phi_cste		
    	 !       	LogLike_GG = LogLike_GG + dlog( Like_GG )
	 !       else
	 !           	q(k)=0
	 !           	Like_GG=0
	 !           	LogLike_GG=-99999
	 !      	end if
	 !end do	
	

	mat_out(ip) = LogLike_GG
	! For reconstructing the marginals afterward
	! Keeping in mind the values explored
	indices_b(ip) = ceiling(mat_in(ip,1)*100)
	indices_mu(ip) = ceiling(mat_in(ip,2)*100)
	indices_sigma(ip) = ceiling(mat_in(ip,3)*100)
    end do
    
    
    maxval2 = maxval(mat_out)
    mat_out = exp(mat_out - maxval2)
    ! Interpolation of B-Mu-Sigma marginals - projection into the respective arrays
    do iq = 1, n_tirs
    	    ib = indices_b(iq) 
    	    im = indices_mu(iq) 
    	    is = indices_sigma(iq) 
	    MLike_b(ib) =  MLike_b(ib) + mat_out(iq)
	    MLike_mu(im) = MLike_mu(im) + mat_out(iq)
	    MLike_sigma(is) = MLike_sigma(is) + mat_out(iq)
    end do
    
    
    ! Double step of normalisation of marginals
    !! Step one : Divides by the number of occurrence in each bin
    do ib =  1, size(bvec) 
    	counts_b(ib) = count(indices_b == ib)
    	MLike_b(ib) =  MLike_b(ib) / counts_b(ib)
    end do 
    
    do im = 1, size(Muvec)
    	counts_mu(im) = count(indices_mu == im)
    	MLike_mu(im) = MLike_mu(im)/counts_mu(im)
    end do
    
    do is = 1, size(Sigmavec)
    	counts_sigma(is) = count(indices_sigma == is)
    	MLike_sigma(is) = MLike_sigma(is)/counts_sigma(is)
    end do
    
    !! Step two : Divides by the maximum likelihood 
    maxvalb = maxval(MLike_b)                 
    maxvalmu = maxval(MLike_mu)           
    maxvalsigma = maxval(MLike_sigma)
    MLike_b = MLike_b/ (B_int * sum(MLike_b))
    MLike_mu = MLike_mu/(Mu_int * sum(MLike_mu))
    MLike_sigma = MLike_sigma/(Sigma_int * sum(MLike_sigma))
    return 
end subroutine   bMuSigma_rand




subroutine likelihood_death_rand(mag_obs,Mmin,  n_tirs, bvec, Muvec, Sigmavec, cpert, tmp_Nc,T_bins, &
				      tmp_model_idata, tmp_model_iprop, & 
				      New_Like_Mu, New_Like_b, New_Like_Sigma, &
				      acc_Like_Mu, acc_Like_b, acc_Like_Sigma, &
				      Likelihood, MaxLikelihood, New_Like, New_MaxLike)
	interface				
	  subroutine bMuSigma_rand(mag_obs, Mmin, n_tirs, bvec, Muvec, Sigmavec, &
				        mat_out, maxval2, MLike_b, MLike_mu, MLike_sigma)
		  real*8, intent(in) :: mag_obs(:)
		  real*8, intent(in) :: Mmin
		  integer, intent(in) :: n_tirs
		  real*8, intent(in) :: bvec(:), Muvec(:), Sigmavec(:)
		  real*8, intent(out) :: mat_out(n_tirs)
		  real*8, intent(out) :: maxval2
		  real*8, intent(out) :: MLike_b(size(bvec)), MLike_mu(size(Muvec)), MLike_sigma(size(Sigmavec))
	  end subroutine bMuSigma_rand 
	end interface
	
	real*8, intent(in) :: mag_obs(:)
	real*8, intent(in) :: Mmin
	integer, intent(in) :: n_tirs
	real*8, intent(in) :: bvec(:), Muvec(:), Sigmavec(:)
	integer,intent(in) :: cpert, tmp_Nc, T_bins
	integer ,intent(in):: tmp_model_idata(:), tmp_model_iprop(:)
	
	real*8 ,intent(in):: Likelihood(:), MaxLikelihood(:)
	real*8 ,intent(in):: acc_Like_Mu(:,:), acc_Like_b(:,:), acc_Like_Sigma(:,:)
	real*8, intent(out):: New_Like(:), New_MaxLike(:)
	real*8, intent(out) :: New_Like_Mu(:,:), New_Like_b(:,:),New_Like_Sigma(:,:)
	
	integer :: cc,n, k_prop
	integer :: Nm, Nb, Ns
	integer :: a_death, b_death, a_prop, b_prop
	real*8 :: mat_out(n_tirs)
	
	real*8 :: maxval2
	real*8 :: MLike_b(size(bvec)), MLike_mu(size(Muvec)), MLike_sigma(size(Sigmavec))
	integer :: flag_ndata(tmp_Nc+1)

	New_Like = 0
	NewMaxLike = 0 

	New_Like_b = 0
	New_Like_Mu =0
	New_Like_Sigma = 0
	
	
	
	flag_ndata = 0
	Nm = size(Muvec(:))
	Ns=size(Sigmavec(:))
	Nb = size(bvec(:))
	write(*,*) ".. Death"
	write(*,*) ".. Id of the concerned discontinuity:", cpert
	
	do cc = 1, (tmp_Nc+1)
		! write(*,*) tmp_model_idata(1:tmp_Nc)	
	      	if (cc==1) then								! 0 - Find indexes
	    		flag_ndata(1) = tmp_model_idata(1) 
	    		
	    		a_death = 1 
	      		b_death = a_death + flag_ndata(1) - 1
	      		
	      		a_prop = 1
	      		b_prop = tmp_model_iprop(cc)

    		else if (cc == tmp_Nc+1) then
    			flag_ndata(cc) = size(mag_obs) - tmp_model_idata(cc-1)  
    		
    			a_death = a_death + flag_ndata(cc-1)
	    		b_death = a_death + flag_ndata(cc) - 1
    		
    			a_prop = tmp_model_iprop(cc-1)+1
      			b_prop = T_bins

    		else
    			flag_ndata(cc) = tmp_model_idata(cc) - tmp_model_idata(cc-1) 
    	
    			a_death = a_death + flag_ndata(cc-1)
	    		b_death = a_death + flag_ndata(cc) - 1

	      		a_prop = tmp_model_iprop(cc-1)+1
	      		b_prop = tmp_model_iprop(cc)
    		end if



		if (cc .lt. cpert) then  							 ! 1 - Avant la discontinuite  CC<CPERT
		
			New_Like(cc) = Likelihood(cc)
			New_MaxLike(cc) = MaxLikelihood(cc)
			
			do k_prop = a_prop, b_prop
					New_Like_b( k_prop , :) =  acc_Like_b(k_prop,  :)
					New_Like_Mu(k_prop,  :) =  acc_Like_Mu(k_prop,  :)
					New_Like_Sigma( k_prop , :) =  acc_Like_Sigma(k_prop,  :)
			end do
		
		else if (cc==cpert) then 							! 2 - A l endroit de la discontinuite CC= CPERT
			if (flag_ndata(cc)==0) then 						! 2.a - Pas d evenement dans le segment
    				write(*,*) '.. No events in the new partition'
    	    			New_Like(cc) = 1
    				New_MaxLike(cc) = 0
    				do k_prop = a_prop, b_prop
						New_Like_Mu(k_prop,  :) =  1/Nm
						New_Like_b( k_prop , :) =  1/Nb
						New_Like_Sigma( k_prop , :) =  1/Ns
				end do	
			else 									! 2.b - Evenements dans le segment
				write(*,*) "... We recompute the likelihood on the segment (Timegrid)",  a_prop, b_prop
				call bMuSigma_rand(mag_obs(a_death:b_death), Mmin, n_tirs, bvec, Muvec, Sigmavec, &
							mat_out, maxval2, MLike_b, MLike_mu, MLike_sigma)
				New_Like(cc) = sum(mat_out)
				New_MaxLike(cc) = maxval2
				


				do k_prop = a_prop, b_prop
					New_Like_b( k_prop , :) = MLike_b
					New_Like_Mu(k_prop,  :) =  MLike_mu
					New_Like_Sigma( k_prop , :) = MLike_sigma
				end do
			end if
			
		else if (cc .gt.cpert) then							! Après la discontinuite
			New_Like(cc) = Likelihood(cc+1)
			New_MaxLike(cc) = MaxLikelihood(cc+1)
			do k_prop = a_prop, b_prop
					New_Like_Mu(k_prop,  :) =  acc_Like_Mu(k_prop,  :)
					New_Like_b( k_prop , :) =  acc_Like_b(k_prop,  :)
					New_Like_Sigma( k_prop , :) =  acc_Like_Sigma(k_prop,  :)
			end do
		end if
		write(*,*) "... Likelihood of the partition", cc, a_prop, b_prop, New_Like(cc)
	end do
	write(*,*) "... Likelihood of the partition :     id     begin(Tgrid) 	end(Tgrid)     New MaxLikelihood"
	if (tmp_Nc+3 >= size(Likelihood)) then 
		New_Like(tmp_Nc+2) = 999999
	else 
		New_Like(tmp_Nc+2) = Likelihood(tmp_Nc+3) ! On associe 999999 pour le segment disparu	
	end if 
	return
	
	
end subroutine likelihood_death_rand





subroutine likelihood_birth_rand(mag_obs, Mmin, n_tirs, bvec, Muvec, Sigmavec, cpert, tmp_Nc,T_bins, &
				 tmp_model_idata, tmp_model_iprop, & 
				 New_Like_Mu, New_Like_b, New_Like_Sigma, &
				 acc_Like_Mu, acc_Like_b, acc_Like_Sigma, &
				 Likelihood, MaxLikelihood, New_Like, New_MaxLike)
	interface
	  subroutine bMuSigma_rand(mag_obs,Mmin, n_tirs, bvec, Muvec, Sigmavec, &
				        mat_out, maxval2, MLike_b, MLike_mu, MLike_sigma)
		  real*8, intent(in) :: mag_obs(:)
		  real*8, intent(in) :: Mmin
		  integer, intent(in) :: n_tirs
		  real*8, intent(in) :: bvec(:), Muvec(:), Sigmavec(:)
		  real*8, intent(out) :: mat_out(n_tirs)
		  real*8, intent(out) :: maxval2
		  real*8, intent(out) :: MLike_b(size(bvec)), MLike_mu(size(Muvec)), MLike_sigma(size(Sigmavec))
	  end subroutine bMuSigma_rand 
	end interface
	
	
	real*8, intent(in) :: mag_obs(:)
	real*8, intent(in) :: Mmin
	integer, intent(in) :: n_tirs
	real*8, intent(in) :: bvec(:), Muvec(:), Sigmavec(:)
	integer,intent(in) :: cpert, tmp_Nc, T_bins
	integer ,intent(in):: tmp_model_idata(:), tmp_model_iprop(:)
	
	real*8 ,intent(in):: Likelihood(:), MaxLikelihood(:)
	real*8 ,intent(in):: acc_Like_Mu(:,:), acc_Like_b(:,:), acc_Like_Sigma(:,:)
	real*8, intent(out):: New_Like(:), New_MaxLike(:)
	real*8, intent(out) :: New_Like_Mu(:,:), New_Like_b(:,:),New_Like_Sigma(:,:)
	
	integer :: cc,n, k_prop
	integer :: Nm, Nb, Ns
	integer :: a_birth, b_birth, a_prop, b_prop
	real*8 :: mat_out(n_tirs)
	real*8 :: maxval2
	integer :: flag_ndata(tmp_Nc+1)
	real*8 :: MLike_b(size(bvec)), MLike_mu(size(Muvec)), MLike_sigma(size(Sigmavec))
	
	
	New_Like = 0
	New_MaxLike = 0 
	
	New_Like_b = 0
	New_Like_Mu = 0
	New_Like_Sigma = 0
	Nm = size(Muvec(:))
	Ns=size(Sigmavec(:))
	Nb = size(bvec(:))
	flag_ndata=0
	write(*,*) ".. Birth"
	write(*,*) ".. Id of the concerned discontinuity:", cpert
	do cc = 1 , (tmp_Nc+1)
	
		if (cc==1) then
	    		flag_ndata(1) = tmp_model_idata(1) 
	    		
	    		a_birth = 1 
	      		b_birth = a_birth + flag_ndata(1) - 1
	      		
	      		a_prop = 1
	      		b_prop = tmp_model_iprop(cc)

    		else if (cc == tmp_Nc+1) then
    			flag_ndata(cc) = size(mag_obs) - tmp_model_idata(cc-1)  
    		
    			a_birth = a_birth + flag_ndata(cc-1)
	    		b_birth = a_birth + flag_ndata(cc) - 1
    		
    			a_prop = tmp_model_iprop(cc-1)+1
      			b_prop = T_bins

    		else
    			flag_ndata(cc) = tmp_model_idata(cc) - tmp_model_idata(cc-1) 
    	
    			a_birth = a_birth + flag_ndata(cc-1)
	    		b_birth = a_birth + flag_ndata(cc) - 1

	      		a_prop = tmp_model_iprop(cc-1)+1
	      		b_prop = tmp_model_iprop(cc)
    		end if
	
		if (cc .lt. cpert) then   										! 1 - Avant la discontinuite
		
			New_Like(cc) = Likelihood(cc)
			New_MaxLike(cc) = MaxLikelihood(cc)
			
			do k_prop = a_prop, b_prop
					New_Like_b( k_prop , :) =  acc_Like_b(k_prop,  :)
					New_Like_Mu(k_prop,  :) =  acc_Like_Mu(k_prop,  :)
					New_Like_Sigma( k_prop , :) =  acc_Like_Sigma(k_prop,  :)
			end do
		
		else if (cc==cpert)  then 
			!write(*,*) 'CC==CPERT',flag_ndata(cc)								! 2- A l endroit de la discontinuite
			if (flag_ndata(cc)==0) then 									! 2.a - Pas d evenement dans le segment
    				!write(*,*) 'Pas d evenements dans segment'
    	    			New_Like(cc) = 1
    					New_MaxLike(cc) = 0
    				do k_prop = a_prop, b_prop
						New_Like_Mu(k_prop,  :) =  1/Nm
						New_Like_b( k_prop , :) =  1/Nb
						New_Like_Sigma( k_prop , :) =  1/Ns
					end do
    			else 												! 2;b - Au moins un evenement dans le segment
				write(*,*) "... We recompute the likelihood on the segment (Timegrid) ",  a_prop, b_prop
				
				!call bMu_func(mag_obs(a_birth:b_birth), Sigma_Gb,  Muvec, bvec,mat1, maxval2 )
				call bMuSigma_rand(mag_obs(a_birth:b_birth),Mmin, n_tirs, bvec, Muvec, Sigmavec, &
							mat_out, maxval2, MLike_b, MLike_mu, MLike_sigma)
				New_Like(cc) = sum(mat_out)
				New_MaxLike(cc) = maxval2
				do k_prop = a_prop, b_prop
					New_Like_b( k_prop , :) = MLike_b/maxval(MLike_b)
					New_Like_Mu(k_prop,  :) =  MLike_mu/maxval(MLike_mu)
					New_Like_Sigma( k_prop , :) = MLike_sigma/maxval(MLike_sigma)
				end do
			end if
		
		
		
		else if (cc==cpert+1) then 
			!write(*,*) 'CC==CPERT+1',flag_ndata(cc)
			if (flag_ndata(cc)==0) then ! Pas d evenement dans le segment
    				write(*,*) '.. No events in the new partition'
    	    			New_Like(cc) = 1
    				New_MaxLike(cc) = 0
    				do k_prop = a_prop, b_prop
						New_Like_Mu(k_prop,  :) =  1/Nm
						New_Like_b( k_prop , :) =  1/Nb
						New_Like_Sigma( k_prop , :) =  1/Ns
				end do
    			else ! Au moins un evenement dans le segment
    
				write(*,*) "... We recompute the likelihood on the segment (Timegrid) :", a_prop, b_prop			
				call bMuSigma_rand(mag_obs(a_birth:b_birth), Mmin, n_tirs, bvec, Muvec, Sigmavec, &
							mat_out, maxval2, MLike_b, MLike_mu, MLike_sigma)
				
				New_Like(cc) = sum(mat_out)
				New_MaxLike(cc) = maxval2
				
				!write(*,*) New_Like(cc)
				do k_prop = a_prop, b_prop
						New_Like_Mu(k_prop,  :) =  MLike_mu
						New_Like_b( k_prop , :) =  MLike_b
						New_Like_Sigma( k_prop , :) =  MLike_sigma
				end do
			end if	
			
		else if (cc .gt.cpert+1) then	! après la discontinuite
			!write(*,*) 'CC>CPERT+1'
			New_Like(cc) = Likelihood(cc-1)
			New_MaxLike(cc) = MaxLikelihood(cc-1)	
			do k_prop = a_prop, b_prop
					New_Like_Mu(k_prop,  :) =  acc_Like_Mu(k_prop,  :)
					New_Like_b( k_prop , :) =  acc_Like_b(k_prop,  :)
					New_Like_Sigma( k_prop , :) =  acc_Like_Sigma(k_prop,  :)
			end do
		end if	

		write(*,*) "... Likelihood of the partition", cc, a_prop, b_prop,New_Like(cc)
	end do
	write(*,*) "... Likelihood of the partition :     id     begin(Tgrid) 	end(Tgrid)     New MaxLikelihood"
	return	
end subroutine likelihood_birth_rand


function gammln(xx)
  !use par_mod, only: dp
  implicit none
  integer :: j
  real*8 :: x,tmp,ser,xx,gammln
  real*8 :: cof(6) = (/ &
       76.18009173, -86.50532033, 24.01409822, &
       -1.231739516, .120858003e-2, -.536382e-5    /)
  real*8 :: stp = 2.50662827465
  real*8 :: half = 0.5, one = 1.0, fpf = 5.5
  x=xx-one
  tmp=x+fpf
  tmp=(x+half)*log(tmp)-tmp
  ser=one
  do j=1,6
    x=x+one
    ser=ser+cof(j)/x
  end do
  gammln=tmp+log(stp*ser)
end function gammln

