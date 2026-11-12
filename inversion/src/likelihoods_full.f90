subroutine bMuSigma_full(mag_obs,Mmin, Svec, Muvec, bvec, mat1, maxval2)
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
! C    # --------------------------
! C     # Prior for the detection function
! C    # Define the prior distribution for the detection function
! C    # The detection function need to go between p1 and p2, otherwise, the prior
! C    # is zero
! C    # --------------------------
! C ----------------------------------------------------------------------
  implicit none
  real*8 , EXTERNAL    :: gammp

  real*8, intent(in) :: mag_obs(:)
  real*8, intent(in) :: Mmin
  real*8, intent(in) :: Muvec(:)
  real*8, intent(in) :: bvec(:)  
  real*8, intent(in) :: Svec(:)
  real*8, intent(out) :: mat1(size(Muvec), size(bvec),size(Svec))
  real*8, intent(out) :: maxval2
  
  real*8 :: Like_GG, LogLike_GG
  real*8 :: q(size(mag_obs))
  real*8 :: sigma_g, b_g, mu_g
  integer :: N, Nb, Nm, Ns
  integer :: im, ib,is,  k
 
  real*8 :: q1, q2, p1, p2, t1, t2
  integer :: count_g
  
  real*8 :: Phi_cste,  gamma_val

  
  N = size(mag_obs) 
  Nb = size(bvec)
  Nm = size(Muvec)
  Ns = size(Svec)
  !write(*,*) "Size", N

    
  count_g=0;
  ! Sample over the 3 directions
    do im = 1, Nm
    
        do ib = 1, Nb
        	do is = 1, Ns
               		b_g     = bvec(ib)
               		mu_g    = Muvec(im)
               		sigma_g = Svec(is)
               		!write(*,*) 'Test', sigma_g
		       	Like_GG = 0
		       	LogLike_GG = 0
			count_g = count_g + 1    
			   
			gamma_val = gammp( b_g * sigma_g ,exp(-(Mmin-mu_g)/sigma_g)) 
			!----integral for q Gumbel function (see David's paper)
			Phi_cste = b_g * sigma_g * gamma(b_g*sigma_g) * exp(b_g*(Mmin-mu_g)) * gamma_val

			do k = 1, N ! # Instead of product over different likelihood for each M(i), take the sum of the log
				!---- q Gumbel function
				q(k)= exp(- exp(-(mag_obs(k)-mu_g)/sigma_g))
				Like_GG = q(k)  * b_g * exp(- b_g * (mag_obs(k)-Mmin)) / Phi_cste
				LogLike_GG = LogLike_GG  + dlog(Like_GG)
			end do
                	mat1(im, ib,is) = LogLike_GG
                enddo
     	enddo
    enddo
      maxval2 = maxval(mat1) 
      mat1 = exp(mat1 - maxval2)
    return 
end subroutine   bMuSigma_full



