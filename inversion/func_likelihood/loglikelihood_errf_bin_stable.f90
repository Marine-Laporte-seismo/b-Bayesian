subroutine LogLikelihood_function_bin_detect(ibin_obs, N_data_subset, N_bins, delta_bin, min_bin, &
                                             beta_i, mu_i, sigma_i, LogLikelihood)
! C ----------------------------------------------------------------------
! C   Compute the log-likelihood of ONE model (beta, mu, sigma) for a subset of magnitudes
! C   considering a detection function (Error function) and a Gutenberg Richter law,
! C   with BINNED magnitudes (same formula as Likelihood_FMD_gridsearch)
! C   Inputs :
! C    . ibin_obs      : magnitude bin index (1..N_bins) of each event of the subset {integer array}
! C    . N_data_subset : number of events in the subset {integer}
! C    . N_bins        : number of magnitude bins {integer}
! C    . delta_bin     : size of the magnitude binning (usually 0.1) {float}
! C    . min_bin       : lower edge of the first bin = lower integration bound {float}
! C    . beta_i, mu_i, sigma_i : proposed model
! C   Output :
! C    . LogLikelihood : log-likelihood of the model (always finite)
! C ----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: ibin_obs(:)
  integer, intent(in) :: N_data_subset, N_bins
  real*8, intent(in) :: delta_bin
  real*8, intent(in) :: min_bin
  real*8, intent(in) :: beta_i, mu_i, sigma_i
  real*8, intent(out) :: LogLikelihood
 
  ! Histogram of the subset
  integer :: mag_bins(N_bins)
  real*8 :: Mmin
 
  ! Detection function
  real*8 :: Q_Mmin, Q_Mc
  real*8 :: Phi_cste, Phi_bin
  real*8 :: terme_m1, terme_m2, terme_3
  real*8 :: m1, m2, q_m1, q_m2, m1_t3, m2_t3, q_m1_t3, q_m2_t3
 
  ! Likelihood/Loglikelihood
  real*8 :: Likelihood
 
  integer :: i, k ! iteratives
 
  Mmin = min_bin
 
  ! Step 0 - Histogram of the magnitudes of the subset
  mag_bins = 0
  do i = 1, N_data_subset
     mag_bins(ibin_obs(i)) = mag_bins(ibin_obs(i)) + 1
  end do
 
  !%----Integral Bounds for q error detection function
  Q_Mmin = q_err(Mmin)
  Q_Mc = q_err(Mmin + sigma_i**2 * beta_i)
 
  Phi_cste = Q_Mmin + exp( ((sigma_i**2) * (beta_i**2) /2) - (beta_i*(mu_i-Mmin)))*(1-Q_Mc)
 
  LogLikelihood = 0.0d0
 
  ! Step 1 - Loop over the number of bins
  do k = 1, N_bins
     !%---- Likelihood of the frequency-magnitude distribution
     !%---- GR(m) x Q(m)
     if (mag_bins(k) .ne. 0) then
        !%- terme M1
        m1 = Mmin + delta_bin * (k-1)
        q_m1 = q_err(m1)
        terme_m1 =  exp(-beta_i*(m1-Mmin))*q_m1
 
        !%- terme M2
        m2 = Mmin + delta_bin * k
        q_m2 = q_err(m2)
        terme_m2 =  exp(-beta_i*(m2-Mmin))*q_m2
 
        !%- terme 3
        m1_t3 = m1 + (sigma_i**2 * beta_i)
        q_m1_t3 = q_err(m1_t3)
        m2_t3 = m2 + (sigma_i**2 * beta_i)
        q_m2_t3 = q_err(m2_t3)
 
        terme_3 = exp(((sigma_i**2)*(beta_i**2))/2 - (beta_i*(mu_i-Mmin))) *( q_m2_t3 - q_m1_t3)
 
        !%---- Likelihood formula (from Marsan private comm see Laporte et al.,2025)
        Phi_bin =  - terme_m2 + terme_m1 + terme_3
 
        if ((Phi_bin > 1.0d-12*max(terme_m1, terme_m2, terme_3)) .and. (Phi_bin > tiny(1.0d0))) then
           Likelihood = Phi_bin/(Phi_cste*delta_bin)
           LogLikelihood = LogLikelihood + mag_bins(k)* log(Likelihood)
        else
           !%--- Bin far in the lower tail of the detection : the formula loses its precision.
           !%--- Approximation by the density at the center of the bin (always finite)
           LogLikelihood = LogLikelihood + mag_bins(k)* ( log(beta_i) - beta_i*(0.5d0*(m1+m2)-Mmin) &
                                                         + log_q_err(0.5d0*(m1+m2)) - log(Phi_cste) )
        end if
     end if
  end do
  return
 
contains
 
  ! q(m) = Phi((m-mu)/sigma), computed with erfc : precise even when q is very small
  real*8 function q_err(m)
    real*8, intent(in) :: m
    q_err = 0.5d0 * erfc( -(m-mu_i)/(sqrt(2.0d0)*sigma_i) )
  end function q_err
 
  ! log(q(m)), finite even far in the lower tail
  real*8 function log_q_err(m)
    real*8, intent(in) :: m
    real*8 :: x
    x = -(m-mu_i)/(sqrt(2.0d0)*sigma_i)
    if (x > 0.0d0) then
       log_q_err = log(0.5d0*erfc_scaled(x)) - x*x
    else
       log_q_err = log(0.5d0*erfc(x))
    end if
  end function log_q_err
 
end subroutine LogLikelihood_function_bin_detect


