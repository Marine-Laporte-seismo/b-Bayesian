subroutine LogLikelihood_function_nobin_detect(mag_obs, N_data_subset, min_mags, beta_prop, mu_prop, sigma_prop, LogLikelihood) ! for inversion of mu, b and sigma
! C ----------------------------------------------------------------------
! C   Compute the log-likelihood for a proposed beta value of the Gutenberg-Richter law according to a subset of magnitudes
! C   and for proposed mu/sigma values for the detection part of the frequency-magnitude distribution
! C   Inputs :
! C    . mag_obs : list of magnitudes {array or list}
! C    . N_data_subset : size of the subset of magnitudes {integer}
! C    . min_mags : minimum magnitude of the TOTAL dataset {real}
! C    . beta_prop : proposed value for beta= b_value*ln(10) {real}
! C    . mu_prop : proposed value for mu parameter of the detection function  {real}
! C    . sigma_prop : proposed value for sigma parameter of the detection function {real}
! C    Returns :
! C    . LogLikelihood : LogLikelihood for the model beta to explain the sub-dataset (according to the GR law)
! C -----------------------------------
  implicit none
  real*8, intent(in) :: mag_obs(:)
  integer, intent(in) :: N_data_subset
  real*8, intent(in) :: min_mags
  real*8, intent(in) :: beta_prop
  real*8, intent(in) :: mu_prop
  real*8, intent(in) :: sigma_prop
  real*8, intent(out) :: LogLikelihood

  real*8 :: tmp_loglikelihood
  real*8 :: M_min
  integer :: N
  integer :: i
  real*8 :: Q_Mmin, Q_Mc, K
  real*8 :: q
  N = N_data_subset
  M_min = min_mags
  tmp_loglikelihood = 0.0d0
  LogLikelihood = 0.0d0
  ! ERROR FUNCTION
  !-----------------------------------------------------------------------------------------
  Q_Mmin = 0.5 + 0.5 * erf( (M_min-mu_prop) / (sqrt(2.0d0)*sigma_prop) )
  Q_Mc = 0.5 + 0.5 * erf( ((M_min+sigma_prop**2 * beta_prop)- mu_prop) / (sqrt(2.0d0)*sigma_prop) )
  K = Q_Mmin + exp( ((sigma_prop**2) * (beta_prop**2) /2) - (beta_prop*(mu_prop-M_min)))*(1-Q_Mc)
  do i = 1, N
	!%---- q error function
	 q = 0.5 + 0.5 * erf( (mag_obs(i)-mu_prop)/(sqrt(2.0d0)*sigma_prop))
	 !%--- Consider all cases where q = 0
	 if (q .le. 0) then
	 	LogLikelihood = -99999999
	 else
	 	tmp_loglikelihood =  log(q) + log( beta_prop) -beta_prop * (mag_obs(i)-M_min) - log(   K)
    	 	LogLikelihood = LogLikelihood + tmp_loglikelihood
    	 end if
  end do
    return
end subroutine LogLikelihood_function_nobin_detect
