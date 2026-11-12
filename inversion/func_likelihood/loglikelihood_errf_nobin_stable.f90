subroutine LogLikelihood_function_nobin_detect(mag_obs, N_data_subset, min_mags, &
                                               beta_prop, mu_prop, sigma_prop, LogLikelihood)
! Log-vraisemblance GR x detection - No binning
! f(m) = Phi((m-mu)/sigma) * beta * exp(-beta(m-Mmin)) / K ,  m >= Mmin
  implicit none
  real*8,  intent(in)  :: mag_obs(:)
  integer, intent(in)  :: N_data_subset
  real*8,  intent(in)  :: min_mags, beta_prop, mu_prop, sigma_prop
  real*8,  intent(out) :: LogLikelihood

  real*8  :: s2, zc, d0, K1, K2, logK, sum_logq, sum_dm
  integer :: i

  s2 = sqrt(2.0d0) * sigma_prop

  ! ---- Normalisation K = K1 + K2 -------------------------------------
  ! K1 = Phi((Mmin-mu)/sigma)
  K1 = 0.5d0 * erfc( -(min_mags - mu_prop) / s2 )
  ! K2 = exp(b^2 s^2/2 - b(mu-Mmin)) * [1 - Phi((Mmin + b s^2 - mu)/s)]
  zc = (min_mags + beta_prop*sigma_prop**2 - mu_prop) / s2
  if (zc > 0.0d0) then
     ! forme stable : exp(a)*erfc(zc) = exp(-d0^2)*erfc_scaled(zc)
     d0 = (min_mags - mu_prop) / s2
     K2 = 0.5d0 * erfc_scaled(zc) * exp(-d0**2)
  else
     K2 = 0.5d0 * erfc(zc) * exp(0.5d0*(beta_prop*sigma_prop)**2 &
                                 - beta_prop*(mu_prop - min_mags))
  end if
  logK = log(K1 + K2)

  ! ---- Somme sur les evenements ---------------------------------------
  sum_logq = 0.0d0
  sum_dm   = 0.0d0
  do i = 1, N_data_subset
     sum_logq = sum_logq + log_Phi( (mag_obs(i) - mu_prop) / s2 )
     sum_dm   = sum_dm   + (mag_obs(i) - min_mags)
  end do

  LogLikelihood = sum_logq + N_data_subset*(log(beta_prop) - logK) - beta_prop*sum_dm

contains

  ! log( Phi ) avec x = (m-mu)/(sqrt(2) sigma), precis meme pour x tres negatif
  pure real*8 function log_Phi(x)
    real*8, intent(in) :: x
    if (x < 0.0d0) then
       log_Phi = log(0.5d0 * erfc_scaled(-x)) - x**2
    else
       log_Phi = log(0.5d0 * erfc(-x))
    end if
  end function log_Phi

end subroutine LogLikelihood_function_nobin_detect


