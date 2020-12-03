! ------------------------------------------------------------------------------
SUBROUTINE gmm_em(n, k, x, mu, sig2, pi, llk, max_iter, eps)
! - Variable declarations ------------------------------------------------------
  IMPLICIT NONE
  ! Input
  INTEGER :: n, k
  DOUBLE PRECISION :: x(n)
  DOUBLE PRECISION :: mu(k), sig2(k), pi(k)
  INTEGER :: max_iter
  DOUBLE PRECISION :: eps
  ! Ouput
  DOUBLE PRECISION :: llk
  ! Local variables
  DOUBLE PRECISION :: p(n, k)
  DOUBLE PRECISION :: prev_llk
  INTEGER :: i, j, l
  DOUBLE PRECISION :: tmp
! - Initialization -------------------------------------------------------------
  llk = -9.9D30
! - EM Algorithm ---------------------------------------------------------------
  DO l=1,max_iter
    ! E Step
    DO j=1,k
      CALL dnorm(n, x, mu(j), sig2(j), p(:, j))
      p(:, j) = p(:, j) * pi(j)
    ENDDO
    ! Check convergence
    prev_llk = llk
    llk = 0.0D0
    DO i=1,n
      tmp = sum(p(i, :))
      llk = llk + log(tmp)
      p(i, :) = p(i, :) / tmp
    ENDDO
    llk = llk / n
    IF (llk - prev_llk < eps) THEN
      max_iter = l
      EXIT
    ENDIF
    ! M Step
    DO j=1,k
      tmp = sum(p(:, j))
      mu(j) = sum(p(:, j) * x) / tmp
      sig2(j) = sum(p(:, j) * (x-mu(j))**2) / tmp
      pi(j) = tmp / n
    ENDDO
  ENDDO
! - Exit -----------------------------------------------------------------------
END SUBROUTINE gmm_em
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
SUBROUTINE dnorm(n, x, mu, sig2, d)
  IMPLICIT NONE
  INTEGER :: n
  DOUBLE PRECISION :: x(n), mu, sig2
  DOUBLE PRECISION :: d(n)
  d = (x - mu)**2 / sig2 + log(2.0D0 * 4.0D0*DATAN(1.D0) * sig2)
  d = exp(-0.5D0 * d)
END SUBROUTINE dnorm
! ------------------------------------------------------------------------------