SUBROUTINE acf(n, x, lag, ac)
    INTEGER :: n, lag, k
    DOUBLE PRECISION :: x(n), ac(lag+1)
    DOUBLE PRECISION :: mu, sig
    DOUBLE PRECISION :: xc(n)
    ac(1) = 1.0D0
    mu = sum(x) / n
    xc = x - mu
    sig = sqrt(sum(xc*xc)/(n-1.0D0))
    xc = xc / sig
    DO k=1,lag
        ac(k+1) = sum(xc(1:(n-k)) * xc((k+1):n)) / (n-1)
    END DO
END SUBROUTINE