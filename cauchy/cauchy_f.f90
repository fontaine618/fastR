SUBROUTINE dcauchy(n, x, loc, scale, d)
    INTEGER :: n 
    DOUBLE PRECISION :: x(n), loc, scale
    DOUBLE PRECISION :: d(n)
    DOUBLE PRECISION :: PI = 4.0D0 * DATAN(1.0D0)
    d = x-loc
    d = scale / (PI * (d**2 + scale**2))
END SUBROUTINE