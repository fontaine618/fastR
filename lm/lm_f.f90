SUBROUTINE lm_f(n, p, X, Y, eps, max_iter, lr, beta)
    INTEGER :: n, p, i, max_iter
    DOUBLE PRECISION :: X(n, p), Y(n, 1)
    DOUBLE PRECISION :: beta(p, 1), g(p, 1)
    DOUBLE PRECISION :: eps, lr
    beta = 0.0D0
    DO i=1,max_iter
        g = matmul(transpose(X), Y - matmul(X, beta)) / n
        IF (sqrt(sum(g*g))*lr < eps) THEN
            exit
        END IF
        beta = beta + lr * g
    END DO
END SUBROUTINE