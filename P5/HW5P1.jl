using LinearAlgebra
using Random

#### Pre proccess input for ploblem before problem
A = randn(10,10)
A = A*A' # make A SPD
C = cholesky(A)
L = C.L
U = C.U
u = randn(10)


## Rank 1 cholesky update
function cholupdate(L,x)
    n = length(x)
    for k = 1:n
        r = sqrt(L[k,k]^2 + x[k]^2)
        c = r/L[k,k];
        s = x[k]/L[k,k];
        L[k,k] = r;
        if k < n
            L[(k+1):n, k] = (L[(k+1):n, k]+ s * x[(k+1):n]) / c;
            x[(k+1):n] = c * x[(k+1):n] - s * L[(k+1):n, k];
        end
    end
    return L
end


Lnew = cholupdate(L,u')

display(Lnew) ## Show update cholesky



