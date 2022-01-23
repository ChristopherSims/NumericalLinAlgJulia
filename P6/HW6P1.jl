using LinearAlgebra, Random,SparseArrays,LinearOperators

include("./cg2.jl")

function lanczos(A,b,k)
    n = size(A,1)
    V = zeros(n,k+1)
    T = Tridiagonal(zeros(k-1), zeros(k), zeros(k-1))
    rho = 0.0
    V[:,1] = b/norm(b)
    for j=1:k
        y = A*V[:,j]
        for i=max(j-1,1):j
        T[i,j] = V[:,i]'*y
        y -= T[i,j]*V[:,i]
        end
        rho = norm(y)
        V[:,j+1] = y/rho
        if j < k
        T[j+1,j] = rho
        end
    end
    return V,T,rho
end
#### Part 1
#### MINRES
function Myminres(A, b,maxiter,tol)
    n = size(b, 1)
    x = zeros(n)
    nb = norm(b)
    res = zeros(maxiter + 1)
    res[1] = nb
    k=10
    V,T,alpha = lanczos(A,b,k)
    c, s = -1, 0
    delta = 0
    
    eps = 0
    pold, p = zeros(n), zeros(n)

    for i = 1:maxiter
        V,T,alpha = lanczos(A,b,k)
        # Solve min ||Tk' y - beta0 e1||
        tau  = c * delta + s * alpha
        gamma = s * delta - c * alpha

        delta = - c * beta
        epsold, eps = eps, s * beta

        # Build rotation
        c, s, eta = symortho(gamma, beta)

        mu = c * phi
        phi = s * phi

        pold, p = p, (qold - tau * p - epsold * pold) / eta

        # Update solution
        x += mu * p

        # Check residual
        res[i + 1] = abs(phi)
        print_residual(i, maxiter, res[i + 1], nb)
        if res[i + 1] / nb < tol
            return (x, res[1:i + 1])
        end
    end
    return (x, res)
end

## part 1 run
# test solution on random set
n = 10;
A = randn(n,n);
A = A*A' # make a SPD
b = ones(n)
#minres_sol = Myminres(A,b)

### Part 2 compare to CG
n = 100;
A = randn(n,n);
A = A'*A;
b = ones(n);
x_cg, hist_cg = cg2(A, b, 1000, 1e-8);



### simple cg
function simple_cg(A,b,maxiter,tol)
    rvec = zeros(maxiter)
    xvec = zeros(length(b),maxiter+1)
    k = 10;
    V,T,rho = lanczos(A,b,k)
    rhs = zeros(k)
    rhs[1] = norm(b)
    y = T \ rhs # solve the tridiagonal system
    xvec[:,1] = V[:,1:k]*y
    for i=2:maxiter
        V,T,rho = lanczos(A,xvec[:,i],k)
        y = T \ rhs # solve the tridiagonal system
        xvec[:,i+1] = V[:,1:k]*y
        rvec[i] = norm(x3 - A*b)/norm(b)
    end
    return x3, rvec
end

x3, rvec3 = simple_cg(A,b,100,1e-8)

plot(hist_cg)
plot(rvec3)








