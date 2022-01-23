function cg2(A, b, max_it, tol)
# conjugate gradient in Julia
# original code from a matlab version written by
#   -- Iterative template routine --
#      Univ. of Tennessee and Oak Ridge National Laboratory
#      October 1, 1993
#      Details of this algorithm are described in "Templates for the
#      Solution of Linear Systems: Building Blocks for Iterative
#      Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
#      Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
#      1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
#
#   x, res = cg(A, b, max_it, tol)
#
#  cg.m solves the symmetric positive definite linear system Ax=b
#  using the Conjugate Gradient method with preconditioning.
#
#  input   A        REAL symmetric positive definite matrix
#          b        REAL right hand side vector
#          max_it   INTEGER maximum number of iterations
#          tol      REAL error tolerance
#
#  output  x        REAL solution vector
#          error    REAL residual norms from each iteration
#
# Code translated to Julia and modified by Huda Nassar for CS 515 - Fall 2017

bnrm2 = norm(b)
if bnrm2 == 0.0
    bnrm2 = 1.0
end

x = copy(b)
x[:] .= 0
r = copy(b)
error = norm(r) / bnrm2
rho_1 = 0
p = 0
res = zeros(max_it)

if error < tol
    return x,res
end

lastiter = 1
for iter = 1:max_it
    z  = r
    rho = dot(r,z)

    if iter > 1
        beta = rho / rho_1
        p = z + beta*p
    else
        p = z
    end

    q = A*p
    alpha = rho / (p'*q)
    x = x + alpha * p                    # update approximation vector

    r = r - alpha*q                      # compute residual
    res[iter] = norm(r) / bnrm2          # check convergence
    if res[iter] <= tol
        break
    end
    rho_1 = rho
	lastiter = iter
end

res = res[1:lastiter]

return x,res
end
