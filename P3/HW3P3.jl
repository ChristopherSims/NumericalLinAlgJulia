## A simple algorithm to solve a linear system
using LinearAlgebra
using Random
Random.seed!(10)
A = randn(4,4); A = A'*A; A = (A + A')/2
b = randn(4)
include("Poisson.jl")
##
function myreduce1(A)
  n = size(A,1)
  D = A[2:end,2:end]
  c = A[1,2:end]
  d = A[2:end,1]
  alpha = A[1,1]
  x = randn(n)
  @show x'*A*x
  @assert(x'*A*x > 0, "A is not positive definite ")
  L = Matrix(1.0I,n,n)
  L[2:end,1] = -d/alpha
  U = Matrix(1.0I,n,n)
  U[1,2:end] = -c/alpha 
  return L*A*U
end
myreduce1(A)

##
# This function saves the elimination structure as a set of matrices.
function myreduce_all(A::Matrix)
  A = copy(A) # save a copy
  n = size(A,1)
  L = Matrix(1.0I,n,n)
  U = Matrix(1.0I,n,n)
  d = zeros(n)
  x = randn(n)
  @show x'*A*x
  @assert(x'*A*x > 0, "A is not positive definite ")
  for i=1:n-1
    alpha = A[i,i]
    d[i] = alpha
    U[i,i+1:end] = A[i,i+1:end]/alpha
    L[i+1:end,i] = A[i+1:end,i]/alpha
    A[i+1:end,i+1:end] -= A[i+1:end,i]*A[i,i+1:end]'/alpha
  end
  d[n] = A[n,n]
  return L,U,d
end
A = randn(4,4); A = A'*A; A = (A + A')/2
#T = Tridiagonal(A)
#A = Matrix(T[:,:])
b = randn(4)
L,U,d = myreduce_all(A)
L*Diagonal(d)*U - A
#@show L
#@show U
#@show A


