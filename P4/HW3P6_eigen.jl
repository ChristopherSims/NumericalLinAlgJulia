##
## A simple algorithm to solve a linear system
using LinearAlgebra, Random
function solve1_pivot2(A::Matrix, b::Vector)
  m,n = size(A)
  @assert(m==n, "the system is not square")
  @assert(n==length(b), "vector b has the wrong length")
  if n==1
    return [b[1]/A[1]]
  else
    # let's make sure we have an equation
    # that we can eliminate!
    # let's try that again, where we pick the
    # largest magnitude entry!
    maxval = abs(A[1,1])
    newrow = 1
    for j=2:n
      if abs(A[j,1]) > maxval
        newrow = j
        maxval = abs(A[j,1])
      end
    end
    if maxval < eps(1.0)
      error("the system is singular")
    end
    #@show newrow
    # swap rows 1, and newrow
    if newrow != 1
      tmp = A[1,:]
      A[1,:] .= A[newrow,:]
      A[newrow,:] .= tmp
      b[1], b[newrow] = b[newrow], b[1]
    end
    D = A[2:end,2:end]
    c = A[1,2:end]
    d = A[2:end,1]
    alpha = A[1,1]
    y = solve1_pivot2(D-d*c'/alpha, b[2:end]-b[1]/alpha*d)
    gamma = (b[1] - c'*y)/alpha
    return pushfirst!(y,gamma)
  end
end
#A = [1e-18 1.0; 1.0 1.0]
n = 5
#Random.seed!(10)
A = randn(n,n)
A = A'*A
b = zeros(n)
e = eigvals(A)
x = e[1]*randn(n)
Anew = (A-e[1]*I)
display(Anew)
display(e)
#b = [5; 6.0]

#xtrue = A\b
#bnew = e[1]*xtrue
xmy = solve1_pivot2(Anew,b)
xtruep2 = A\b
@show norm(xtruep2-xmy)


function solve_eigen(A::Matrix, e::Float64)
  A = (A-e*I)
  m,n = size(A)
  b=zeros(n)
  b[end] = 1
  if n==1
    return [b[1]/A[1]]
  else
    # let's make sure we have an equation
    # that we can eliminate!
    # let's try that again, where we pick the
    # largest magnitude entry!
    maxval = abs(A[1,1])
    newrow = 1
    for j=2:n
      if abs(A[j,1]) > maxval
        newrow = j
        maxval = abs(A[j,1])
      end
    end
    if maxval < eps(1.0)
      error("the system is singular")
    end
    #@show newrow
    # swap rows 1, and newrow
    if newrow != 1
      tmp = A[1,:]
      A[1,:] .= A[newrow,:]
      A[newrow,:] .= tmp
      b[1], b[newrow] = b[newrow], b[1]
    end
    D = A[2:end,2:end]
    c = A[1,2:end]
    d = A[2:end,1]
    alpha = A[1,1]
    y = solve1_pivot2(D-d*c'/alpha, b[2:end]-b[1]/alpha*d)
    gamma = (b[1] - c'*y)/alpha
    return pushfirst!(y,gamma)
  end
end

x2 = solve_eigen(A,e[1])
xrev = x2/(minimum(x2))
display(x2/(minimum(x2)))
V = eigvecs(A)
display(V)
xrev./V[:,1]

