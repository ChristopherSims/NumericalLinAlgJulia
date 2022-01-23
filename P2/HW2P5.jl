using LinearAlgebra
using Plots
#pyplot()

A = [3 -1.0; -1.0 4]
b = [-1, 2.0]
x0 = randn(2)
alpha = 0.1
x = copy(x0)
quadratic(A,b,x) = 0.5*x'*A*x - x'*b

x = y = range(-1, stop = 1, length = 25)
p = plot(x,y,(x,y) -> quadratic(A,b,[x,y]), camera=(150,30),st=:surface, xlim=(-1,1), ylim=(-1,1), zlim=(0,5))
plot!([x0[1]], [x0[2]], [quadratic(A,b,x0)])

x = copy(x0)
r = b - A*x
niter = 1000
n = size(A,1)
iter = 1
while norm(r) >= n*1e-5 && iter <= niter # set error to 1e5
  x .= x + alpha*r
  global iter += 1
  global r = b - A*x
  @show r
end
print("\n\n\n\n\n\n\n")

function Gdes(A::Matrix,b::Vector,niter::Integer,alpha::Float64,x0::Vector,res::Float64)
  x = copy(x0)
  r = b - A*x
  n = size(A,1)
  iter = 1
  while norm(r) >= n*res && iter <= niter # set error to 1e5
    x .= x + alpha*r
    iter += 1
    r = b - A*x
    @show r
  end
  return x
end
print("\n")

x = Gdes(A,b,niter,alpha,x0,1e-5)


