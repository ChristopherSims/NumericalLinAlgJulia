using LinearAlgebra
include("Poisson.jl")
using Latexify
using Plots
#using Plot

n = 10
L = KLaplacian(10)

function Rich(A::Matrix,b::Vector,niter::Integer,alpha::Float64)
    n = size(A,1)
    @assert(n == size(A,2), "The matrix is not square")
    if niter <= 0
        niter = 10*size(A,1)
    end
    x = copy(b) # make a copy of the right hand side
    r = b - A*x # compute the residual
    iter = 0
    while norm(r) >= n*1e-5 && iter <= niter # set error to 1e5
        x = x + alpha*r
        iter += 1
        r = b - A*x
    end
    #@show r
    @show iter
    return x
end
function Rich_res(A::Matrix,b::Vector,niter::Integer,alpha::Float64)
    n = size(A,1)
    @assert(n == size(A,2), "The matrix is not square")
    if niter <= 0
        niter = 10*size(A,1)
    end
    x = copy(b) # make a copy of the right hand side
    r = b - A*x # compute the residual
    iter = 0
    while norm(r) >= n*1e-5 && iter <= niter # set error to 1e5
        x = x + alpha*r
        iter += 1
        r = b - A*x
    end
    return r
end

### Using test case from lecture
#A = Matrix(SymTridiagonal(ones(11), -0.5*ones(10)))
#A[1,2] = 0
#A[end,end-1] = 0


# Set alhpa = 1
A = -L #flip sign to converge
#b = randn(size(L,1))
b = ones(size(L,1))
x = Rich(A,b,1000,1.0)

U = A\b
display(x)
x_tex = latexify(x)
@show x_tex # print for latex
S = reshape(U,n,n)
surface(S)
savefig("laplaceM")

x_soln =zeros(11,size(L,1))
x_err =zeros(11,size(L,1))
a = 0:0.1:1
for ii=1:10
    U = A\b
    x_soln[ii,:] = Rich(A,b,1000,a[ii])
    x_err[ii,:] = U - x_soln[ii,:]
end

xnew = sum(x_err[:,:],dims = 2)

@show xnew
#plot(a,xnew)
#### For part 4
a = 0:0.01:0.3
x_soln =zeros(length(a),size(L,1))
x_err =zeros(length(a),size(L,1))
iter_vec = length(a)'
for ii=1:length(a)
    U = A\b
    x_soln[ii,:] = Rich(A,b,1000,a[ii])
    x_err[ii,:] = U - x_soln[ii,:]
end

xnew = sum(x_err[:,:],dims = 2)

plot(a,xnew)
#display(U)
#display(U-x)


#### Problem 1 part 3 Set N=20
n = 20
L = KLaplacian(n)
A = -L #flip sign to converge
#b = randn(size(L,1))
b = ones(size(L,1))
x = Rich(A,b,1000,0.3)
U = A\b
error = U-x
display(error) #0.3 does not converge
##############################
n = 20
L = KLaplacian(n)
A = -L #flip sign to converge
#b = randn(size(L,1))
b = ones(size(L,1))
x = Rich(A,b,1000,0.25) #0.25 converges better
U = A\b
error = U-x

#### For part 4
a = 0:0.01:0.3
x_soln =zeros(length(a),size(L,1))
x_err =zeros(length(a),size(L,1))
iter_vec = length(a)'
for ii=1:length(a)
    U = A\b
    x_soln[ii,:] = Rich(A,b,1000,a[ii])
    x_err[ii,:] = U - x_soln[ii,:]
end

xnew = sum(x_err[:,:],dims = 2)

#plot(a,xnew)

#### Part 5
a = Array((0:0.01:0.3))

res_soln =zeros(length(a),size(L,1))
for ii=1:length(a)
    res_soln[ii,:] = Rich_res(A,b,1000,a[ii])
end

resnew = sum(res_soln[:,:],dims = 2)
#plot(a,res_soln[:,5])
plot(a,resnew)
ylims!((0,100))
savefig("P120.png")
display(a)
display(res_soln[:,5])