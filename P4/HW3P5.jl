using LinearAlgebra
using Random
using TimerOutputs

#Random.seed!(10)
n = 10 # create n for inexing 
U = triu(randn(n,n)); # for prob 2
L = tril(randn(n,n));
b = randn(n); # for prob 2

## part 1 backsolve
function backsolve(A::Matrix,b::Vector)
    m,n = size(A) #index variables
    #@show n # for debug
    x = zeros(n) #create x variable
    for i = n:-1:1 #iterate down N-> 1
        x[i] = b[i] #diagonal condition
        for j = i+1:n #along the upperdiagonal i->N
            x[i] = x[i] -A[i,j]*x[j] #B - A[21]X[1] for every element
        end
        x[i] = x[i]/A[i,i] #compute (C/A[i,i])
    end
    return x
end
@time xx = backsolve(U,b)
#@show xx
#@show U
xtrue = U\b
#### PROBLEM 2
#### Compare backsolve to julia backsolve
@show norm(xtrue-xx) #1e-13


#### Part 1 forward solve
function ForwardSolve(L::Matrix,b::Vector)
    m,n = size(L)
    x = zeros(n) #create x variable
    # for i = 1:n
    #     x[i] = b[i]/A[i,i] # same as before
    #     for j=1:i-1
    #         x[i] = x[i]-A[i,j]*x[j]
    #     end
    #     x[i] = x[i]/A[i,i] #compute (C/A[i,i])
    # end
    for i = 1:n
        tmp = b[i]
        for j = 1:i-1
            tmp -= L[i,j]*x[j]
        end
        x[i] = tmp/L[i,i]
    end
    return x
end
yy = ForwardSolve(L,xtrue)
xtruey = L\xtrue
#### PROBLEM 2
#### Compare backsolve to julia backsolve
@show norm(xtruey-yy) #1e-13
#@show yy



### Part 3 compare linear solvers
A = randn(n,n) # create random matrix 
A = A'*A; A = (A + A')/2
b = randn(n); # for prob 2
x = randn(n)
#@show x'*A'x

function mylinsolv(A::Matrix,b::Vector)
    L,U,p = lu(A)
    #@show p
    y = ForwardSolve(L,b)
    Z = backsolve(U,y)
    xt = Z[p]
    #display(Z)
    #display(Zt)
    return Z
end
@time xmy = mylinsolv(A,b)
### Compare

@time xtrue2 = A\b
#display(xtrue2)
@show norm(xtrue2-xmy) 

L,U,p = lu(A)
y2 = L\b
z2 = U\y2
#z2=z2[p]
@show norm(xtrue2-z2) 

