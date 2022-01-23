using LinearAlgebra
using Random

n = 10;
A = randn(n,n);
A = A*A' # make A SPD
F = factorize(A);
#display(F)
b = randn(n)
x = A\b
V = randn(n,n)
U = V
function update(x::Vector, b::Vector, F::Cholesky, U::Matrix, V::Matrix)
    nmax = size(U,2);
    y = 0
    for i = 1:nmax
        w = inv(A)*U[i,:];
        r = V[i,:]'*x./(1 .- V[i,:]'*w);
        y = x + r*w
        x = y
    end
    return y
end

abc = update(x,b,F,U,V)