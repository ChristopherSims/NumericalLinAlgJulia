using LinearAlgebra
using Random

n = 10;
A = randn(n);
#A = A*A'; # make A SPD
#D = Diagonal(A)

display(A)
b = randn(10);

function alg(b::Vector)
    nt = length(b)
    a = zeros(nt)
    for i = 1:n
        #a[i] = 
    end
end

