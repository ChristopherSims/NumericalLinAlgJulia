using LinearAlgebra
using Random

Random.seed!(10)
n = 10;
A = randn(n,n);
A = A*A'; # make A SPD
B = randn(n)
D = Diagonal(B)
DI = inv(D)
AE = eigen(A);

C = inv(D)*A;

CE = eigen(C);

M1 = inv(D*A*D');
M1E = eigen(M1);
#display(AE.values)
display(CE.values)

F = A + D;
FE = eigen(F)
display(FE.values)
@show issymmetric(M1)
