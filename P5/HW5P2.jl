using LinearAlgebra
using Random

Random.seed!(10)
n = 10;
B = randn(n);
A = Diagonal(B); # make A SPD
M2 = randn(n,n)
M2 = M2*M2'

T = cos(A)^2 + sin(A)^2;

display(T)


M2T = cos(M2)^2 + sin(M2)^2;

display(M2T)

display(diag(M2T))

C = cos(M2)
S = sin(M2)

C2 = C*C
S2 = S*S

