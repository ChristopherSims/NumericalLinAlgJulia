### KRYLOV conflicts with linalg so I am using another file
## Generate files for homework
using DelimitedFiles
using SparseArrays
data = readdlm("poisson2D-data.csv",',')
A = sparse(Int.(data[:,1]), Int.(data[:,2]),(data[:,3]))
A = (A + A')./2
b = vec(readdlm("poisson2D-rhs.csv"))
mytol = 1e-6
mymaxiter = 100
myrestart = 30
#A = Matrix(A)
## Part 3
using Krylov # ] clone https://github.com/JuliaSmoothOptimizers/Krylov.jl.git
using LinearOperators
P1 = opInverse(LowerTriangular(ichol(A)))
Mprecond = P1'*P1;
x, hist_pmin = Krylov.minres(A,b,M=Mprecond,itmax=mymaxiter, rtol=mytol,history=true);
pmin1 = hist_pmin.residuals;
plot(1:length(pmin1),pmin1,label = "pre-MINRES", lw = 3)
savefig("P5-3pmin.png") ## 88 iterations

Evec = eigvecs(Matrix(A))
Evals = eigvals(Matrix(A))
eigen_x = Evec[1,:]
eigen_y = Evec[2,:]
gr()
#scatter(eigen_x,eigen_y,label = "A")
scatter(Evals,label = "Avals")

Evec2 = real(eigvecs(Matrix(P1*A)))
Evals2 = real(eigvals(Matrix(P1*A)))
eigen_x2 = Evec2[1,:]
eigen_y2 = Evec2[2,:]
#scatter!(eigen_x2,eigen_y2,label ="Preconditioned")
scatter!(Evals2,label = "precond vals")
savefig("P5-3scatter.png") ## 88 iterations
