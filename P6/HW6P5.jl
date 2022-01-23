using Random, SparseArrays, IterativeSolvers, Plots, LinearAlgebra

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
sigma = 1.7e-2
A = (A-sigma*I)

## Part 1
x, hist_min = minres(A,b,abstol = mytol, maxiter = mymaxiter, log=true);
x, hist_gm = gmres(A,b,restart = myrestart, abstol = mytol, maxiter = mymaxiter, log=true);
#hist_min.data[:resnorm]
min1 = hist_min.data[:resnorm];
gm1 = hist_gm.data[:resnorm];
plot(1:length(min1),min1,label = "MINRES", lw = 3)
savefig("P5-1minres.png")
plot(1:length(gm1),gm1,label = "GMRES", lw = 3)
savefig("P5-1gmres.png")

plot(1:length(min1),min1,label = "MINRES", lw = 3)
plot!(1:length(gm1),gm1,label = "GMRES", lw = 3)
savefig("P5-1both.png")


## Part 2
include("precond.jl")
#Lchol = ichol(A)
M1 = LowerPrecond(ichol(A))
x, hist_pgm = gmres(A,b,restart = myrestart, abstol = mytol, maxiter = mymaxiter,Pl=M1, log=true)

pgm2 = hist_pgm.data[:resnorm];
plot(1:length(pgm2),pgm2,label = "pre-GMRES", lw = 3)
savefig("P5-2pgm.png") ## 88 iterations

L,U = ilu(A)
M2 = LUPrecond(L,U)

x, hist_pgm2 = gmres(A,b,restart = myrestart, abstol = mytol, maxiter = mymaxiter,Pl=M2, log=true)
pgm2lu = hist_pgm2.data[:resnorm];
plot(1:length(pgm2lu),pgm2lu,label = "pre-GMRES LU", lw = 3)
savefig("P5-2pgmLU.png") ## 88 iterations

