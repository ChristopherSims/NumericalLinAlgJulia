using Random, SparseArrays, IterativeSolvers, LinearAlgebra, Plots

using BenchmarkTools
include("multigrid_functions.jl")

n = [31, 63, 127, 255, 511, 1023]
nx = n
ny = n
bruntime = zeros(length(n))

function runpoi(nx,ny)
    x = solve_poisson_direct(poisson_setup(nx,ny, (x,y) -> 1))
    return x
end
for i = 1:length(n)
     bruntime[i] = @elapsed runpoi(nx[i],ny[i])
     #@time runpoi(nx[i],ny[i])
end
# using Plots
# pyplot()
# #surface(X)

# plot(n,bruntime,xaxis=:log, yaxis=:log)
# savefig("P6-1.png")
# for i = 1:length(bruntime)-1
#     print(bruntime[i+1]/bruntime[i] )
#     print("\n")
# end
######################
### part 2 and 3
#####################
# maxiter = 10000
# R = poisson_setup(nx[1],ny[1], (x,y) -> 1);
# b = randn(size(R,1),size(R,2))
# b = R
# U = zeros(size(R,1),size(R,2),maxiter)
# #U[:,:,1] = R
# ek1 = zeros(maxiter)
# ek = zeros(maxiter)
# tot_error = zeros(maxiter)
# for i = 1:maxiter-1
#     U[:,:,i+1]  = apply_poisson_jacobi(U[:,:,i],b);
#     ek1[i] = norm(U[:,:,i+1] - R)/norm(R)
#     ek[i] = norm(U[:,:,i] - R)/norm(R)
#     tot_error[i] = ek1[i]/ek[i]
# end
# plot(10:maxiter-4,tot_error[10:maxiter-4].-1)
# ylims!(0,1e-6)
# savefig("P6-3.png")


######################
### part 6
#####################
x = solve_poisson_direct(poisson_setup(31,31, (x,y) -> 1))
xtrue = solve_poisson_direct(poisson_setup(63,63, (x,y) -> 1))
x2 = interpolate(x)
surface(x2-xtrue)
enorm = norm(x2-xtrue)
######################
### part 7
#####################
maxiter = 10000
R = poisson_setup(nx[2],ny[2], (x,y) -> 1);
#b = randn(size(R,1),size(R,2))
b = R
U = zeros(size(R,1),size(R,2),maxiter)
U[:,:,1] = x2
ek1 = zeros(maxiter)
ek = zeros(maxiter)
tot_error = zeros(maxiter)
for i = 1:maxiter-1
    U[:,:,i+1]  = apply_poisson_jacobi(U[:,:,i],b);
    ek1[i] = norm(U[:,:,i+1] - R)/norm(R)
    ek[i] = norm(U[:,:,i] - R)/norm(R)
    tot_error[i] = abs(ek1[i]/ek[i])
end
#display(tot_error)
plot(10:maxiter-4,tot_error[10:maxiter-4].-1, lw = 3)
ylims!(-1e-8,1e-8)
savefig("P6-7.png")
maxiter = 10000
R = poisson_setup(nx[2],ny[2], (x,y) -> 1);
b = randn(size(R,1),size(R,2))
b = R
U = zeros(size(R,1),size(R,2),maxiter)
#U[:,:,1] = R
ek1 = zeros(maxiter)
ek = zeros(maxiter)
tot_error = zeros(maxiter)
for i = 1:maxiter-1
    U[:,:,i+1]  = apply_poisson_jacobi(U[:,:,i],b);
    ek1[i] = norm(U[:,:,i+1] - R)/norm(R)
    ek[i] = norm(U[:,:,i] - R)/norm(R)
    tot_error[i] = ek1[i]/ek[i]
end
plot(10:maxiter-4,tot_error[10:maxiter-4].-1)
ylims!(0,1e-5)
savefig("P6-7ref.png")
######################
### part 10
#####################
n2 = [31, 63, 127, 255, 511, 1023]
nx = n2
ny = n2
niter = 25
bruntime2 = zeros(length(n2))
for i = 1:length(n2)
     bruntime2[i] = @elapsed poisson_multigrid_residual(nx[i],ny[i],niter);
     #@time runpoi(nx[i],ny[i])
end
plot(n2,bruntime2,label ="multigrid", xaxis=:log, yaxis=:log)
plot!(n,bruntime,label = "poissen-direct", xaxis=:log, yaxis=:log)
ylabel!("time seconds")
ylabel!("matrix size")
savefig("P6-10.png")



