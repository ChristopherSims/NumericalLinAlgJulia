using DelimitedFiles
using SparseArrays
using Plots, IterativeSolvers,LinearAlgebra
#using LinearAlgebra.minres 


######################
#prreprosecing
################
coord = readdlm("chutes-and-ladders-coords.csv",',')
xc = coord[:,1]
yc = coord[:,2]

data = readdlm("chutes-and-ladders-matrix.csv",',')
TI = Int.(data[:,1])
TJ = Int.(data[:,2])
TV = data[:,3]
T = sparse(TI,TJ, TV, 101,101)
#print(T[100,100])
start = 101
endstate = 100

k = 100
x0 = zeros(k+1) #create starting state
x0[start] = 1 # load stating state into vector
b = zeros(k+1)
b[endstate] = 1
##########################
# END preprocesing
#####################


##########################
# Part 2
###########################
### MINRES
nmax = 100
X = zeros(k+1,nmax)
X[:,1] = x0
#xnew = minres!(X[:,1], T'*T, b)
rvec = zeros(nmax)
for i = 1:nmax-1
    rvec[i] = norm(T'*T*X[:,i] - T'b)/norm(b) # relative 2 norm risd
    X[:,i+1] = minres!(X[:,i], T'*T, T'*b)
end

# MINRES algorithms have M matrix-vector products
#plot(1:nmax,rvec)
#xlims!(0,10)
#ylims!(0,1e-7)

##########################
# CG
##########################
nmax = 100
X = zeros(k+1,nmax)
X[:,1] = x0
CG_rvec = zeros(nmax)
for i = 1:nmax-1
    CG_rvec[i] = norm(T'*T*X[:,i] - T'b)/norm(b) # relative 2 norm risd
    X[:,i+1] = cg!(X[:,i], T'*T, T'*b)
end

plot(1:nmax,rvec,label = "MINRES", lw = 3)
plot!(1:nmax,CG_rvec, label = "CG", lw = 3)
#xlims!(0,10)
ylims!(0,1e-7)
savefig("P4cgmin.png")

# for i = 1:100
#     r1 = rvec[i]
#     r2 = CG_rvec[i]
#     print("$i $r1 $r2 \n")
# end



