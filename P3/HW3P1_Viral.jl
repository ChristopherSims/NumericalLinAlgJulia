## Generate a simple spatial graph model to look at.
using NearestNeighbors, Distributions, SparseArrays
using Latexify, LinearAlgebra
using Statistics
 
function spatial_graph_edges(n::Integer,d::Integer;degreedist=LogNormal(log(4),1))
xy = rand(d,n)
T = BallTree(xy)
# form the edges for sparse
ei = Int[]
ej = Int[]
for i=1:n
    deg = min(ceil(Int,rand(degreedist)),n-1)
    idxs, dists = knn(T, xy[:,i], deg+1)
    for j in idxs
        if i != j
            push!(ei,i)
            push!(ej,j)
            end
        end
    end
    return xy, ei, ej
end
function spatial_network(n::Integer, d::Integer; degreedist=LogNormal(log(3),1))
    xy, ei, ej = spatial_graph_edges(n, d;degreedist=degreedist)
    A = sparse(ei,ej,1,n,n)
    return max.(A,A'), xy
end
using Random
#Random.seed!(10) # ensure repeatable results...
#A,xy = spatial_network(10, 2)
""" return a new "social" graph where we have implemented
social distancing by removing an f-fraction of your neighbors
based on spatial proximity. So f=0 is the original network
and f=1 is the empty network."""
function social_distance(A::SparseMatrixCSC, xy::Matrix, f::Real)
# access the CSC arrays directly, see help on nzrange
    rowval = rowvals(A)
    n = size(A,1)
    new_ei = Vector{Int}()
    new_ej = Vector{Int}()
    for j = 1:n
        neighbors = Vector{Int}()
        dists = Vector{Float64}()
        myxy = @view xy[:,j] # don't make a copy
        for nzi in nzrange(A, j)
        # edge from (i,j)
        i = rowval[nzi]
        push!(neighbors, i)
        push!(dists, norm(xy[:,i]-myxy))
        end
        p = sortperm(dists) # sort distances
        nkeep = ceil(Int, (1-f)*length(dists))
        for i=1:nkeep
            push!(new_ei, neighbors[p[i]])
            push!(new_ej, j)
        end
    end
    A = sparse(new_ei,new_ej,1, size(A,1),size(A,2))
    return max.(A,A')
end

#A = social_distance(A,xy,0.6)
# #print(A)
# #print("\n")
# #print(xy)
# LatA = latexify(A)
# #@show LatA
# Latxy = latexify(xy')
# #@show Latxy
#### 
# Plot Results (6.2)
####
using Plots
function plotgraph(A::SparseMatrixCSC,xy::AbstractArray{T,2};kwargs...) where T
    px,py = zeros(T,0),zeros(T,0)
    P = [px,py]
    rows = rowvals(A)
    skip = NaN.*xy[:,begin] # first row
    for j=1:size(A,2) # for each column
        for nzi in nzrange(A, j)
            i = rows[nzi]
            if i > j
                push!.(P, @view xy[:,i])
                push!.(P, @view xy[:,j])
                push!.(P, skip)
            end
        end
    end
    plot(px,py;framestyle=:none,legend=false,kwargs...)
end

# plotgraph(A,xy,alpha=0.25)
# scatter!(xy[1,:],xy[2,:],
#     markersize=2, markerstrokewidth=0, color=1)
# savefig("Viral2.png")

##### Viral Pt 3
function evolve(x::Vector, p::Real, A::AbstractMatrix)
    log_not_infected = log.(1 .- p.*x)
    y = 1 .- exp.(A*log_not_infected)
    y = max.(y, x)
    end

    """
    Run k steps of the evolution and return the results as a matrix.
    Each column of the matrix has the probabilities that the node
    is infected under the `wrong` evolve function.
    The first column of X is the initial vector x0.
    At each iteration, we make sure the probabilities are at least x0 and these
    are fixed.
    """
function evolve_steps(x0::Vector, p::Real, A::AbstractMatrix, k::Int)
    X = zeros(length(x0),k+1)
    X[:,1] = x0
    for i=1:k
        X[:,i+1] = max.(evolve(X[:,i], p, A), X[:,1]) # fix the initial probability x0
    end
    return X
end
# p = 1
# x0 = zeros(size(A,1))
# x0[1] = 1
# k = 10
# XF = evolve_steps(x0, p, A, k)
# #@show XF
# LatXF = latexify(XF)
# #@show LatXF

# C = sum(XF[2:end,:],dims=1)
# #@show C
# plot(C',legend=false)
# savefig("Viral3.png")

""" Run k steps of the approximate evolution and return the results
as a matrix. Each column of the matrix has the probabilities that
the node is infected under the wrong evolve function. The first column of X is the initial vector x0. At each iteration, we make sure
the probabilities are at least x0 and these are fixed. """ 
function approx_evolve(x::Vector, p::Real, A::AbstractMatrix)
    y = p.*(A*x)
    end    
function approx_evolve_steps(x0::Vector, p::Real, A::AbstractMatrix, k::Int) 
    X = zeros(length(x0),k+1) # fill in the code using the approximation evolution
    X[:,1] = x0
    for i=1:k
        X[:,i+1] = max.(approx_evolve(X[:,i], p, A), X[:,1]) # fix the initial probability x0
    end
    return X
end
# p = 0.10
# x0 = zeros(size(A,1))
# x0[1] = 1
#Y = approx_evolve_steps(x0, p, A, 15)
#X = evolve_steps(x0, p, A, 15)
#plot(sum(X[2:end,:],dims=1)',legend=false,color ="red")
#savefig("ViralEvolve.png")
#plot(sum(Y[2:end,:],dims=1)',legend=false,color="Blue")
#savefig("Viralapprox.png")
#################################
# Code for HW3 starts here
################################
Random.seed!(10) # ensure repeatable results...
A,xy = spatial_network(10, 2)
p = 0.2
x0 = zeros(size(A,1))
x0[1] = 1
Y = approx_evolve_steps(x0, p, A, 10)
#display(Y[:,end])
RQT = Y[:,end]'*p*A*Y[:,end]
#display(RQT)
RQB = Y[:,end]'*Y[:,end]
#display(RQB)
RQ = RQT/RQB
print("10 steps\n")
display(RQ)

###### 50 steps
Y = approx_evolve_steps(x0, p, A, 50)
#display(Y[:,end])
RQT = Y[:,end]'*p*A*Y[:,end]
#display(RQT)
RQB = Y[:,end]'*Y[:,end]
#display(RQB)
RQ = RQT/RQB
print("50 steps\n")
display(RQ)

###### 100 steps
Y = approx_evolve_steps(x0, p, A, 100)
#display(Y[:,end])
RQT = Y[:,end]'*p*A*Y[:,end]
#display(RQT)
RQB = Y[:,end]'*Y[:,end]
#display(RQB)
RQ = RQT/RQB
print("100 steps\n")
display(RQ)


#######################
# part 1 n=1000
####################
#####1000x1000 network
Random.seed!(10) # ensure repeatable results...
A,xy = spatial_network(1000, 2)
x0 = zeros(size(A,1))
x0[end] = 1
Y = approx_evolve_steps(x0, p, A, 10)
#display(Y[:,end])
RQT = Y[:,end]'*p*A*Y[:,end]
#display(RQT)
RQB = Y[:,end]'*Y[:,end]
#display(RQB)
RQ = RQT/RQB
print("n=1000 10 steps\n")
display(RQ)

###### 50 steps
Y = approx_evolve_steps(x0, p, A, 50)
#display(Y[:,end])
RQT = Y[:,end]'*p*A*Y[:,end]
#display(RQT)
RQB = Y[:,end]'*Y[:,end]
#display(RQB)
RQ = RQT/RQB
print("n=1000 50 steps\n")
display(RQ)

###### 100 steps
Y = approx_evolve_steps(x0, p, A, 100)
#display(Y[:,end])
RQT = Y[:,end]'*p*A*Y[:,end]
#display(RQT)
RQB = Y[:,end]'*Y[:,end]
#display(RQB)
RQ = RQT/RQB
print("n=1000 100 steps\n")
display(RQ)

####################
# part 2 power method
###################
function power_steps(x0::Vector, p::Real, A::AbstractMatrix, k::Int) 
    X = zeros(length(x0),k+1) # fill in the code using the approximation evolution
    X[:,1] = x0
    for i=1:k
        X[:,i+1] = A*X[:,i]
    end
    return X
end
################
# power 10
###############
A,xy = spatial_network(10, 2)
p = 0.2
x0 = zeros(size(A,1))
#x0[end] = 1
#### 10x10 network
x0 = ones(size(A,1))
Y = power_steps(x0, p, A, 100)
X_p = Y[:,end]./(minimum(Y[:,end]))
#@show Y[:,end]
#@show minimum(Y[:,end])
#display(X_p)
RQT = X_p'*p*A*X_p
#display(RQT)
RQB = X_p'*X_p
#display(RQB)
RQ = RQT/RQB
print("part 2 power method:\n")
display(RQ)

#####1000x1000 network
Random.seed!(10) # ensure repeatable results...
A,xy = spatial_network(1000, 2)




##### Part 3 
### Social Distancing
Random.seed!(10) # ensure repeatable results...
A,xy = spatial_network(1000, 2)
A = social_distance(A,xy,0.85)
p=0.2
x0 = randn(size(A,1))
Y = power_steps(x0, p, A, 100)
X_p = Y[:,end]/(minimum(Y[:,end]))
#display(X_p)
RQT = X_p'*p*A*X_p
#display(RQT)
RQB = X_p'*X_p
#display(RQB)
RQ = RQT/RQB
print("SD 10x10 100 steps:\n")
display(RQ)
############# Set p lower
p=0.065
Random.seed!(10) # ensure repeatable results...
A,xy = spatial_network(1000, 2)
#A = social_distance(A,xy,0.85)
x0 = randn(size(A,1))
Y = power_steps(x0, p, A, 100)
X_p = Y[:,end]/(minimum(Y[:,end]))
#display(X_p)
RQT = X_p'*p*A*X_p
#display(RQT)
RQB = X_p'*X_p
#display(RQB)
RQ = RQT/RQB
print("set P lower, 100 steps\n")
display(RQ)



####### part 4
## Vacinne
p=0.2
#display(B)
Mvec = rand(1:2000,25)
RQ_vec = zeros(25)
for jj = 1:25
    Random.seed!(10) # ensure repeatable results...
    A,xy = spatial_network(1000, 2)
    A = social_distance(A,xy,0.1)
    B = A[A .!= 0] ### Find non zero elements
    maxn = Int(round(length(B.nzval)*(0.72))) #max number of non zeros * ones we want to be zero
    rng = MersenneTwister(Mvec[jj])
    idx = shuffle!(rng,B.nzind) #randomly permute indicies
    idx = idx[1:maxn]
    xynz = zeros(length(B.nzval),2)
    for i = 1:length(idx)
        A.nzval[idx[i]] = 0
    end
    x0 = randn(size(A,1))
    Y = power_steps(x0, p, A, 100)
    X_p = Y[:,end]/(minimum(Y[:,end]))
    #display(X_p)
    RQT = X_p'*p*A*X_p
    #display(RQT)
    RQB = X_p'*X_p
    #display(RQB)
    RQ = RQT/RQB
    #display(RQ)
    RQ_vec[jj] = RQ
end
RQ_mean = mean(RQ_vec)
print("vaccine mean\n")
display(RQ_mean)
RQ_std = std(RQ_vec)
print("vaccine std\n")
display(RQ_std)