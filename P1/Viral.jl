## Generate a simple spatial graph model to look at.
using NearestNeighbors, Distributions, SparseArrays
using Latexify
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
Random.seed!(10) # ensure repeatable results...
A,xy = spatial_network(10, 2)
#print(A)
#print("\n")
#print(xy)
LatA = latexify(A)
#@show LatA
Latxy = latexify(xy')
#@show Latxy
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

plotgraph(A,xy,alpha=0.25)
scatter!(xy[1,:],xy[2,:],
    markersize=2, markerstrokewidth=0, color=1)
savefig("Viral2.png")

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
p = 1
x0 = zeros(size(A,1))
x0[1] = 1
k = 10
XF = evolve_steps(x0, p, A, k)
#@show XF
LatXF = latexify(XF)
#@show LatXF

C = sum(XF[2:end,:],dims=1)
#@show C
plot(C',legend=false)
savefig("Viral3.png")

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
p = 0.10
x0 = zeros(size(A,1))
x0[1] = 1
Y = approx_evolve_steps(x0, p, A, 15)
X = evolve_steps(x0, p, A, 15)
plot(sum(X[2:end,:],dims=1)',legend=false,color ="red")
savefig("ViralEvolve.png")
plot(sum(Y[2:end,:],dims=1)',legend=false,color="Blue")
savefig("Viralapprox.png")


    
