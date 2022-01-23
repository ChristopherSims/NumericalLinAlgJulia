using DelimitedFiles
using SparseArrays
using Plots
using LinearAlgebra
using BenchmarkTools
include("SparseMatrixCSR.jl")
coord = readdlm("chutes-and-ladders-coords.csv",',')
xc = coord[:,1]
yc = coord[:,2]

data = readdlm("chutes-and-ladders-matrix.csv",',')
TI = Int.(data[:,1])
TJ = Int.(data[:,2])
TV = data[:,3]
T = sparse(TI,TJ,TV, 101,101)
@time PR = SparseMatrixCSR(T) #### CSR is better for Lin alg
rowptr,colval,nzval = PR.rowptr, PR.colval, PR.nzval
#print(T[100,100])
start = 101
endstate = 100
k = 100
x0 = zeros(k+1) #create starting state
display(x0)

x0[start] = 1 # load stating state into vector
#AA = (I-T')\x0'
#show AA
function chutesnlads(x0::Vector, T::AbstractMatrix, k::Int) 
    X = zeros(length(x0),k+1) # fill in the code using the approximation evolution
    X[:,1] = x0 #starting state
    for i=1:k
        X[:,i+1] = (T)*X[:,i] # (I-T')X mod to solve
        #@show i
    end
    return X
end
@btime X = chutesnlads(x0,T,k)


function matvec(rowptr::Vector,colval::Vector,nzval::Vector,x::Vector)
y = zeros(maximum(colval))
#@show size(y,1)
#@show size(y,2)
for i=1:size(y,1)
    for k=rowptr[i]:rowptr[i+1]-1
            y[i] += nzval[k] * x[colval[k]]
    end
end
return y
end

alphar = 0.25
rel = 1e-5
niter = 100
function chutesnlads_rich(alphar::Float64,rel::Float64,rowptr::Vector,colval::Vector,nzval::Vector,niter::Integer,x0::Vector)
    #X = zeros(length(x0),k+1) # fill in the code using the approximation evolution
    x = x0 #starting state
    r = matvec(rowptr,colval,nzval,x)
    n = maximum(colval)
    iter = 1
    while norm(r) >= n*rel && iter <= niter # set error to 1e5
        #X[:,i+1] = T*X[:,i] # fix the initial probability x0
        r = matvec(rowptr,colval,nzval,x)
        x += alphar*r
        iter +=  1 
    end
    return x
end

@btime X3 = chutesnlads_rich(alphar,rel,rowptr,colval,nzval,niter,x0)
#display(X3)





#p = scatter(xc, yc, zcolor=sum(X[:,1:end-1],dims=2),
p = scatter(xc, yc, zcolor=X3,
    markersize=16, label="", marker=:square, markerstrokewidth=0, size=(400,400),
    xlims=(0.5,10.5),ylims=(0.5,10.5),aspect_ratio=:equal)
function draw_chutes(p)
    CL = [ 1 4 9 21 28 36 51 71 80 98 95 93 87 64 62 56 49 48 16
        38 14 31 42 84 44 67 91 100 78 75 73 24 60 19 53 11 26 6]
    for col=1:size(CL,2)
        i = CL[1,col]
        j = CL[2,col]
        if i > j # this is a chute
            plot!(p,[xc[i],xc[j]],[yc[i],yc[j]],color=2,label="")
        else
            plot!(p,[xc[i],xc[j]],[yc[i],yc[j]],color=1,label="")
        end
    end
    map(i->annotate!(p,xc[i],yc[i], text("$i", :white, 8)), 1:100)
    p
end
draw_chutes(p)
savefig("GameState.png")

