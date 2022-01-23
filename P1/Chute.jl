using DelimitedFiles
using SparseArrays
using Plots
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


function chutesnlads(x0::Vector, T::AbstractMatrix, k::Int) 
    X = zeros(length(x0),k+1) # fill in the code using the approximation evolution
    X[:,1] = x0 #starting state
    for i=1:k
        X[:,i+1] = T*X[:,i] # fix the initial probability x0
    end
    return X
end
k = 100
x0 = zeros(k+1) #create starting state
x0[start] = 1 # load stating state into vector
#@show x0 
X = chutesnlads(x0,T,k)
#x = T*y0
#plot(sum(X[2:end,:],dims=1)',legend=false,color ="red")
plot(X[1,:],legend=false)
for i =1:k+1
    plot!(X[i,:],legend=false)
end
ylims!((0,0.2))
title!("All end state prob")
savefig("allstart.png")

plot(X[100,:],legend=false,color ="red")
title!("100 state prob")
savefig("lads.png")

P = zeros(size(X,1))
for i=1:k
    P[i+1] += i*X[100,i] + P[i]
end
evec = P[2:end] - P[1:end-1]
for ii = 1:length(evec)
    if evec[ii] < 0.15 && ii > 10
        print("\n")
        print(evec[ii])
        print(" ")
        print(P[ii])
    end 
end
#@show evec
#@show P
#@show size(X)
plot(P,legend=false,color="black")
savefig("total_prob.png")
p = scatter(xc, yc, zcolor=sum(X[:,1:end-1],dims=2),
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