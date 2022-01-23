using Plots, IterativeSolvers,LinearAlgebra  
using ZipFile, DelimitedFiles,SparseArrays
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
mytol = 1e-6
mymaxiter = 100

x, hist_min = minres(T'*T,T'*b,abstol = mytol, maxiter = mymaxiter, log=true);
x = reverse(sort(x))
display(x[1:5])
##lanczos
function load_data()
    r = ZipFile.Reader("wikipedia-2005.zip")
    try
        @assert length(r.files) == 1
        f = first(r.files)
        data = readdlm(f,'\n',Int)
        n = data[1]
        colptr = data[2:2+n] # colptr has length n+1
        rowval = data[3+n:end] # n+2 elements before start of rowval
        A = SparseMatrixCSC(n,n,colptr,rowval,ones(length(rowval))) |>
            A->begin ei,ej,ev = findnz(A); d = sum(A;dims=2);
        return sparse(ej,ei,ev./d[ei], size(A)...) end
    finally
        close(r)
    end
end
P = load_data()
nmax = size(P,1)
v = ones(nmax)/nmax;
x2, hist_min2 = minres(P,v,abstol = mytol, maxiter = mymaxiter, log=true);
x2 = reverse(sort(x2))
display(x2[1:5])


