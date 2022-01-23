using LinearAlgebra
using SparseArrays
include("SparseMatrixCSR.jl")
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
function csc_lookup(colptr, rowval, nzval, m, n, i, j)
    B=0
    if (i <= n) || (j <= m)
        for jj=1:length(colptr)-1 # for each column
            for nzi=colptr[jj]:colptr[jj+1]-1 # for each entry in the column
                ii = rowval[nzi]
                v = nzval[nzi]
                if (ii == i) && (jj == j)
                    B = v
                end
            end
        end
        
    end
    return B
end
m = 6
n = 6
M = rand(m,n)*5
S = sparse(M)
niter = 100
PR = SparseMatrixCSR(S) #### CSR is better for Lin alg
rowptr,colval,nzval = PR.rowptr, PR.colval, PR.nzval
sb = ones(m)
colptr = S.colptr
rowval = S.rowval
nzval = S.nzval
V = range(1,step=2,stop=n-1) # Create a vector 1:2:n
#V = range(1,step=1,stop=n-1)
x0 = randn(n)
x = copy(x0)
display(x)
############ Iterate over n
for i = V
    global r = matvec(rowptr,colval,nzval,x) - sb
    global iter = 1
    while norm(r[i]) >= 1e-5 && iter <= niter # set error to 1e5
        global r = matvec(rowptr,colval,nzval,x) - sb
        global temp = csc_lookup(colptr,rowval,nzval,m,n,i,i)
        global temp2 = csc_lookup(colptr,rowval,nzval,m,n,i+1,i+1) 
        c1 = -r[i]/temp
        c2 = -r[i+1]/temp2
        C = [c1 c2]'
        #@show (norm(r[i]) + norm(r[i+1]))
        #@show C
        global ei = zeros(length(x))
        global ei[i] = 1
        global ej = zeros(length(x))
        global ej[i+1] = 1
        global E_ij = hcat(ei,ej)
        #display(E_ij)
        #display(E_ij*C)
        global x += vec(E_ij*C)
        #display(x)
        global iter += 1
    end
    @show iter
    #@show i
    #@show x
end
@show x