using LinearAlgebra
using BenchmarkTools

########## Load Test File
using ZipFile, DelimitedFiles,SparseArrays
include("SparseMatrixCSR.jl")
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
#P = I-P
colptr,rowval,nzval = P.colptr, P.rowval, P.nzval
PR = SparseMatrixCSR(P) #### CSR is better for Lin alg
rowptr,colval,nzval = PR.rowptr, PR.colval, PR.nzval ## cast vars to vectors
############ End load test File

########### load secondary test file
# using SparseArrays, Random
# Random.seed!(0)
# P = sprand(1_000_000,1_000_000, 15/1_000_000) |>
#     A->begin fill!(A.nzval,1); A; end |>
#     A->begin ei,ej,ev = findnz(A); d = sum(A;dims=2);
#         return sparse(ej,ei,ev./d[ei], size(A)...); end
# P = I-P
# @time colptr,rowval,nzval = P.colptr, P.rowval, P.nzval
# @time PR = SparseMatrixCSR(P) #### CSR is better for Lin alg
# rowptr,colval,nzval = PR.rowptr, PR.colval, PR.nzval ## cast vars to vectors
########## Testing array #2


#### csc multiply
function cscmatvec(TI,TJ,TV,x,m)
    y = zeros(m)
    for nzi in 1:length(TV)
        i = TI[nzi]
        j = TJ[nzi]
        v = TV[nzi]
        if i == j
            y[i] += 1 - v*x[i]
            @show i
        else
            y[i] += v*x[i]
        end
    end
    return y
end

#### Sparse Matrix, vector multiplication
function matvec(rowptr::Vector,colval::Vector,nzval::Vector,x::Vector)
    y = zeros(maximum(colval))
    #@show size(y,1)
    #@show size(y,2)
    for i=1:maximum(colval)
        for k=rowptr[i]:rowptr[i+1]-1
            if rowptr[i] == colval[k]
                y[i] += (1 - nzval[k]) * x[colval[k]] ### Identity matrix
            else
                y[i] += nzval[k] * x[colval[k]]
            end
        end
    end
    return y
end

####
# Initialize variables
####
niter = 100 #max iterations
rel = 1e-5 #convergence value
v = ones(size(P,1))
v0 = v/sum(v)
@show length(colptr)
#@time Test = matvec(rowptr,colval,nzval,v0)
#@time Test = cscmatvec(rowval,colptr,nzval,v0,length(v))
#@show sum(v0) # test sum(x)=1
alphar=0.25
y = zeros(size(P,1))
#matvec((I-alpha*P),v)
#r = (1-alpha)*v + matvec((I-alpha*P),v)
####

##### 
# Page rank with richardson iteration
#####
##### Beta is veriable for pagerank
##### Alpha is the variable for richardson
function RankCSR(alphar::Float64,beta::Float64,rel::Float64,rowptr::Vector,colval::Vector,nzval::Vector,niter::Integer,v0::Vector)
    v = v0
    x = randn(size(v,1))
    n = maximum(colval)
    nzval = beta*nzval ### Add beta term
    r = (1-beta)*v - matvec(rowptr,colval,nzval,x)
    iter = 1
    while norm(r) >= n*rel && iter <= niter # set error to 1e5
        r = (1-beta)*v - matvec(rowptr,colval,nzval,x)
        x += alphar*r
        iter+=1
        @show norm(r)
        @show iter
    end
    return x
end
alphar = 1.0
beta = 0.85
@time BB = RankCSR(alphar,beta,rel,rowptr,colval,nzval,niter,v0)



#### page rank CSC
#colptr,rowval,nzval
function RankCSC(alphar::Float64,beta::Float64,rel::Float64,colptr::Vector,rowval::Vector,nzval::Vector,niter::Integer,v0::Vector)
    v = v0
    x = randn(size(v,1))
    n = maximum(rowval)
    nzval = beta*nzval ### Add beta term
    r = (1-beta)*v - cscmatvec(colptr,rowval,nzval,x,maximum(rowval))
    iter = 1
    while norm(r) >= n*rel && iter <= niter # set error to 1e5
        r = (1-beta)*v - cscmatvec(colptr,rowval,nzval,x,maximum(rowval))
        x += alphar*r
        iter+=1
        @show norm(r)
        @show iter
    end
    return x
end
alphar = 0.5
beta = 0.85
#@time CC = Richardson(alpha,beta,rel,colptr,rowval,nzval,niter,v0)

# Normal pagerank function for comparison
function RRank(alphar::Float64,beta::Float64,rel::Float64,P::SparseMatrixCSC,niter::Integer,x0::Vector)
    v = v0
    x = randn(size(v,1))
    r = (1-beta)*v - (I-beta*P)*x
    iter = 1
    n = size(P,1)
    while norm(r) >= n*rel && iter <= niter # set error to 1e5
        r = (1-beta)*v - (I-beta*P)*x
        x += alphar*r
        iter+=1
        #@show norm(r)
        #@show iter
    end
    return x
end

#@time AA = RRank(alpha,beta,rel,P,niter,v0)


