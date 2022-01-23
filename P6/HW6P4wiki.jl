using DelimitedFiles
using SparseArrays
using Plots, LinearAlgebra,Krylov
##########################
# PART 3
##########################
using SparseArrays, Random
Random.seed!(0)
P = sprand(1_000,1_000, 15/1_000) |>
    A->begin fill!(A.nzval,1); A; end |>
    A->begin ei,ej,ev = findnz(A); d = sum(A;dims=2);
        return sparse(ej,ei,ev./d[ei], size(A)...); end
### Pre proc
using ZipFile, DelimitedFiles,SparseArrays
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
nmax = size(P,1)
v = ones(nmax)/nmax;
alpha = 0.5
A = (I-alpha*P);
b = (1-alpha)*v;
### END preproc

c, hist_cmin = Krylov.lsqr(A, b,rtol = 1e-8,itmax = 1000, history = true)
#x, hist_pmin = Krylov.minres(A,b,rtol = 1e-8,itmax = 1000, history = true)
#CG1 = Krylov.cg(A,b,rtol = 1e-8, history = true)

#pmin1 = hist_pmin.residuals;
cmin1 = hist_cmin.residuals;
# maxiter = 100
# X = zeros(nmax,maxiter)
# #xnew = minres!(X[:,1], T'*T, b)
# rvec = zeros(maxiter)
# for i = 1:maxiter-1
#     rvec[i] = norm(A*X[:,i] - b)/norm(b) # relative 2 norm risd
#     X[:,i+1] = Krylov.lsqr(X[:,i], A, b)
#     print(i) 

# end
# plot(1:maxiter,rvec)