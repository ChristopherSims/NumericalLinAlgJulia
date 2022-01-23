using SparseArrays
using BlockBandedMatrices
using BandedMatrices
using Latexify
using Plots
#using Kron
# Part 1
n = 10
D2 = BandedMatrix(0 => Fill(-2,n), 1 => Fill(1,n-1), -1 => Fill(1,n-1))
D_xx = BandedBlockBandedMatrix(kron(D2, Eye(n)))
D_yy = BandedBlockBandedMatrix(kron(Eye(n), D2))
L = D_xx + D_yy
L = convert(Matrix,L)
#print("\n")
#print("Returns =  \n")
#show(stdout, "text/plain", L)
# ULax = latexify(L)
# @show ULax
function laplacian(n::Integer, f::Function)
    N = (n+1)^2
    nz = n+1
    I = zeros(N,nz)
    J = zeros(N,nz)
    V = zeros(nz)
    fvec = zeros(N)
    # the transpose mirrors the row indexing we had before.
    G = reshape(1:N, n+1, n+1)' # index map, like we saw before;
    h = 1.0/(n)
    index = 1
    for i=0:n
        for j=0:n
            row = G[i+1,j+1]
            if i==0 || j == 0 || i == n || j == n
                # we are on a boudnary
                fvec[row] = 0.0
                # fill in entries in I,J,V and update index
                I[row,j+1] = index
                J[row,j+1] = index
                V[j+1] = 0.0
                index = index + 1
            else
                fvec[row] = f(i*h, j*h)*h^2
                # fill in entries in I,J,V and update index
                I[row,j+1] = index
                J[row,j+1] = index
                V[j+1] = 0.0
                index = index + 1
            end
        end
    end
    A = sparse(I,J,V,N,N)
    return A, fvec
end

function testf(x,y)
    A = 1
    return A
end

A,fvec = laplacian(10, testf)


print("\n")
fvec = ones(n*n)
print(size(fvec))
u = L\fvec
S = reshape(u,n,n)
surface(S)
savefig("laplaceM")