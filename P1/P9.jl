using SparseArrays
using Arpack
m = 18
n = 5
I = [1, 4, 3, 5] 
J = [4, 7, 18, 9] 
V = [1, 2, -5, 3]
xn = randn(n)
xm = randn(m)
S = sparse(I,J,V)
#@show S
#print(varinfo(r"S"))


print("\n")
print("\n")
print(S.colptr)
print("\n")
print(S.rowval)
print("\n")
print(S.nzval)
print("\n")
#ST = sparse(S.rowval[1:end-1],S.colptr,S.nzval,n,m)
show(stdout, "text/plain", S')
print("\n")
# function csc2mat(colptr,rowval,nzval,m,n)
#     y = zeros(m,n) 
#     for j=1:length(colptr)-1 # for each column ...
#         for nzi=colptr[j]:colptr[j+1]-1 # for each entry in the column
#             i = rowval[nzi]
#             v = nzval[nzi]
#             y[i,j] = v
#             #y[i] += v*x[j]
#         end
#     end
#     return y
# end

# A = csc2mat(S.colptr,S.rowval,S.nzval,m,n)
#@show A
#show(stdout, "text/plain", A)

#show(stdout, "text/plain", A)




# Problem 9.1  Done 
""" Returns y = A'*x where A is given by the CSC arrays
colptr, rowval, nzval, m, n and x is the vector. """
function csc_transpose_matvec(colptr, rowval, nzval, m, n, x)
    if (length(x) == n)
        S1 = sparse(colptr,rowval,nzval)
        S2 = S1'
        y = S2*x
        return y
    else
        return "Error Matrix mismatch"
    end
end
y = csc_transpose_matvec(I,J,V,m,n,xn)
print("y = A'*x \n")
show(stdout, "text/plain", y)


#Problem 9.2 Done
""" Returns = A[:,i]'*x where A is given by the CSC arrays
colptr, rowval, nzval, m, n and x is the vector. """
function csc_column_projection(colptr, rowval, nzval, m, n, i, x)
    S = sparse(colptr,rowval,nzval)
    if (length(x) == n)
        SF = S[:,i]'*x
        return SF
    else
        return "Error: Dimension Mismatch"
    end
    
end
y2 = csc_column_projection(I,J,V,m,n,3,xn)
print("\n")
print("Returns = A[:,i]'*x \n")
show(stdout, "text/plain", y2)
print("\n")


# Problem 9.3 Done
""" Returns rho = A[:,i]'*A[:,j] where A is given by the CSC arrays
colptr, rowval, nzval, m, n and i, and j are the column indices. """
# Is this the norm or just the product????
function csc_col_col_prod(colptr, rowval, nzval, m, n, i, j)
    S = sparse(colptr,rowval,nzval) 
    B = S[:,j] #A[:,j]
    A = S[:,i]' # A[:,i]
    T = A*B
    return T
end
y3 = csc_col_col_prod(I,J,V,m,n,2,3)
print("\n")
print("rho = A[:,i]'*A[:,j]\n")
show(stdout, "text/plain", y3)
print("\n")

#Problem 9.4 Done
""" Returns rho = A[i,j] where A is given by the CSC arrays
colptr, rowval, nzval, m, n and i, and j are the column indices. """
function csc_lookup(colptr, rowval, nzval, m, n, i, j)
    if (i < m) || (j < m)
        B=0
        for jj=1:length(colptr)-1 # for each column
            for nzi=colptr[jj]:colptr[jj+1]-1 # for each entry in the column
                ii = rowval[nzi]
                v = nzval[nzi]
                if (ii == i) && (jj == j)
                    B = v
                end
            end
        end
        return B
    else
        return "Error out of bounds"
    end
    
end
# Test CSC Lookup
B = csc_lookup(I,J,V,m,n,4,1)
print("\n")
show(stdout, "text/plain", B)
print("\n")


# Problem 9.5 Done
""" Returns x = A[i,:] where A is given by the CSC arrays
colptr, rowval, nzval, m, n and i is the row index . """
function csc_lookup_row(colptr, rowval, nzval, m, n, i)
    A = sparse(colptr,rowval,nzval) 
    if i < m
        X = A[i,:]
        
        return Vector(X)
    else
        return "Error: out of bounds"
    end
end
x = csc_lookup_row(I,J,V,m,n,2)
print("\n")
show(stdout, "text/plain", x)
print("\n")