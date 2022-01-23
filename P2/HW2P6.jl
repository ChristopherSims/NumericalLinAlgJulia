using LinearAlgebra
using Plots
using SparseArrays
using Random
#include("HW1P9.jl")
include("SparseMatrixCSR.jl")


#####################
# Functions for the problem
##################
function cscmatvec(TI,TJ,TV,x,m)
  y = zeros(m) # allocate output
  for nzi in 1:length(TI) # for each entry in T
    i = TI[nzi]
    j = TJ[nzi]
    v = TV[nzi]
    y[i] += v*x[j]
  end
  return y
end
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

############### 
# Create Variables
###################
m = 5
n = 5
M = rand(m,n)*5
S = sparse(M)
@time PR = SparseMatrixCSR(S) #### CSR is better for Lin alg
rowptr,colval,nzval = PR.rowptr, PR.colval, PR.nzval
sb = ones(m)
colptr = S.colptr
rowval = S.rowval
nzval = S.nzval
display(colptr)
display(rowval)
display(nzval)
niter = 100
x0 = randn(n)
x = copy(x0)
i=2
#ABC = csc_lookup(colptr,rowval,nzval,m,n,i,i)
#@show ABC
##### Gradient Desent, Part 6
#@time r = S*x - sb
#@show r

#ABC = cscmatvec(rowval,colptr,nzval,x,m)



############ Iterate over n
# Part 1
#######################
for i = 1:n
  global r = matvec(rowptr,colval,nzval,x) - sb
  global iter = 1
  while norm(r[i]) >= 1e-5 && iter <= niter # set error to 1e5
    global r = matvec(rowptr,colval,nzval,x) - sb
    temp = csc_lookup(colptr,rowval,nzval,m,n,i,i)
    #@show temp
    global x[i] += -r[i]/temp
    #@show x
    #@show temp
    global iter += 1
    #@show iter
  end
  #@show i
end
@show x


x0 =  ones(n)
x = copy(x0)
#@show x

############ random n
# Part 2
#######################

perm = randperm(n)[1:n]
for i = perm
  global r = matvec(rowptr,colval,nzval,x) - sb
  global iter = 1
  while norm(r[i]) >= 1e-5 && iter <= niter# set error to 1e-55
    global r = matvec(rowptr,colval,nzval,x) - sb
    temp = csc_lookup(colptr,rowval,nzval,m,n,i,i)
    #@show temp
    global x[i] += -r[i]/temp
    #@show x
    #@show temp
    global iter += 1
  end
  #@show i
end
@show x



#@show randperm(10)[1:10]
# ##### Gradient Desent, Part 6 Check function
# r = S*x - sb
# iter = 1
# while norm(r) >= n*1e-5 || iter <= niter # set error to 1e5
#   i = mod(iter,n)+ 1
#   global r = S*x -sb
#   global x[i] += -r[i]/nzval[i]
#   global iter+=1
# end
# @show x

############# Pick random I
# function randCD(colptr::Vector,rowval::Vector,rowptr::Vector,colval::Vector,nzval::Vector,x::Vector,niter::Integer,n::Integer)
#   i = rand(1:n)
#   r = matvec(rowptr,colval,nzval,x) - sb
#   iter = 1
#   if norm(r[i]) > 1e-5 
#     while norm(r[i]) >= 1e-5 && iter <= niter# set error to 1e-55
#       r = matvec(rowptr,colval,nzval,x) - sb
#       temp = csc_lookup(colptr,rowval,nzval,m,n,i,i)
#       x[i] += -r[i]/temp
#       #@show temp
#       iter += 1
#     end
#   end
#   @show norm(r)
#   if norm(r) > n*1e-5
#     randCD(colptr,rowval,rowptr,colval,nzval,x,niter,n)
#   else
#     return x
#   end
# end