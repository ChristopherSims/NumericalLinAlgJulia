using LinearAlgebra, Random,SparseArrays, Plots
#### Constant values
L1 = 0.1;
Ln = 100;
rho = 0.9;
n = 30;
#### Function taken from class notes
function lanczos(A,b,k)
    n = size(A,1)
    V = zeros(n,k+1)
    T = Tridiagonal(zeros(k-1), zeros(k), zeros(k-1))
    rho = 0.0
    V[:,1] = b/norm(b)
    for j=1:k
      y = A*V[:,j]
      for i=max(j-1,1):j
        T[i,j] = V[:,i]'*y
        y -= T[i,j]*V[:,i]
      end
      rho = norm(y)
      V[:,j+1] = y/rho
      if j < k
        T[j+1,j] = rho
      end
    end
    return V,T,rho
  end
  #V,T,rho = lanczos(A,b,4)

n = 100
D = zeros(n)
A = zeros(n,n)
e1 = ones(n)

for i = 1:100
    D[i] = L1 + (i-1)/(n-1)*(Ln - L1)*(rho ^ (n-1))
    A[i,i] = D[i]
end

v1 = e1./sqrt(n);

V,T,rho = lanczos(A,v1,61)
soln1 = zeros(30) # for part 1
soln2 = zeros(30)# for part 2
soln3 = zeros(30)# for part 2
for i = 1:30
    soln1[i] = log(norm(V[:,i]'*V[:,i] - I) + 1e-20)
    soln2[i] = log(norm(V[:,1]'*V[:,i] - I) + 1e-20)
    if i >= 3
        soln3[i] = log(norm(V[:,i-2]'*V[:,i] - I) + 1e-20)
    end
end
soln4 = zeros(60) # for part 4
T = Matrix(T)
for i = 1:60
    soln4[i] = log(norm(A*V[:,i] - V[:,i+1]*T[i+1,i+1]) + 1e-20)
end

plot(1:30,soln1, label = "log(V^TV-I)", lw = 3)
xlabel!("K iter")
savefig("P2soln1.png")

plot(1:30,soln2, label = "log(V1^TVk-I)", lw = 3)
plot!(1:30,soln3, label = "log(V(k-2)^TVk-I)", lw = 3)
xlabel!("K iter")
savefig("P2soln2.png")

plot(1:60,soln4, label = "log(A*V[:,i]' - V[:,i+1]T[:,i+1])", lw = 3)
xlabel!("K iter")
savefig("P2soln4.png")


