# Problem 1-1

# M1 = [[2, 3, 5]
#     [7, 11, 13]
#     [17, 19, 23]]
# M2 = [[1, 2, -3]
#     [-4, -5, -6]
#     [7, -8, 9]]

M1 = [2 3 5; 7 11 13; 17 19 23]
M2 = [1 2 -3; -4 -5 -6; 7 -8 9]
#M3 = zeros(3,3)
M3 = M1'*M2
print("P1-1:\n ")
print(M3)
print("\n")
#problem 1-2
x = ones(1000,1)
y = 1:1000
V1 = x'*y
print("The Answer to P1-2 is:\n ")
print(V1)
print(size(x))
print(size(y))
print("\n")
#Problem 1-3
x = [1.5 2 -3]'
e = [1 1 1 1]'
Me1 = e*x'
#@show Me1
Me2 = x*e'
#@show Me2
print("Problem 1-3\n")
print("ex':\n")
print(Me1)
print("\n")
print("xe':\n")
print(Me2)
print("\n")
#Problem 1-4
x = [-5 4 2]'
e1 = [1 0 0]'
e2 = [0 1 0]'

Me3 = e2*x'
Me4 = x*e1'
print("Problem 1-4\n")
print("ex':\n")
print(Me3)
print("\n")
print("xe':\n")
print(Me4)
print("\n")
#Problem 1-5
M51 = [2 3 5; 7 11 13; 17 19 23]
e3 = [0;0;1]
e1 = [1 0 0]'
M5A1 = M51*e3
M5A2 = e1'*M51
#@show M5A1
print("Problem 1-5\n")
print(M5A1)
print("\n")
print(M5A2)
print("\n")
#Problem 1-6
M4 = [2 0 0; 0 0.5 0; 0 0 -1]
M5 = [1 2 3; 10 8 6; -2 -2 -2]
M6 = M4*M5
print("Problem 1-6\n")
print(M6)


