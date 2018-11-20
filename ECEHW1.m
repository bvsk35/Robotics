clc;
clear all;

%QUESTION 1
v1 = [-1;2;4];
v2 = [-2;4;-3;3];
v3 = [-1;4;-3];
v4 = [0;-2;1;5];
%v1+v2 Can't be computed
v1+v3
%dot(v2,v3) Can't be computed
dot(v1,v3)
cross(v1,v3)
norm(v4)
acos(dot(v2,v4)/(norm(2)*norm(v4)))

%QUESTION 2
A = [-4,5,4;-2,-4,-2;-3,4,2];
B = [-4,-5,3;5,-5,-4;-3,1,-1;3,2,-4];
C = [3,1,-2;1,-2,3;0,-2,-3];
%A*B Can't be computed
%A+B Can't be computed
A+C
B*C
%C*B Can't be computed
A*A
%B*B Can't be computed

%QUESTION 3
%Simple calculation of eigen vectors and values of given matrix
%E eigen value matirx and V eigen vector matrix 
D = [3,-1,-4;6,-2,-6;2,-1,-3];
[V,E] = eig(D)

%QUESTION 4
%See pdf
%Simply calculate Laplace inverse of sI-A to get
%solution of the ODE with given initial conditions 
%e^A = L^-1((sI-A)^-1)
syms s;
AA = [-3,-2;4,3];
IA = ilaplace(inv(s*eye(2)-AA)) %Matrix computation for e^A
X0 = [-3;1];
X = IA*X0

%QUESTION 5
%See pdf
%Below Calculation is for C part where we need to add
%Matrices A and B
%Find Eigen values and vectors to diagonalize it
%Then e^(A+B) will be exp raised to the power of 
%eigen values e^(A+B) = P*D*(P^-1) D eigen values diagonal matrix
syms s;
F = [0,1;-4,-4];
G = [1,0;-4,-1];
H = F+G;
[Z,I] = eig(H)
% Another wayt to calculate of e^(A+B) using method from method from question 4
IH = ilaplace(inv(s*eye(2)-H))                               
%Approach 2 of pdf
%Simply calculate Laplace inverse like 4 question
%then substitute t=1 to get e^A and e^B
syms s;
AC = [0,1;-4,-4];
IC = ilaplace(inv(s*eye(2)-AC))
AD = [1,0;-4,-1];
ID = ilaplace(inv(s*eye(2)-AD))
IC*ID %e^A * e^B
ID*IC %e^B * e^A

