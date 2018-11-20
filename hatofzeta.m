function [a] = hatofzeta(A)
v = A(1:3,1);
w = A(4:6,1);
w_hat = [0 -w(3,1) w(2,1);w(3,1) 0 -w(1,1);-w(2,1) w(1,1) 0];
a = [w_hat v;0 0 0 0];