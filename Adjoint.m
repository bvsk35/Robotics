function [A] = Adjoint(g)
    R = g(1:3,1:3);
    P = g(1:3,4);
    P_hat = [0 -P(3,1) P(2,1);P(3,1) 0 -P(1,1);-P(2,1) P(1,1) 0];
    A = [R P_hat*R;zeros(3,3) R];

