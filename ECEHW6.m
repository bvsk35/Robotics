clc;
clear all;

% Question 1
% Hat of zeta is function which computes hat of the zeta, where  twist = zeta = [-cross(w,q);w]
% zeta are the twist which are sent as argument to the function 
% here I sent zeta*theta to compute hat of it 
% Adjoint is a function written by me which computes Adjoint of a given
% 4x4 matrix using simple formula given in textbook. Similarly IAdjoint is
% a function which computes Inverse of adjoint matrix using direct formula
% rather than relying on inv(.) function of MATLAB
syms L1 L2 L3 theta1 theta2 theta3 alpha real
z11 = [0;0;0;0;0;1];
z12 = [0;-L1;0;-1;0;0];
z13 = [(L2*sin(alpha)-L1*cos(alpha));0;0;0;cos(alpha);sin(alpha)];
zeta11 = hatofzeta(z11*theta1);
zeta12 = hatofzeta(z12*theta2);
zeta13 = hatofzeta(z13*theta3);
P2 = [0;L2+L3;L1];
gST1_0 = [eye(3,3) P2;0 0 0 1];
gST1_t = expm(zeta11)*expm(zeta12)*expm(zeta13)*gST1_0;
z_12 = Adjoint(expm(zeta11))*z12;
z_13 = Adjoint(expm(zeta11)*expm(zeta12))*z13;
Js_st1 = [z11 z_12 z_13];
Jb_st1 = IAdjoint(gST1_t)*Js_st1; % Answer 
% Last column of body jacobian can be computed to verify if our answer is
% correct or wrong by applying text book formula
% K = IAdjoint(expm(zeta13)*gST1_0)*z13; 

% Question 2 & 3
syms L1 L2 L3 theta1 theta2 theta3 real
z21 = [0;0;0;0;0;1];
z22 = [0;L1;0;1;0;0];
z23 = [0;L1+L2;0;1;0;0];
zeta21 = hatofzeta(z21*theta1);
zeta22 = hatofzeta(z22*theta2);
zeta23 = hatofzeta(z23*theta3);
z_22 = Adjoint(expm(zeta21))*z22;
z_23 = Adjoint(expm(zeta21)*expm(zeta22))*z23;
Js_st2 = [z21 z_22 z_23]; % Answer for Q2
Det2 = simplify(det(Js_st2'*Js_st2)); 
solve(Det2 == 0, theta2); 
r2 = rank(Js_st2);
a2 = -2*pi:0.01:2*pi;
A2 = subs(Det2,{L1,L2},[1 1]);
A2 = subs(A2,{theta2},a2);
plot(a2,A2); 
grid on; 
% No singularity exist since equation Det2 can't be solved - Answer
% Graph between Det2 v/s theta2, Det2 only depends on theta2 and L1 = L2 = 1
% We can observe from the graph that Det2 never becomes zero
% hence no singularites - Answer for Q3

% Question 4
syms L0 L1 L2 D3 theta1 theta2 theta3 theta4 theta5 theta6 real
z41 = [0;0;0;0;0;1];
z42 = [0;L0;0;1;0;0];
z43 = [0;1;0;0;0;0];
z44 = [-L0;0;L1;0;1;0];
z45 = [0;L0;-D3;1;0;0];
z46 = [-L0;0;L1;0;1;0];
zeta41 = hatofzeta(z41*theta1);
zeta42 = hatofzeta(z42*theta2);
zeta43 = hatofzeta(z43*theta3);
zeta44 = hatofzeta(z44*theta4);
zeta45 = hatofzeta(z45*theta5);
zeta46 = hatofzeta(z46*theta6);
z_42 = Adjoint(expm(zeta41))*z42;
z_43 = Adjoint(expm(zeta41)*expm(zeta42))*z43;
z_44 = Adjoint(expm(zeta41)*expm(zeta42)*expm(zeta43))*z44;
z_45 = Adjoint(expm(zeta41)*expm(zeta42)*expm(zeta43)*expm(zeta44))*z45;
z_46 = Adjoint(expm(zeta41)*expm(zeta42)*expm(zeta43)*expm(zeta44)*expm(zeta45))*z46;
Js_st4 = [z41 z_42 z_43 z_44 z_45 z_46]; % Answer with entire Jacobian

% Question 5
% This solution is considering only last three joints i.e wrist and
% computing the Jacobians and the singularites
% Frame S & T are placed at the wrist and both are coincedent
syms theta1 theta2 theta3 real
z51 = [0;0;0;0;0;1];
z52 = [0;0;0;-1;0;0];
z53 = [0;0;0;0;1;0];
zeta51 = hatofzeta(z51*theta1);
zeta52 = hatofzeta(z52*theta2);
zeta53 = hatofzeta(z53*theta3);
z_52 = Adjoint(expm(zeta51))*z52;
z_53 = Adjoint(expm(zeta51)*expm(zeta52))*z53;
Js_st5 = [z51 z_52 z_53];
Det5 = simplify(det(Js_st5'*Js_st5));
solve(Det5==0,theta2);
a5 = -2*pi:0.01:2*pi;
A5 = subs(Det5,{theta2},a5);
figure;
plot(a5,A5);
grid on;
% Graph between Det5 v/s theta2, Det5 only depends on theta2
% We can observe from the graph that Det5 becomes zero whenver theta2 is
% equal to pi/2, hence there is singularity - Answer
% See below for even more possible singularity

% % Question 5
% % This soultion is where I consider all the twist and compute the Jacobian
% % and the frame S is placed at the base and the frame T is placed at the
% % wrist and the refernce configuration is as shown in figure also computed
% % the det of the spatial Jacobian for refernce. Here I also show how all
% % the text book examples of the singularites are possible
% % A61 is the singularity occuring becasue of the two collinear axis, this
% % happens when joint 5 is rotate by pi/2 rads then joint 4 & joint 6 become
% % collinear, so I kept all joints at reference config and rotated just
% % joint 5 and we come observe how joint 4 & joint 6 axis are collinear.
% % Here we are observing the Jacobian
% % A62 is the singularity occuring becasue of the three parallel coplanar
% % axis this occurrs at refernce config so I kept everything at refernce
% % config and we can observe axis of joint 2,3& 5 are coplanar. Here we
% % are observing the Jacobian
% % A63 is singularity when the joint 3 is rotated by pi radians then all the
% % wrist and joint 1,2 axis become 4 intersecting axis. Here we are
% % observing the Jacobian
% % Det6 = simplify(det(Js_st6)); is equal to the value given below
% % -L2*L3*cos(theta5)*(L2*cos(theta2)*sin(theta3) - L3*sin(theta2) + L3*cos(theta3)^2*sin(theta2) + L3*cos(theta2)*cos(theta3)*sin(theta3))
syms L1 L2 L3 theta1 theta2 theta3 theta4 theta5 theta6 real
z61 = [0;0;0;0;0;1];
z62 = [0;-L1;0;-1;0;0];
z63 = [0;-L1;L2;-1;0;0];
z64 = [L2+L3;0;0;0;0;1];
z65 = [0;-L1;L2+L3;-1;0;0];
z66 = [-L1;0;0;0;1;0];
zeta61 = hatofzeta(z61*theta1);
zeta62 = hatofzeta(z62*theta2);
zeta63 = hatofzeta(z63*theta3);
zeta64 = hatofzeta(z64*theta4);
zeta65 = hatofzeta(z65*theta5);
zeta66 = hatofzeta(z66*theta6);
z_62 = Adjoint(expm(zeta61))*z62;
z_63 = Adjoint(expm(zeta61)*expm(zeta62))*z63;
z_64 = Adjoint(expm(zeta61)*expm(zeta62)*expm(zeta63))*z64;
z_65 = Adjoint(expm(zeta61)*expm(zeta62)*expm(zeta63)*expm(zeta64))*z65;
z_66 = Adjoint(expm(zeta61)*expm(zeta62)*expm(zeta63)*expm(zeta64)*expm(zeta65))*z66;
Js_st6 = [z61 z_62 z_63 z_64 z_65 z_66];
A61 = subs(simplify(Js_st6),{theta1,theta2,theta3,theta4,theta5,theta6},[0 0 0 0 pi/2 0]);
A62 = subs(simplify(Js_st6),{theta1,theta2,theta3,theta4,theta5,theta6},[0 0 0 0 0 0]);
A63 = subs(simplify(Js_st6),{theta1,theta2,theta3,theta4,theta5,theta6},[0 0 pi 0 0 0]);


