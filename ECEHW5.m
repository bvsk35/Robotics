clc;
clear all;

% Question 1
gab1 = [0.528383 0.758351 -0.381725 3.32152;
            -0.150539 0.526175 0.836945 3.00205;
            0.835553 -0.384763 0.392184 2.63034;
            0 0 0 1];
zeta1b = [-3;5;-1;-4;2;3];        
W1b = [4;-3;5;-1;3;2];
Rab1 = gab1(1:3,1:3);
Pab1 = gab1(1:3,4);
Pab1_hat = [0 -Pab1(3,1) Pab1(2,1);Pab1(3,1) 0 -Pab1(1,1);-Pab1(2,1) Pab1(1,1) 0];
Adgab1 = [Rab1 Pab1_hat*Rab1;zeros(3,3) Rab1];
Adgba1 = [Rab1' -Rab1'*Pab1_hat;zeros(3,3) Rab1'];
zeta1a = Adgab1*zeta1b; % Answer
W1a = Adgba1'*W1b; % Answer

% Question 2
% subs(PST2_t,{theta1,theta2,theta3},[0,0,0])
% Hat of zeta is function which computes hat of the zeta, where  twist = zeta = [-cross(w,q);w]
% zeta are the twist which are sent as argument to the function 
% here I sent zeta*theta to compute hat of it 
syms L1 L2 L3 theta1 theta2 theta3 alpha real
zeta21 = hatofzeta([0;0;0;0;0;1]*theta1);
zeta22 = hatofzeta([0;-L1;0;-1;0;0]*theta2);
zeta23 = hatofzeta([(L2*sin(alpha)-L1*cos(alpha));0;0;0;cos(alpha);sin(alpha)]*theta3);
P2 = [0;L2+L3;L1];
gST2_0 = [eye(3,3) P2;0 0 0 1];
gST2_t = expm(zeta21)*expm(zeta22)*expm(zeta23)*gST2_0;
RST2_t = simplify(gST2_t(1:3,1:3),2); % Answer
PST2_t = vpa(simplify(gST2_t(1:3,4)),2); % Answer

% Question 3
syms L1 L2 L3 theta1 theta2 theta3 real
zeta31 = hatofzeta([0;0;0;0;0;1]*theta1);
zeta32 = hatofzeta([0;L1;0;1;0;0]*theta2);
zeta33 = hatofzeta([0;L1+L2;0;1;0;0]*theta3);
P3 = [0;0;L1+L2+L3];
gST3_0 = [eye(3,3) P3;0 0 0 1];
gST3_t = expm(zeta31)*expm(zeta32)*expm(zeta33)*gST3_0;
RST3_t = vpa(gST3_t(1:3,1:3),2); % Answer
PST3_t = vpa(simplify(gST3_t(1:3,4)),2); % Answer

% Question 5
syms L1 L2 L3 theta1 theta2 theta3 theta4 theta5 theta6 real
zeta51 = hatofzeta([0;0;0;0;0;1]*theta1);
zeta52 = hatofzeta([0;-L1;0;-1;0;0]*theta2);
zeta53 = hatofzeta([0;-L1;L2;-1;0;0]*theta3);
zeta54 = hatofzeta([-L1;0;0;0;1;0]*theta4);
zeta55 = hatofzeta([0;-L1;L2+L3;-1;0;0]*theta5);
zeta56 = hatofzeta([-L1;0;0;0;1;0]*theta6);
P5 = [0;L2+L3;L1];
gST5_0 = [eye(3,3) P5;0 0 0 1];
gST5_t = expm(zeta51)*expm(zeta52)*expm(zeta53)*expm(zeta54)*expm(zeta55)*expm(zeta56)*gST5_0;
RST5_t = vpa(gST5_t(1:3,1:3),2); % Answer
PST5_t = vpa(simplify(gST5_t(1:3,4)),2); % Answer

