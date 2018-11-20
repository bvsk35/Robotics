clc;
clear all;

%Question 1
zeta1 = [-0.2;0.5;0.1;-0.4;0.3;-0.1];
a1 = [zeta1(4,1);zeta1(5,1);zeta1(6,1)];
w1 = a1/norm(a1);
b1 = [zeta1(1,1);zeta1(2,1);zeta1(3,1)];
v1 = b1/norm(a1);
w1_hat = [0 -w1(3,1) w1(2,1);w1(3,1) 0 -w1(1,1);-w1(2,1) w1(1,1) 0];
zeta1_hat = [w1_hat v1;0 0 0 0];
theta1 = norm(a1);
R1 = eye(3) + (w1_hat*sin(theta1)) + ((w1_hat^2)*(1-cos(theta1))); %Rodrigues formula for Rotational Matrix
P1 = (eye(3)-R1)*(cross(w1,v1)) + (w1*w1'*v1*theta1); %Comes from expansion of exp(zeta_hat*theta)
g1 = [R1 P1;0 0 0 1]; %ANS - Matrix exponential of zeta_hat is nothing but Rigid Body Transformation Matrix 

%Question 2
g2 = [-0.573973 0.312409 -0.756938 1.57652;
        0.520871 -0.573973 -0.631861 -1.28446;
        -0.631861 -0.756938 0.16672 0.875234;
        0 0 0 1];
R2 = [g2(1,1) g2(1,2) g2(1,3);g2(2,1) g2(2,2) g2(2,3);g2(3,1) g2(3,2) g2(3,3)];
P2 = [g2(1,4);g2(2,4);g2(3,4)]; %Because every RBT g can be written as [R P;0 0 0 1]
theta2 = acos((trace(R2)-1)/2); %For calulting w and theta if R or exp(w_hat*theta) is given
a2(1,1) = (0.5/sin(theta2))*(R2(3,2)-R2(2,3));
a2(2,1) = (0.5/sin(theta2))*(R2(1,3)-R2(3,1));
a2(3,1) = (0.5/sin(theta2))*(R2(2,1)-R2(1,2));
w2 = a2/norm(a2);%ANS - part of the screw
w2_hat = [0 -w2(3,1) w2(2,1);w2(3,1) 0 -w2(1,1);-w2(2,1) w2(1,1) 0];
A2 = ((eye(3) - R2)*w2_hat) + (w2*w2'*theta2);
v2 = A2\P2; %Comes from the expanison of exp(zeta_hat*theta)
h2 = (w2'*v2)/(norm(w2)^2); %ANS - All the below formulas come from screw attached to the zeta 
q2 = (cross(w2,v2))/(norm(w2)^2); %ANS - Screw S = ({q,w},h,theta) theta is magnitude M
M2 = theta2;%ANS
d2 = h2*theta2;
zeta2 = [v2;w2];
zeta2_hat = [w2_hat v2;0 0 0 0]; %Twist
exp_coordinates2 = zeta2*theta2; %ANS - Exponential coordiantes of g is nothing but zeta*theta

%Question 3
h3 = 1.4;
theta3 = 140*(pi/180);
a3 = [-12;20;9];
w3 = a3/norm(a3);
w3_hat = [0 -w3(3,1) w3(2,1);w3(3,1) 0 -w3(1,1);-w3(2,1) w3(1,1) 0];
q3 = [-2;3;1];
R3 = eye(3) + (w3_hat*sin(theta3)) + ((w3_hat^2)*(1-cos(theta3)));
P3 = ((eye(3) - R3)*q3) + (h3*theta3*w3);
g3 = [R3 P3;0 0 0 1]; %ANS - Required RBT
v3 = (-cross(w3,q3)) + (h3*w3);
zeta3 = [v3;w3];
zeta3_hat = [w3_hat v3;0 0 0 0];
P31 = (eye(3)-R3)*(cross(w3,v3)) + (w3*w3'*v3*theta3); %To show that v = -wxq + hw
g31 = [R3 P31;0 0 0 1]; %Coming from the defination of screw or twist will lead to same RBT g

%Question 4
% a4 = sin(pi/4);
% b = pi/4;
% g4_45 = [(4*a4+2)/9 -(5*a4+4)/9 (-7*a4+4)/9 (-12*b+29*a4+13)/9;
%             (8*a4-2)/9 (-a4+4)/9 -(5*a4+4)/9 (12*b+13*a4-31)/9;
%             (8*a4+1)/9 (8*a4-2)/9 (4*a4+2)/9 -(6*b+32*a4-20)/9;
%             0 0 0 1];
b4 = pi/2;
g4_90 = [8/9 -4/9 1/9 (19-12*b4)/9;
            4/9 7/9 -4/9 (-40+12*b4)/9;
            1/9 4/9 8/9 -(10+6*b4)/9;
            0 0 0 1];
g4_0 = [0 -1 0 4;0 0 -1 -1;1 0 0 2;0 0 0 1];
g4 = g4_90/(g4_0); %g4_45 is for theta = pi/4 will also lead to same screw description 
R4 = [g4(1,1) g4(1,2) g4(1,3);g4(2,1) g4(2,2) g4(2,3);g4(3,1) g4(3,2) g4(3,3)];
P4 = [g4(1,4);g4(2,4);g4(3,4)];
theta4 = acos((trace(R4)-1)/2);    
a4(1,1) = (0.5/sin(theta4))*(R4(3,2)-R4(2,3));
a4(2,1) = (0.5/sin(theta4))*(R4(1,3)-R4(3,1));
a4(3,1) = (0.5/sin(theta4))*(R4(2,1)-R4(1,2));
w4 = a4/norm(a4); %ANS - part of the screw i.e. axis of rotation 
w4_hat = [0 -w4(3,1) w4(2,1);w4(3,1) 0 -w4(1,1);-w4(2,1) w4(1,1) 0];
A4 = ((eye(3) - R4)*w4_hat) + (w4*w4'*theta4);
v4 = A4\P4;
h4 = (w4'*v4)/(norm(w4)^2); %ANS - Screw S = ({q,w},h,theta) theta is magnitude M
q4 = (cross(w4,v4))/(norm(w4)^2); %ANS
M4 = theta4; %ANS
d4 = h4*theta4;
zeta4 = [v4;w4];
zeta4_hat = [w4_hat v4;0 0 0 0];

%Question 5
g5 = [0.327697 0.229144 0.916574 0.632756;
        -0.229144 0.960453 -0.158189 -1.26669;
        -0.916574 -0.158189 0.367244 1.23325;
        0 0 0 1];
R5 = [g5(1,1) g5(1,2) g5(1,3);g5(2,1) g5(2,2) g5(2,3);g5(3,1) g5(3,2) g5(3,3)];
P5 = [g5(1,4);g5(2,4);g5(3,4)];
theta5 = acos((trace(R5)-1)/2);
a5(1,1) = (0.5/sin(theta5))*(R5(3,2)-R5(2,3));
a5(2,1) = (0.5/sin(theta5))*(R5(1,3)-R5(3,1));
a5(3,1) = (0.5/sin(theta5))*(R5(2,1)-R5(1,2));
w5 = a5/norm(a5); %ANS - part of the screw i.e. axis of rotation 
w5_hat = [0 -w5(3,1) w5(2,1);w5(3,1) 0 -w5(1,1);-w5(2,1) w5(1,1) 0];
A5 = ((eye(3) - R5)*w5_hat) + (w5*w5'*theta5);
v5 = A5\P5;
h5 = (w5'*v5)/(norm(w5)^2); %ANS - Screw S = ({q,w},h,theta) theta is magnitude M
q5 = (cross(w5,v5))/(norm(w5)^2); %ANS
M5 = theta5; %ANS
zeta5 = [v5;w5];
zeta5_hat = [w5_hat v5;0 0 0 0]; %Twist 
exp_coordinates5 = zeta5*theta5; %Exponential coordinates with screw combined describe the geometric description of the RBT

