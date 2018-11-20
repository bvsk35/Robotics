clc;
clear all;

%Question 1
v = transpose([5 -1 -2]);
v_hat = [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0] %textbook formula for hat operator 
a = norm(v_hat); %not equal to one hence using general formula of Rodrigues
R11 = eye(3) + (v_hat*sin(a)/a) + ((v_hat^2)*(1-cos(a))/(a^2)) %angle of rotation is theta*norm(v_hat)
R12 = expm(v_hat); %not necessary but R12 and R11 are same, above formula is more efficent 
axis_of_rotation = v/a %unit vector of axis of rotation 
angle_of_rotation = a %in radians

%Question 2
w = transpose([4 1 -3]); %axis of rotation
w_hat = [0 -w(3) w(2);w(3) 0 -w(1);-w(2) w(1) 0]; %textbook formula for hat operator 
b = norm(w_hat); %not equal one hence using general formula of Rodrigues
p = transpose([5 2 -4]); %intial vector 
theta = 75*(pi/180); %conerting it to radians 
R21 = eye(3) + (w_hat*sin(b*theta)/b) + ((w_hat^2)*(1-cos(b*theta))/(b^2)); %angle of rotation is theta*norm(w_hat)
R22 = expm(w_hat*theta); %not necessary but R22 and R21 are same, above formula is more efficent 
p_new = R21*p %final vector after rotation 

%Question 3
c = pi/4;
d = pi/3;
Rx = [1 0 0;0 cos(d) -sin(d);0 sin(d) cos(d)];
Rz = [cos(c) -sin(c) 0;sin(c) cos(c) 0;0 0 1];
R31 = Rx*Rz %Euler angle of rotation, gives Rab Frame B relative to Frame A 
R32 = Rz*Rx %Fixed angle of rotation, gives Rab Frame B relative to Frame A but NOT same as above  

%Question 4
R4 = [0.4619 -0.1189 -0.8790;-0.5615 -0.8063 -0.1860;-0.6866 0.5794 -0.4392];
theta1 = acos((trace(R4)-1)/2) %angle can be chosen to be theta1+2*pi*n or -theta+2*pi*n, Angle is in radians 
e = (R4(3,2)-R4(2,3))/(2*sin(theta1));
f = (R4(1,3)-R4(3,1))/(2*sin(theta1));
g = (R4(2,1)-R4(1,2))/(2*sin(theta1));
w1 = [e;f;g] %norm of w1 is 1, hence it is a unit vector of angle of rotation 
exponential_coordinates = w1*theta1 %exp coordinates are components of vetor w*theta1 

%Question 5
R5 = [0.6325 0.2533 -0.7319;-0.7074 0.5737 -0.4128;0.3154 0.7789 0.5421];
qA = [-4;3;5];
qB = transpose(R5)*qA %using direct formula taught from class 



