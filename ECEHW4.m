clc;
clear all;

% %Question 1
syms t real;
theta1 = (3*(t^2)) + (2*t);
theta1_d = 6*t + 2;
q1 = [2;-5;3];
a1 = [20;-9;-12];
w1 = a1/norm(a1); % normalizing for the w direction of the axis
w1_hat = [0 -w1(3,1) w1(2,1);w1(3,1) 0 -w1(1,1);-w1(2,1) w1(1,1) 0];
h1 = 1.3;
gab1_0 = [0.510799 0.147375 0.846974 0.606065;
            0.767512 -0.522024 -0.372044 1.90451;
            0.38731 0.840102 -0.379762 1.12323;
            0 0 0 1];
v1 = -cross(w1,q1) + (h1*w1); % comes from the screw theory 
R1 = eye(3) + (w1_hat*sin(theta1)) + ((w1_hat^2)*(1-cos(theta1))); % Rodriguez formula
P1 = (eye(3)-R1)*(cross(w1,v1)) + (w1*w1'*v1*theta1); % comes from expansion of exp(zeta*theta)
g1 = [R1 P1;0 0 0 1];
gab1_t = vpa(g1*gab1_0,4); % ANSWER vpa is used to repersent the final answer in decimal form otherwise output will look garbage
zeta1_hat = [w1_hat v1;0 0 0 0];
zeta1 = [v1;w1];
Vs1_ab = vpa(zeta1*theta1_d); % ANSWER comes from screw theory spatial vel = zeta*theta_dot
Rab0 = [gab1_0(1,1) gab1_0(1,2) gab1_0(1,3);gab1_0(2,1) gab1_0(2,2) gab1_0(2,3);gab1_0(3,1) gab1_0(3,2) gab1_0(3,3)];
Pab0 = [gab1_0(1,4);gab1_0(2,4);gab1_0(3,4)];
Pab0_hat = [0 -Pab0(3,1) Pab0(2,1);Pab0(3,1) 0 -Pab0(1,1);-Pab0(2,1) Pab0(1,1) 0];
IAd_gab0 = [Rab0' -Rab0'*Pab0_hat;zeros(3) Rab0']; 
Vb1_ab = vpa(IAd_gab0*zeta1*theta1_d,2); % ANSWER comes from screw theory body vel = (Adgb(0))^-1*zeta*theta_dot
% Below calculation is for body velocity using Vb_ab = (Adgab(0))^-1*Vs_ab
% above used and this formula will get same result
% Rabt = [gab1_t(1,1) gab1_t(1,2) gab1_t(1,3);gab1_t(2,1) gab1_t(2,2) gab1_t(2,3);gab1_t(3,1) gab1_t(3,2) gab1_t(3,3)];
% Pabt = [gab1_t(1,4);gab1_t(2,4);gab1_t(3,4)];
% Pabt_hat = [0 -Pabt(3,1) Pabt(2,1);Pabt(3,1) 0 -Pabt(1,1);-Pabt(2,1) Pabt(1,1) 0];
% Ad_gabt = [Rabt' -Rabt'*Pabt_hat;zeros(3) Rabt'];
% Vb1_abt = vpa(Ad_gabt*Vs1_ab,2); 
% This part of the code is remove the negliable part of the final answer
% which are bascially zero but MATLAB fails to remove
% for i = 1:1:6
%     [C,T] = coeffs(simplify(Vb1_abt(i,1)));
%     for j = 1:1:length(C)
%         if(abs(C(j))<0.00000001)
%             C(j) = 0;
%         end
%     end
%     V(i,1) = dot(C,T);
% end
% Vb1_abt = vpa(V,2);
 
% %Question 2
syms t0 real;
gab2_t0 = [0.930114 -0.209659 -0.301549 1.35887;
            -0.209659 0.371024 -0.904646 -2.7234;
            0.301549 0.904646 0.301137 0.738995;
            0 0 0 1];
Vs2_ab = [-1.2;0.2;0.4;-0.6;-1.0;0.8];
b2 = [-1.2;0.2;0.4];
a2 = [-0.6;-1.0;0.8];
w2 = a2/norm(a2); % normalizing both v and w obtained from spatial velocity 
v2 = b2/norm(a2);
q2 = cross(w2,v2)/(norm(w2)^2); % ANSWER point on axis, formula comes from screw theory
h2 = w2'*v2/(norm(w2)^2); % ANSWER pitch, formula comes from screw theory
w2_hat = [0 -w2(3,1) w2(2,1);w2(3,1) 0 -w2(1,1);-w2(2,1) w2(1,1) 0];
zeta2_hat = [w2_hat v2;0 0 0 0];
zeta2 = [v2;w2];
R2 = eye(3) + (w2_hat*sin(t0)) + ((w2_hat^2)*(1-cos(t0))); % Rodriguez formula
P2 = (eye(3)-R2)*(cross(w2,v2)) + (w2*w2'*v2*t0); % comes from expansion of exp(zeta*theta)
g2 = [R2' -R2'*P2;0 0 0 1];
gab2_0 = vpa(simplify(g2*gab2_t0),2); % computing initial config frame from given conditions 
% This part of the code is remove the negliable part of the final answer
% which are bascially zero but MATLAB fails to remove
for i = 1:1:4
    for j = 1:1:4
        [C,T] = coeffs(gab2_0(i,j));
        for m = 1:1:length(C)
            if(abs(C(m))<0.00000001)
                C(m) = 0;
            end
        end
        V(i,j) = dot(C,T);
    end
end
gab2_0 = vpa(V,2); % ANSWER inital configruration

% %Question 3
syms t real;
a3 = [-1;-3;1];
ws3_ab = a3; % Saptial angualr velocity 
w3 = a3/norm(a3); % Normalising is it because it is axis of rotation also 
P3 = [-cos(t)*t^2+t*sin(t)+3;sin(t)*t^2+2*t*cos(t)-3*cos(t);t^2-2*sin(t)+3*t*cos(t)-1]; % This is nothing but Pab
P3_d = [(t^2*sin(t))-(t*cos(t))+(sin(t));(t^2*cos(t))+(2*cos(t))+(3*sin(t));(2*t)-(3*t*sin(t))+(cos(t))];
Vs3_ab = [-cross(ws3_ab,P3)+P3_d;ws3_ab]; % ANSWER Defination of spatial velocity 
w3_hat = [0 -w3(3,1) w3(2,1);w3(3,1) 0 -w3(1,1);-w3(2,1) w3(1,1) 0];
R3 = eye(3) + (w3_hat*sin(t)) + ((w3_hat^2)*(1-cos(t))); % Rodriguez formula
Vb3_ab = [R3'*P3_d;R3'*ws3_ab]; % ANSWER Defination of body velocity
R31 = subs(R3,2);
P31 = subs(P3,2);
g3 = vpa([R31 P31;0 0 0 1],6); % ANSWER Subsituting t = 2 
% One more way to calculate Body velocity using Vb = (Adgab)^-1*Vs but
% still will lead same result as above
% P3_hat = [0 -P3(3,1) P3(2,1);P3(3,1) 0 -P3(1,1);-P3(2,1) P3(1,1) 0];
% IADg = [R3' -R3'*P3_hat;zeros(3) R3'];
% Vb31_ab = vpa(simplify(IADg*Vs3_ab),2);

% %Question 4
Vs4_ab = [0.5;-0.2;-0.3;-0.7;-0.4;0.6];
gab4 = [0.544224 0.756286 -0.363114 1.24676;
            -0.821625 0.392997 -0.412899 -0.593889;
            -0.169567 0.523052 0.835262 -1.16116;
            0 0 0 1];
R4 = [gab4(1,1) gab4(1,2) gab4(1,3);gab4(2,1) gab4(2,2) gab4(2,2);gab4(3,1) gab4(3,2) gab4(3,3)];
P4 = [gab4(1,4);gab4(2,4);gab4(3,4)];
P4_hat = [0 -P4(3,1) P4(2,1);P4(3,1) 0 -P4(1,1);-P4(2,1) P4(1,1) 0];
IAd_gab4 = [R4' -R4'*P4_hat;zeros(3) R4'];
Vb4_ab = IAd_gab4*Vs4_ab; % ANSWER Direct formula Vb = (Adgab)^-1*Vs

% %Question 5
syms t real;
Vs5_ab = [-2;7;5;-1;-4;8];
Vb5_bc = [6;9;-5;4;7;-1];
gab5_0 = [-0.43432 -0.583051 -0.686599 0.608676;
            -0.876071 0.0961901 0.472491 1.37076;
            -0.209442 0.806721 -0.552571 -1.28077;
            0 0 0 1];
gbc5_0 = [0.191896 0.929215 0.31581 1.15862;
            0.623914 -0.363899 0.691599 0.385134;
            0.757567 0.0643232 -0.649581 1.30157;
            0 0 0 1];
a51 = [Vs5_ab(4,1);Vs5_ab(5,1);Vs5_ab(6,1)];
b51 = [Vs5_ab(1,1);Vs5_ab(2,1);Vs5_ab(3,1)];
w51 = a51/norm(a51); % Calculating zeta and normalizing it from Vs saptial velocity
v51 = b51/norm(a51); % because we know Vs = zeta*theta_dot for screw motion 
R5 = [gbc5_0(1,1) gbc5_0(1,2) gbc5_0(1,3);gbc5_0(2,1) gbc5_0(2,2) gbc5_0(2,3);gbc5_0(3,1) gbc5_0(3,2) gbc5_0(3,3)];
P5 = [gbc5_0(1,4);gbc5_0(2,4);gbc5_0(3,4)];
P5_hat = [0 -P5(3,1) P5(2,1);P5(3,1) 0 -P5(1,1);-P5(2,1) P5(1,1) 0];
Adg5 = [R5 P5_hat*R5;zeros(3) R5];
zeta52 = Adg5*Vb5_bc;
a52 = [zeta52(4,1);zeta52(5,1);zeta52(6,1)];
b52 = [zeta52(1,1);zeta52(2,1);zeta52(3,1)];
w52 = a52/norm(a52); % Calculating zeta and normalizing it from Vb body velocity
v52 = b52/norm(a52); % because we know Vb = (Adgbc(0))^-1*zeta*theta_dot
w51_hat = [0 -w51(3,1) w51(2,1);w51(3,1) 0 -w51(1,1);-w51(2,1) w51(1,1) 0];
w52_hat = [0 -w52(3,1) w52(2,1);w52(3,1) 0 -w52(1,1);-w52(2,1) w52(1,1) 0];
R51 = eye(3) + (w51_hat*sin(t)) + ((w51_hat^2)*(1-cos(t)));
R52 = eye(3) + (w52_hat*sin(t)) + ((w52_hat^2)*(1-cos(t)));
P51 = (eye(3)-R51)*(cross(w51,v51)) + (w51*w51'*v51*t);
P52 = (eye(3)-R52)*(cross(w52,v52)) + (w52*w52'*v52*t);
P51_hat = [0 -P51(3,1) P51(2,1);P51(3,1) 0 -P51(1,1);-P51(2,1) P51(1,1) 0];
P52_hat = [0 -P52(3,1) P52(2,1);P52(3,1) 0 -P52(1,1);-P52(2,1) P52(1,1) 0];
g51 = [R51 P51;0 0 0 1]; % exp(zeta1*theta1) = gs1
g52 = [R52 P52;0 0 0 1]; % exp(zeta2*theta2) = gs2
gab5_t = vpa(g51*gab5_0,2); % ANSWER Calculating gab(t) = exp(zeta1*theta1) * gab(0)
gbc5_t = vpa(g52*gbc5_0,2); % ANSWER Calculating gbc(t) = exp(zeta2*theta2) * gbc(0)
IAd_gbc5 = [R52' -R52'*P52_hat;zeros(3) R52'];
IAd_gab5 = [R51' -R51'*P51_hat;zeros(3) R51'];
Ad_gab5 = [R51 P51_hat*R51;zeros(3) R51];
Ad_gbc5 = [R52 P52_hat*R52;zeros(3) R52];
Vb5_ac1 = (IAd_gbc5*IAd_gab5*Vs5_ab) + (Vb5_bc); % Direct formula
Vb5_ac = vpa(simplify(Vb5_ac1),2);
% This part of the code is remove the negliable part of the final answer
% which are bascially zero but MATLAB fails to remove
for i = 1:1:6
    [C,T] = coeffs(Vb5_ac(i,1));
    for j = 1:1:length(C)
        if(abs(C(j))<0.00000001)
            C(j) = 0;
        end
    end
    V(i,1) = dot(C,T);
end
Vb5_ac = vpa(V,2); % ANSWER Body velocity b/w frames A,C
Vs5_ac1 = (Vs5_ab) + (Ad_gab5*Ad_gbc5*Vb5_bc); % Direct formula
Vs5_ac = vpa(simplify(Vs5_ac1),2);
% This part of the code is remove the negliable part of the final answer
% which are bascially zero but MATLAB fails to remove
for i = 1:1:6
    [C,T] = coeffs(Vs5_ac(i,1));
    for j = 1:1:length(C)
        if(abs(C(j))<0.00000001)
            C(j) = 0;
        end
    end
    V(i,1) = dot(C,T);
end
Vs5_ac = vpa(V,2); % ANSWER Saptial velocity b/w frames A,C