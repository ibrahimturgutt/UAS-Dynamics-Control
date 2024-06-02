%---------------Control----------
clear all;
close all;
clc;
%% 
%---------------A--------------
A = [-0.0345, 0.847, -2.950, -9.698; -0.513, -3.249, 19.340, -1.479; 0.278, -1.586, -3.233, 0; 0, 0, 1, 0];
B = [0.465; -3.050; -28.910; 0];
eigen = eig(A);
figure;
zplane(eigen);

%% 
%%---------------B--------------
% sX(s) - x(0) = AX(s) + BU(s)
% Y(s) = CX(s) + DU(s)
% sX(s) - AX(s) = BU(s)
% X(s) = (sI-A)^(-1)BU(ss)
% Y(s) = [C(sI-A)^(-1)B+D]U(s)
% G(s) = C(sI-A)^(-1)B + D
syms s

C_angle = [0,0,0,1];
D = [0];

%sys1 is for pitch angle
[num, den] = ss2tf(A,B,C_angle,D);
poles_num = zpk(num);
poles_den = zpk(den);
sys1 = tf(num,den)
K = dcgain(sys1)
figure;
pzmap(sys1);

%sys2 is for pitch rate
C_rate = [0,0,1,0];
[num_rate,den_rate] = ss2tf(A,B,C_rate,D);
sys2 = tf(num_rate,den_rate)
L = dcgain(sys2)
figure;
pzmap(sys2);
%%
%---------------C---------------

figure;
margin(sys1)
figure;
margin(sys2)

figure;
nyquist(sys1)
figure;
nyquist(sys2)
%%
%---------------D-------------

[G_s,G_f] = slowfast(sys1,2)
damp(sys1)
damp(G_s) % long period
damp(G_f) % short period
figure;
margin(G_s);
figure;
margin(G_f);

%%
%---------------E---------------
K = [0.823, 8.84, 16.89];
Cpid = tf(K,[1,0])
CL_sys = feedback(Cpid*sys1,sys1); %cl for pitch angle
OL_Bode = bode(Cpid*sys1);
CL_sys2 = feedback(Cpid*sys2,1); %cl for pitch rate
figure;
margin(Cpid*sys1)
figure;
step(CL_sys)
figure;
step(CL_sys2)
figure;
pzmap(CL_sys)

%%
%---------------F---------------
CL_sys_delay1 = feedback(Cpid*sys1*0.05,sys1);
CL_sys_delay2 = feedback(Cpid*sys1*0.1,sys1);
CL_sys_delay3 = feedback(Cpid*sys1*0.2,sys1);
CL_sys_delay4 = feedback(Cpid*sys1*0.3,sys1);
CL_sys_delay5 = feedback(Cpid*sys1*0.5,sys1);
figure;
step(CL_sys_delay1)
figure;
step(CL_sys_delay2)
figure;
step(CL_sys_delay3)
figure;
step(CL_sys_delay4)
figure;
step(CL_sys_delay5)

%%
%---------------G---------------
t = 0:0.01:2;
figure
step(sys1, t)

Q = C_angle'*C_angle;
R = [0.1];

Co = ctrb(A,B);
unCo = length(Co)-rank(Co);
if unCo == 0
    fprintf('the system is completely controllable!\n')
else
    fprintf('the system has %d uncontrollable states!\n', unCo);
end


Ob = obsv(A,C_angle);
unOb = length(Ob)-rank(Ob);
if unOb == 0
    fprintf('the system is completely observable!\n')
else
    fprintf('the system has %d unobservable states!\n', unCo);
end

p = 50;
Q = p*C_angle'*C_angle;
R = 1;
[K] = lqr(A,B,Q,R)
Nbar = -inv(C_angle*inv(A-B*K)*B);
sys3 = ss(A-B*K,B*Nbar,C_angle,D);
step(0.1571*sys3)
ylabel('pitch angle (rad)');
title('CL Step Response (LQR)');
