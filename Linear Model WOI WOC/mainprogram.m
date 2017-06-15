 %linearised decoupled model 
%using ode45 to solve
%load function sys.m before executing

clear
  global A B  C D
    
close all;
clc;

syms dtheta1 dtheta2 dtheta3 dtheta4 phi phi_dot
syms theta theta_dot
syms    psi psi_dot x x_dot y y_dot z z_dot 
syms omega1 omega2 omega3 omega4 

syms Ix Iy Iz g m kt2 kq a1 a2 a3 a4 a5 b1 b2 b3
syms thetafix vi1 vi2 vi3 vi4 k2 R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ix=4.856*10^(-3);
Iy=4.856*10^(-3);
Iz=8.801*10^(-3);
l=0.225; %% arm length 


Jr=3.357*10^(-6); %% Rotor Inertia
a1=(Iy-Iz)/Ix;
a2=Jr/Ix;
a3=(Iz-Ix)/Iy;
a4=Jr/Iy;
a5=(Ix-Iy)/Iz;
b1=l/Ix;
b2=(l)/Iy;
b3=1/Iz;


%Mass
m=.468;
g=9.81;
%Aerodynamic force and Moment constant
kt=2.98*10^(-6);
kq=1.14*10^(-7);
omega=sqrt(m*g/(4*kt));
%{

omega4 = omega, omega3 = omega4, omega2 = omega3, omega1 = omega2 ;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state vectors
u = [omega1; omega2;omega3;omega4];
x = [phi;phi_dot;theta;theta_dot;
    psi;psi_dot;x;x_dot;y;y_dot;z;z_dot;]; 
y = [phi;phi_dot;theta;theta_dot;
    psi;psi_dot;x;x_dot;y;y_dot;z;z_dot];



%ct=3.13*10^(-5)
T1=kt*omega1^(2);
T2=kt*omega2^(2);
T3=kt*omega3^(2);
T4=kt*omega4^(2);





% non-linear system
F1=phi_dot;
F2=(theta_dot*psi_dot*a1)-(theta_dot*(omega1-omega2+omega3-omega4)*a2)+(b1*(T1-T3));
F3=theta_dot;
F4=(phi_dot*psi_dot*a3)+(phi_dot*(omega1-omega2+omega3-omega4)*a4)+(b2*(T4-T2));
F5=psi_dot;
F6=phi_dot*theta_dot*a5 + b3*kq*(T1-T2+T3-T4);
F7=x_dot;
F8=(-(T1+T2+T3+T4)/m)*(cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi));
F9=y_dot;
F10=(-(T1+T2+T3+T4)/m)*(cos(phi)*sin(psi)*sin(theta)-cos(psi)*sin(phi));
F11=z_dot;
F12=(g)-( (T1+T2+T3+T4)/m )*cos(phi)*cos(theta);



F=[F1;F2;F3;F4;F5;F6;F7;F8;F9;F10;F11;F12];

A.symbolic = jacobian(F, x);
B.symbolic = jacobian(F, u);
A.symbolic=subs(A.symbolic,sym('phi'),0);
A.symbolic=subs(A.symbolic,sym('phi_dot'),0);
A.symbolic=subs(A.symbolic,sym('theta'),0);
A.symbolic=subs(A.symbolic,sym('theta_dot'),0);
A.symbolic=subs(A.symbolic,sym('psi'),0);
A.symbolic=subs(A.symbolic,sym('psi_dot'),0);
A.symbolic=subs(A.symbolic,sym('x'),0);
A.symbolic=subs(A.symbolic,sym('x_dot'),0);
A.symbolic=subs(A.symbolic,sym('y'),0);
A.symbolic=subs(A.symbolic,sym('y_dot'),0);
A.symbolic=subs(A.symbolic,sym('z'),0);
A.symbolic=subs(A.symbolic,sym('z_dot'),0);
A.symbolic=subs(A.symbolic,sym('R'),.125);
A.symbolic=subs(A.symbolic,sym('omega1'),omega);
A.symbolic=subs(A.symbolic,sym('omega2'),omega);
A.symbolic=subs(A.symbolic,sym('omega3'),omega);
A.symbolic=subs(A.symbolic,sym('omega4'),omega);

A=eval(A.symbolic);


B.symbolic=subs(B.symbolic,sym('phi'),0);
B.symbolic=subs(B.symbolic,sym('phi_dot'),0);
B.symbolic=subs(B.symbolic,sym('theta'),0);
B.symbolic=subs(B.symbolic,sym('theta_dot'),0);
B.symbolic=subs(B.symbolic,sym('psi'),0);
B.symbolic=subs(B.symbolic,sym('psi_dot'),0);
B.symbolic=subs(B.symbolic,sym('x'),0);
B.symbolic=subs(B.symbolic,sym('x_dot'),0);
B.symbolic=subs(B.symbolic,sym('y'),0);
B.symbolic=subs(B.symbolic,sym('y_dot'),0);
B.symbolic=subs(B.symbolic,sym('z'),0);
B.symbolic=subs(B.symbolic,sym('z_dot'),0);
B.symbolic=subs(B.symbolic,sym('R'),.125);
B.symbolic=subs(B.symbolic,sym('omega1'),omega);
B.symbolic=subs(B.symbolic,sym('omega2'),omega);
B.symbolic=subs(B.symbolic,sym('omega3'),omega);
B.symbolic=subs(B.symbolic,sym('omega4'),omega);

B=eval(B.symbolic);
C=eye(12);
D=0*rand(12,4);

tspan = [0 30];
iniCon = [0.25 .05 0.35 0.06 0 0 0 0 0 0 0 0];
[t, x] = ode45(@func, tspan, iniCon);


%Plotting
figure
subplot(2,3,1);
plot(t,x(:,1:2));
legend('phi','phidot');
xlabel('Time(sec)');
ylabel('Amplitude');
subplot(2,3,2);
plot(t,x(:,3:4));
legend('theta','thetadot');
xlabel('Time(sec)');
ylabel('Amplitude');

subplot(2,3,3);
plot(t,x(:,5:6));
legend('psi','psidot');
xlabel('Time(sec)');
ylabel('Amplitude');

subplot(2,3,4);
plot(t,x(:,7:8));
legend('x','x dot');
xlabel('Time(sec)');
ylabel('Amplitude');
subplot(2,3,5);
plot(t,x(:,9:10));
legend('y','y dot');
xlabel('Time(sec)');
ylabel('Amplitude');

subplot(2,3,6);
plot(t,x(:,11:12));
legend('z','z dot');
xlabel('Time(sec)');
ylabel('Amplitude');

