
clear ; close all
clc;
clf;
global Ixx Iyy Izz  Jr m l g kt kq  omega1 omega2 omega3 omega4 omegar
global omega12 omega22 omega32 omega42 
%Inertia matrix Components Initialization
Ixx=4.856*10^(-3);
Iyy=4.856*10^(-3);
Izz=8.801*10^(-3);
J=[Ixx 0 0;0 Iyy 0;0 0 Izz];

%Body Dimensions
l=0.225;
%%Rotor inertia
Jr=3.357*10^(-6);
%Mass
m=.468;
g=9.81;
%Aerodynamic force and Moment constant
kt=2.98*10^(-6);
kq=1.14*10^(-7);

%Rotor speed
omegas=(m*g/(4*kt));
omegas=1*omegas ;
omega12=1*omegas ;
omega22=omegas ;
omega32=1*omegas ;
omega42=omegas ;

omega1=sqrt(omega12);
omega2=sqrt(omega22);
omega3=sqrt(omega32);
omega4=sqrt(omega42);

omegar=-omega1+omega2-omega3+omega4;

%Initial Condition
xi_b= [0.25 .05 0.35 0.06 0 0 0 0 0 0 0 0];

%Funktion Call 
t = 0:.01:30;
[t,xa] = ode45(@statespacesolution,t,xi_b);
 
%Plotting
subplot(2,3,1);
plot(t,xa(:,1:2));
legend('phi','phidot');
xlabel('Time(sec)');
ylabel('Amplitude');
subplot(2,3,2);
plot(t,xa(:,3:4));
legend('theta','thetadot');
xlabel('Time(sec)');
ylabel('Amplitude');

subplot(2,3,3);
plot(t,xa(:,5:6));
legend('psi','psidot');
xlabel('Time(sec)');
ylabel('Amplitude');

subplot(2,3,4);
plot(t,xa(:,7:8));
legend('x','x dot');
xlabel('Time(sec)');
ylabel('Amplitude');
subplot(2,3,5);
plot(t,xa(:,9:10));
legend('y','y dot');
xlabel('Time(sec)');
ylabel('Amplitude');

subplot(2,3,6);
plot(t,xa(:,11:12));
legend('z','z dot');
xlabel('Time(sec)');
ylabel('Amplitude');


%{
subplot(1,2,1);
plot(t,xa(:,1:6));
legend('phi','phidot','theta','thetadot','psi','psidot');
xlabel('Time(sec)');
ylabel('Amplitude');

subplot(1,2,2);
plot(t,xa(:,7:12));
legend('x','xdot','y','ydot','z','zdot');
xlabel('Time(sec)');
ylabel('Amplitude');
%}
