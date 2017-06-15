 function dx = statespacesolution(t,x)
global Ixx Iyy Izz  Jr m l g kt kq thetafix k2
global R omega1 omega2 omega3 omega4 omegar omega12 omega22 omega32 omega42
%t

dx = zeros(12,1);
a1=(Iyy-Izz)/Ixx;
a2=Jr/Ixx;
a3=(Izz-Ixx)/Iyy;
a4=Jr/Iyy;
a5=(Ixx-Iyy)/Izz;
b1=l/Ixx;
b2=(l)/Iyy;
b3=1/Izz;


%{
 Gyrospopic Moments (using body angular rates)
t_phigyro=Jr*(3.14/30)*x_b(4)*(omega1-omega2+omega3-omega4); 
t_thetagyro=Jr*(3.14/30)*x_b(2)*(-omega1+omega2-omega3-omega4); 
%}


T1=k2*omega1^(2)*((thetafix/3) -x(13)/(2*omega1*R));
T2=k2*omega2^(2)*((thetafix/3) -x(14)/(2*omega2*R));
T3=k2*omega3^(2)*((thetafix/3) -x(15)/(2*omega3*R));
T4=k2*omega4^(2)*((thetafix/3) -x(16)/(2*omega4*R));




% State Equations

dx(1)=x(2);
dx(2)=(x(4)*x(6)*a1)-(x(4)*omegar*a2)+(b1*(T1-T3));
dx(3)=x(4);
dx(4)=(x(2)*x(6)*a3)+(x(2)*omegar*a4)+(b2*(T4-T2));
dx(5)=x(6);
dx(6)=x(2)*x(4)*a5 + b3*kq*(T1-T2+T3-T4);
dx(7)=x(8);
dx(8)=(-(T1+T2+T3+T4)/m)*(cos(x(1))*sin(x(3))*cos(x(5))+sin(x(1))*sin(x(5)));
dx(9)=x(10);
dx(10)=(-(T1+T2+T3+T4)/m)*(cos(x(1))*sin(x(5))*sin(x(3))-cos(x(5))*sin(x(1)));
dx(11)=x(12);
dx(12)=(g)-((T1+T2+T3+T4)/m)*cos(x(1))*cos(x(3));
dx(13)=-5.22*x(13);
dx(14)=-5.22*x(14);
dx(15)=-5.22*x(15);
dx(16)=-5.22*x(16);
end

