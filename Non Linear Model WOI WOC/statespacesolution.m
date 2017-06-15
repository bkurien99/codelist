function dx = statespacesolution(t,x)
global Ixx Iyy Izz  Jr m l g kt kq 
global omega1 omega2 omega3 omega4 omegar omega12 omega22 omega32 omega42
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


 %Thrust Matrix
u1=kt*(omega12 +omega22 +omega32+omega42);
u2=kt*(-omega22+omega42); %rolling left --> positive
u3=kt*(omega12-omega32); % pitching up --> positive
u4=kq*(omega12-omega22+omega32-omega42);





% State Equations

dx(1)=x(2);
dx(2)=(x(4)*x(6)*a1)-(x(4)*omegar*a2)+(b1*u2);
dx(3)=x(4);
dx(4)=(x(2)*x(6)*a3)+(x(2)*omegar*a4)+(b2*u3);
dx(5)=x(6);
dx(6)=x(2)*x(4)*a5 + b3*u4;
dx(7)=x(8);
dx(8)=(-u1/m)*(cos(x(1))*sin(x(3))*cos(x(5))+sin(x(1))*sin(x(5)));
dx(9)=x(10);
dx(10)=(-u1/m)*( (cos(x(1))*sin(x(5))*sin(x(3)) )-(cos(x(5))*sin(x(1))) );
dx(11)=x(12);
dx(12)=(g)-(u1/m)*cos(x(1))*cos(x(3));




end

