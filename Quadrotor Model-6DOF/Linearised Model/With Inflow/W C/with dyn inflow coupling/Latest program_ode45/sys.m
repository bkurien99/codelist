function  dx = sys(t, x)
      
global A B K C D
    
u=-K*x ;
dx=A*x+B*u ;

       
       
end

