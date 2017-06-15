function  dx = sys(t, x)
      
global A B k C D
    
u=-k*x ;
[dx]=A*x+B*u ;

       
       
end

