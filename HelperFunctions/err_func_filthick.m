function [res] = err_func_filthick( x )

global M; % skel mask
global O; % original
global r; % radius of the vicinity
global xc; 
global yc; 
global X;
global Y;
global angle;

t       = angle;
sigma   = x(1);
A       = x(2);
B       = x(3);
    
  X1 = ( cos(t)*X + sin(t)*Y ) / (sqrt(2)*sigma);
  f = exp( - X1.^2 );
  
  f = A*f/sum(f(:)) + B;        
  % f = A*f + B; % don't forget show_fitthick!       
  
%diff = (O(xc+r,yc+r) - f);
diff = (O(xc+r,yc+r) - f).*M;

res = diff(:)/numel(diff(:));


