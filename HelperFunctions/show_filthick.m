function icyvol = show_filthick( x )

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
    f = exp( -X1.^2 );
    f = A*f/sum(f(:)) + B;
        
icyvol(:,:,1,1,1) = O(xc+r,yc+r);
icyvol(:,:,2,1,1) = f;
% icyvol(:,:,3,1,1) = double(M(xc+r,yc+r));

end

