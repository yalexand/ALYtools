function z = box_average_3D(U,d) % diameter

z = [];

if ~3==numel(size(U)), return, end;

r = round(d);
mask = ones(r,r,r)/r^3;

% CONVN SLOW
% z = convn(U,mask,'same'); 
% norm_image = convn(ones(size(U)),mask,'same'); 
% z = z./norm_image;

% USE CONVOLUTION IN FREQUENCY DOMAIN INSTEAD
z = conv3Dfreq(U,mask); 

