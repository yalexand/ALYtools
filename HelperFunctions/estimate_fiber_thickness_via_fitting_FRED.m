function [sigIm,chi2Im] = estimate_fiber_thickness_via_fitting_FRED(sgm_filaments,mask,original,radius)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global M; % skel mask
global O; % original
global r; % radius of the vicinity
global xc; 
global yc; 
global X;
global Y;
global angle;

sigma_min   = 0.5;
sigma_max   = 16;

A_min = 1; 
A_max = Inf;

B_min = 0;
B_max  = Inf;

lb = [sigma_min A_min B_min ];
ub = [sigma_max A_max B_max ];

options = optimset('Jacobian','off','Display','off');
options.MaxFunEvals = options.MaxFunEvals*10;


O = original;
r = (-radius:radius);
[Y,X] = meshgrid(r,r);

M = zeros(size(Y));
M(radius+1,radius+1)=1;
M = bwdist(M);
M = M<=radius;

[sX,sY] = size(sgm_filaments);
sigIm=zeros(sX, sY);
chi2Im=zeros(sX, sY);

sklbl = bwlabel(sgm_filaments); % labelled skeleton

N = max(sklbl(:)); % number of fibers

stats = regionprops(sklbl,'PixelList');

hw = waitbar(0,'fittings.. please wait');
for k=1:N    
    list = stats(k).PixelList;
     y = list(:,1);
     x = list(:,2);

        for p=1:length(x)        
            if x(p)-radius >=1 && y(p)-radius >=1 && x(p)+radius <=sX && y(p)+radius <=sY
                img = sgm_filaments(x(p)-radius:x(p)+radius,y(p)-radius:y(p)+radius);
                msk = mask(x(p)-radius:x(p)+radius,y(p)-radius:y(p)+radius);
                img(~msk)=0;
                stats2 = regionprops(img,'Orientation');
                
                xc = x(p);
                yc = y(p);
                angle = stats2(1).Orientation*pi/180; % cos reported in degrees;
                
                    vic = O(xc+r,yc+r);
                    
                    % B_guess = min(vic(:));
                    B_guess = 0;
                                        
                    A_guess = max(vic(:)) - B_guess;
                    
                    sigma_guess = 3;
                    %
                    guess = [sigma_guess A_guess B_guess ];
                    %
                    [res,errval] = lsqnonlin('err_func_filthick',guess,lb,ub,options);
                    
%                    icyvol = show_filthick(res);                                        
%                    if errval < ...
%                       icy_imshow(icyvol);
%                    end
                    
                    sigIm(x(p),y(p)) = res(1);
                    chi2Im(x(p),y(p))= errval;                    
                                
                    %sigIm(x(p),y(p)) = stats2(1).Orientation;
                                                
            end
        end
if ~isempty(hw), waitbar(k/N,hw); drawnow, end;
end
if ~isempty(hw), delete(hw), drawnow; end;

end