function [N,avr_val,std_val] = estimate_fiber_thickness_via_fitting(sgm_filaments,sgm_junctions,original)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global M; % skel mask
global O; % original
global r;
global xc; 
global yc; 
global X;
global Y;

f = 2;

loci = imresize(sgm_filaments & sgm_junctions,f);
u = imresize(original,f,'bicubic');
loci = bwmorph(loci,'thin',Inf);

%icyvol(:,:,1,1,1) = u;
%icyvol(:,:,2,1,1) = loci;
%icy_imshow(icyvol);

O = u;

vR = ceil(3*f);
r = (-vR:vR);
[Y,X] = meshgrid(r,r);

M = imdilate(loci,strel('disk',vR));

t_guess     = pi/2;
t_min       = 0; 
t_max       = 2*pi;

sigma_guess = 2;
sigma_min   = 0.5;
sigma_max   = 16;

A_guess = 1;
A_min = 0; 
A_max = Inf;

B_guess = 1;
B_min = 0;
B_max  = Inf;

lb = [t_min sigma_min A_min B_min ];
ub = [t_max sigma_max A_max B_max ];

guess = [t_guess sigma_guess A_guess B_guess ];
options = optimset('Jacobian','off','Display','off');
options.MaxFunEvals = 6000*10;                                
options.MaxIter = 900000;
 
data = [];
cnt = 0;
hw = waitbar(0,'fittings.. please wait');
[sX,sY] = size(loci);
for x=1:sX,
    if ~isempty(hw), waitbar(x/sX,hw); drawnow, end;
    for y=1:sY,
            if 0~=loci(x,y) && x > vR+1 && x < sX-vR-1 && y > vR+1 && y < sY-vR-1 
                cnt = cnt + 1;
                xc = x;
                yc = y;
                [res,resnorm] = lsqnonlin('err_func_filthick',guess,lb,ub,options);                
                icyvol = show_filthick(res);
                theorvals = icyvol(:,:,2,1,1);
                theorvals = theorvals(theorvals > 0);
                if max(theorvals(:)) - min(theorvals(:)) > 3                
                    rec = {resnorm res(1) res(2) res(3) res(4)};
                    data = [data; rec];
                    % disp(num2str([resnorm res(1) res(2) res(3) res(4)]));                
                end
            end
    end
end
if ~isempty(hw), delete(hw), drawnow; end;

err1 = cell2mat(data(:,1));
r1 = cell2mat(data(:,3));
s1 = r1(intersect(find( r1<3 ),find( err1<5 )));
N = numel(s1);
avr_val = mean(s1)/f;
std_val = std(s1)/f;

end

