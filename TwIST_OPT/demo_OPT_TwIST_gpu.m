function demo_OPT_TwIST
%% EDIT BY TC
% OPT demo: Reconstruction from noisy random projections
%
% This demo illustrates the twIST algorithm computing the solution  of  
%
%     xe = arg min 0.5*||A x-y||^2 + tau TV(x),
%             x
%  where x is an image, Ax is the Fourier Transform of x  computed in a 
%  "small" set of frequencies, y is a noisy observation and TV(x) is
%  the total variation of x.
% 
% The proximal operator, i.e., the solution of
%
%     xe = arg min 0.5 ||x-y||^2 + tau TV(x)
%             x
% is  given by the Chambolle algorithm
% 
% A. Chambolle, �An algorithm for total variation minimization and
% applications,� Journal of Mathematical Imaging and Vision, vol. 20,
% pp. 89-97, 2004.
%
%
% For further details about the TwIST algorithm, see the paper:
%
% J. Bioucas-Dias and M. Figueiredo, "A New TwIST: Two-Step
% Iterative Shrinkage/Thresholding Algorithms for Image 
% Restoration",  IEEE Transactions on Image processing, 2007.
%
%
% (available at   http://www.lx.it.pt/~bioucas/publications.html)
%
%
% Authors: Jose Bioucas-Dias and Mario Figueiredo, 
% Instituto Superior T�cnico, October, 2007
%
%
addpath('/cs/research/medim/projects2/projects/tcorreia/Matlab/IST/TwIST_v2')

dangle = 8;
angles = 0:dangle:(180-dangle);
beams = length(angles);
recon =zeros(430,430,430);
imangles = 1:dangle:180;

% load data
load fluor
im= images(16:end-17,:,:);
N=size(im,1);

% The number of rows of the sinogram given by radon is
% not the size of the image and we need to zero pad the images.
zpt = ceil((2*ceil(norm(size(im(:,:,1))-floor((size(im(:,:,1))-1)/2)-1))+3 - N)/2);
zpb = floor((2*ceil(norm(size(im(:,:,1))-floor((size(im(:,:,1))-1)/2)-1))+3 - N)/2);
st = abs(zpb - zpt); 

for j = 120 %size(images,2)
    count = 0;
    clear R
    for i = 1:beams
        count = count + 1;
        Rpad = padarray(squeeze(im(:,j,imangles(i))),zpt,0,'pre');
        R(:,count) = padarray(Rpad,zpb,0,'post');
    end
    
hR = @(x)  radon(x, angles);
hRT = @(x) iradon(x, angles,'linear','Ram-Lak',1,N);

% compute the  filtered backprojection image
[FB,H] = hRT(R);
figure(1); colormap jet; 
imagesc(FB); axis off;
title('Filtered backprojection image - Hann filter')

y=gpuArray(R./max(R(:)));
% denoising function;
tv_iters = 25;
Psi = @(x,th)  tvdenoise(x,2/th,tv_iters);

% set the penalty function, to compute the objective
Phi = @(x) TVnorm_gpu(x);

% regularization parameters (empirical)
tau = gpuArray(0.0008);

tolA = 1e-4;
% -- TwIST ---------------------------
% stop criterium:  the relative change in the objective function 
% falls below 'ToleranceA'         
 [x_twist,dummy,obj_twist,...
    times_twist,dummy,mse_twist]= ...
         TwIST_gpu_OPT(y,hR,tau,...
         'Lambda', 1e-4, ...
         'AT', hRT, ...
         'Psi', Psi, ...
         'Phi',Phi, ...
         'Monotone',1,...
         'MaxiterA', 10000, ...
         'Initialization',0,...
         'StopCriterion',1,...
       	 'ToleranceA',tolA,...
         'Verbose', 1);
   


figure(11);
set(11,'Position',[400 385 300 300])
imagesc(x_twist);
axis image
title('Estimate')


recon(:,:,j) = gather(x_twist);
if (j==430 )
save recon_TwIST_22 recon
end
end

