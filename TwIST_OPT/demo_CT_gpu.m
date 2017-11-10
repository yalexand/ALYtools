function demo_CT
%% EDIT BY TC
% CT demo: Reconstruction from noisy random projections
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

n =1300;

beams = 45;
angle = 180/beams;
angles =[0:angle:180-angle];

% Measurement noise standard deviation
sigma = 1e-1;

% generate phantom
f = gpuArray(phantom(n));

colormap(gray);
figure(1)
subplot(221);imagesc(f);
title('Original image')


hR = @(x)  radon(x, angles);
hRT = @(x) iradon(x, angles,'linear','ram-lak',1,n);

% Compute Radon transform 
R = hR(f);
R = R(:,:,1);
subplot(222);imagesc(R);
title('Radon transform')

% vector of observations
R_noise = gpuArray(R); %+ sigma*randn(size(R));
subplot(223);imagesc(R_noise);
title('Radon transform with noise')

% compute the  filtered backprojection image
[FB,H] = hRT(R_noise);
figure(1); colormap gray; 
subplot(224);imagesc(FB); axis off;
title('Filtered backprojection image - Hann filter')

y=R_noise;
% denoising function;
tv_iters = gpuArray(5);
Psi = @(x,th)  tvdenoise(x,2/th,tv_iters);

% set the penalty function, to compute the objective
Phi = @(x) TVnorm_gpu(x);

% regularization parameters (empirical)
tau = gpuArray(0.08);

tolA = 1e-6;
% -- TwIST ---------------------------
% stop criterium:  the relative change in the objective function 
% falls below 'ToleranceA'         
 [x_twist,dummy,obj_twist,...
    times_twist,dummy,mse_twist]= ...
         TwIST_gpu(y,hR,tau,...
         'Lambda', 1e-4, ...
         'AT', hRT, ...
         'Psi', Psi, ...
         'Phi',Phi, ...
         'True_x', f, ...
         'Monotone',1,...
         'MaxiterA', 10000, ...
         'Initialization',0,...
         'StopCriterion',1,...
       	 'ToleranceA',tolA,...
         'Verbose', 1);
   
% [x_twist,dummy,obj_twist,...
%     times_twist,dummy,mse_twist]= ...
%          SpaRSA(y,hR,tau,...
%          'AT', hRT, ...
%          'Psi', Psi, ...
%          'Phi',Phi, ...
%          'True_x', f, ...
%          'BB_factor', 0.8, ...
%          'Monotone',1,...
%          'Initialization',0,...
%          'StopCriterion',4,...
%        	 'ToleranceA',obj_twist(end),...
%          'Verbose', 1);


figure(11);
%set(11,'Position',[400 385 300 300])
subplot(131);imagesc(x_twist);
axis image
title('Estimate')
colormap(gray)


%figure(12);
%set(12,'Position',[800 0 300 300])
subplot(132);semilogy(times_twist, (mse_twist*prod(size(f))).^0.5/norm(f,'fro'),'LineWidth',2);axis square
title('error ||x^{t}-x||_2/||x||_2')
xlabel('CPU time (sec)')
grid on

%figure(13)
%set(13,'Position',[800 400 300 300])
subplot(133);semilogy(times_twist,obj_twist,'b','LineWidth',2); axis square
title('Objective function')
xlabel('CPU time (sec)')
grid on

