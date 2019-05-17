%-------------------------------------------------------------------------%     
function DX_DY_DRIFT = AI_Powered_2D_SMLM_Reconstruction_calculate_XY_drift_increments(SX,SY,XY,F,n_frames,upscale_fac)

    DX_DY_DRIFT = zeros(size(F,1),2);

    frames = 1:n_frames;    
    
    nf = zeros(n_frames,1);
    for f=1:n_frames
        ind = F(F==f);
        ind = ind(ind~=0);
        nf(f) = numel(ind);
    end        
%     figure;
%     plot(frames,nf,'b.-');
%     grid on;

    z = cumsum(nf)/sum(nf(:));
    
    q = [0 0.2 0.4 0.6 0.8 1];
    thresh = zeros(1,numel(q)-1);
    for k=1:numel(q)-1
        selection = frames(q(k)<z&z<q(k+1));
        selection = selection(0~=selection);
        thresh(k) = max(selection(:));
    end
%     figure;
%     plot(1:numel(z),z,'r.-',thresh,zeros(size(thresh)),'bo-');
%     grid on;
    
    thresh = [1 thresh];    
          
    img = zeros(SX,SY,numel(q)-1,1,1);
    for k=1:numel(q)-1
        range_k = thresh(k):thresh(k+1);
        frames_k = ismember(F,range_k);
        XY_k = XY(frames_k,:);
        % size(XY_k,1)
        img(:,:,k,1,1) = AI_Powered_2D_SMLM_Reconstruction_ASH_2d(SX,SY,XY_k,max(3,round(upscale_fac/2)));                
    end
    %icy_imshow(img);
%     
%     
        u1 = squeeze(img(:,:,1,1,1)); % first image
        z1 = u1>0;
        z1 = bwareaopen(z1,upscale_fac*upscale_fac);
        u1 = u1.*z1;
        
        dx = 0;
        dy = 0;
        
        %parfor k=2:size(img,3) % possible but causeed memory faults on large FOVs
        for k=2:size(img,3)
            u2 = squeeze(img(:,:,k,1,1));    
            z2 = u2>0;
            z2 = bwareaopen(z2,upscale_fac*upscale_fac);
            u2 = u2.*z2;
            %
            z = xcorr2_fft(u1,u2); 
            %
            s1=3;
            s2=7;
            g1 = gsderiv(z,s1,0);
            g2 = gsderiv(z,s2,0);
            z = (g1-g2);
            %
            z(z<0)=0;
            %
            [wc,hc] = size(z);
            wc=fix(wc/2);
            hc=fix(hc/2);
            %
            [x,y] = find(z==max(z(:)));
            d_shift = [x-wc y-hc] - 1; % !!!!!!!!
            dx = [dx; d_shift(1)];
            dy = [dy; d_shift(2)];
%            edges_q_k = z(x,y)/mean(z(:)); % quality                                                    
%            icy_imshow(z,num2str(k));
        end
                
          f = [];
          for k=2:length(thresh)
              range_k=thresh(k-1):thresh(k);
              f_k = sum(range_k'.*nf(range_k))/sum(nf(range_k));
              f = [f; f_k];
          end
             
          degree = 2;
          [px,Sx,mux] = polyfit(f,dx,degree);
          [ix,~]=polyval(px,1,Sx,mux);
          [dx_fitted,~]=polyval(px,frames,Sx,mux);
          res_dx = polyval(px,frames,Sx,mux)-ix;
          %
          [py,Sy,muy] = polyfit(f,dy,degree);
          [iy,~]=polyval(py,1,Sy,muy);
          [dy_fitted,~]=polyval(py,frames,Sy,muy);
          res_dy = polyval(py,frames,Sy,muy)-iy;
          
%           figure
%           plot(f,dx-ix,'bo',f,dy-iy,'ro',frames,dx_fitted-ix,'b.-',frames,dy_fitted-iy,'r.-');
                    
          DX_DY_DRIFT = [res_dx' res_dy'];
end