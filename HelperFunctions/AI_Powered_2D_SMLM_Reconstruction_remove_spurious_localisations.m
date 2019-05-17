%-------------------------------------------------------------------------%     
function [XY_corr,F_corr] = AI_Powered_2D_SMLM_Reconstruction_remove_spurious_localisations(sx,sy,XY,F,R)

    r = max(1,fix(R));

    x = XY(:,1);
    y = XY(:,2);
    
    N = numel(x);
    
    se=strel('disk',r,0);
    mask = double(se.Neighborhood);
    %
    D = 2*r+2;
    uext = zeros(sx+2*D,sy+2*D);
    shift = D;
    for k=1:N
        uext(shift+x(k),shift+y(k)) = 1;
    end

    z=conv2(uext,mask,'same');
    
%     iv=zeros(sx+2*D,sy+2*D,2,1,1,1);
%     iv(:,:,1,1,1) = uext;
%     iv(:,:,2,1,1) = z;
%     icy_imshow(iv);

    single_vic_sum = sum(mask(:));
  
%     tic
%     vics = zeros(size(mask,1),size(mask,2),N);
%     parfor k=1:N
%         xc = shift + x(k);
%         yc = shift + y(k);
%         z_=z;
%         vics(:,:,k) = z_(xc-r:xc+r,yc-r:yc+r);
%     end
%     vics = vics.*repmat(mask,[1 1 N]);
%     vic_sums = squeeze(sum(vics,[1 2]));
%     non_separate_pixel = vic_sums~=single_vic_sum.*ones(N,1);
%     all_inds=1:N;
%     non_separate_pixel = all_inds(0~=non_separate_pixel);
%     toc

    non_separate_pixel = 1:N;  % indices    
    for k=1:N
        xc = shift + x(k);
        yc = shift + y(k);
        vic = z(xc-r:xc+r,yc-r:yc+r).*mask;
        sum_cur_vic = sum(vic(:));
        if sum_cur_vic==single_vic_sum
            non_separate_pixel(k) = 0;
        end
    end
    non_separate_pixel = non_separate_pixel(0~=non_separate_pixel);
    
    x_=x(non_separate_pixel);
    y_=y(non_separate_pixel);
    
    F_corr = F(non_separate_pixel);
    
%     u_ = zeros(sx,sy);
%     for k=1:length(x_)
%         u_(x_(k),y_(k))=1;
%     end
%     z=conv2(u_,mask,'same');
%     iv=zeros(sx,sy,2,1,1,1);
%     iv(:,:,1,1,1)=u_;
%     iv(:,:,2,1,1)=z;
%     icy_imshow(iv);
    
    XY_corr = [x_ y_];
end