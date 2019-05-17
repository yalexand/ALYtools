%-------------------------------------------------------------------------%     
function I = AI_Powered_2D_SMLM_Reconstruction_ASH_2d(sX,sY,XYs_in,d)
% Averaged Shifted Histograms in 2d (image), depending on window size d
% attention - size and coordinates should be consistent, scale-wise!

    I = zeros(sX,sY);

    suitable_for_plotting = XYs_in(:,1) >= ceil(d/2) & XYs_in(:,2) >= ceil(d/2) & ... 
                            XYs_in(:,1) <= sX-ceil(d/2) & XYs_in(:,2) <= sY-ceil(d/2);

    XYs = round(XYs_in(suitable_for_plotting,:));

    w = pyramid_weights(d);
    r = max(1,floor(d/2)); % vicinity radius

    for k = 1:size(XYs,1)
        x = XYs(k,1);
        y = XYs(k,2);
        I(x-r:x+r,y-r:y+r) = I(x-r:x+r,y-r:y+r) + w;
    end

end
%-------------------------------------------------------------------------%     
function w = pyramid_weights(d)
% pyramid_weights(d) - creates special normalized wiights mask of the size d
%   used in Averaged Shifted Histograms visualization of
%   emitters localisations in SMLM

    r = max(1,floor(d/2)); % vicinity radius

    w = zeros(2*r+1,2*r+1);
    
    xc = r+1;
    yc = r+1;

    for x=-r:r
        for y=-r:r
            w(xc+x,yc+y) = (r - abs(x)+1)*(r - abs(y)+1);
        end
    end

    w = w/sum(w(:));
    
end