function gtxt = filaments_detection(u,sigma1,sigma2,nmax,n)

    angles = ((1:n)-1)*pi/n;
        
%     figure(2);
%     for k = 1:16
%         subplot(4,4,k);
%         dk = compute_directional_kernel(sigma1,sigma2,angles(k),nmax);
%         k
%         imagesc(dk); daspect([1 1 1]);
%     end
    
    gtxt = [];
    for k = 1:n
        dk = compute_directional_kernel(sigma1,sigma2,angles(k),nmax);
        res = abs(conv2(u,dk,'same'));
        if isempty(gtxt)
            gtxt = res;
        else
            gtxt = pixelwise_max(gtxt,res);
        end;    
    end
    
end