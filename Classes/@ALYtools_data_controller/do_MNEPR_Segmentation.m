        function sgm = do_MNEPR_Segmentation(obj,send_to_Icy,~) 
   
            u = single(obj.imgdata);
                    
            puncta = zeros(size(u)); 
            nuclei = zeros(size(u)); 
            %
            [sX,sY,sC,sZ,sT] = size(u);
            %
            for tind = 1:sT
                img = squeeze(u(:,:,1,1,tind));
                p = segment_puncta(img);
                puncta(:,:,1,1,tind) = p; 
                n = segment_nuclei(img,p);
                nuclei(:,:,1,1,tind) = n; 
            end
            %
        % nuclear correction - starts
            ref = zeros(sX,sY); % sum of frames
            for tind = 1:sT  
                n = squeeze(nuclei(:,:,1,1,tind));
                ref = ref + n;
            end
            %
            refmax = zeros(size(ref));
            L = bwlabel(ref);
            for l=1:max(L(:))
                img = ref(L==l);
                refmax(L==l) = max(img(:));
            end                        
            %
            for tind = 1:sT
                n = squeeze(nuclei(:,:,1,1,tind));
                L = bwlabel(n);
                for l=1:max(L(:))
                    % condition on max
                    img = refmax(L==l);
                    m = mean(img(:));
                    if m < sT-6 %remove this object
                        n(L==l)=0;
                    end
                        %condition on average
                        t = 0.2;
                        img = ref(L==l);
                        m = mean(img(:));
                        if m < t*sT %remove this object
                            n(L==l)=0;
                        end                    
                    nuclei(:,:,1,1,tind) = n;
                end
            end                                   
            % nuclear correction - ends
            %
            sgm = zeros(sX,sY,2,sZ,sT);            
            sgm(:,:,1,1,:) = puncta;
            sgm(:,:,2,1,:) = nuclei;
            %                                    
                if send_to_Icy                
                    try
                        icyvol = zeros(sX,sY,3,sZ,sT,'uint16');
                        %
                        icyvol(:,:,1,1,:) = uint16(u);
                        icyvol(:,:,2,1,:) = uint16(puncta);
                        icyvol(:,:,3,1,:) = uint16(nuclei);
                        %
                        notification = [obj.current_filename ' - Segmentation: MNEPR_Segmentation'];
                        if isempty(obj.h_Icy_segmentation_adjustment)
                            obj.h_Icy_segmentation_adjustment = icy_imshow(icyvol,notification);                    
                        else
                            icy_imshow(obj.h_Icy_segmentation_adjustment,icyvol,notification);                    
                        end
                    catch
                        errordlg('problem with Icy, - might be not running');
                    end                
                end
        %
        end                
%-------------------------------------------------------------------------%        
        function p = segment_puncta(u)
        
                K = 3;
                S1 = 2;
                S2 = 4;
                a = 0.5;
                t = 0.1;
                 
                nth1 = nonlinear_tophat(u,S1,K)-1;
                nth1(nth1<0)=0;
                nth2 = nonlinear_tophat(u,S2,K)-1;
                nth2(nth2<0)=0;
                
                str1 = zeros(size(u));
                str2 = zeros(size(u));
                %
                sigma1 = fix(S1/2);
                sigma2 = fix(S2/2);
                [uxx1,uxy1,uyy1] = gsderiv(u,sigma1,2);
                [uxx2,uxy2,uyy2] = gsderiv(u,sigma2,2);
                %
                mode = 'Peak';
                %
                % dirty solution
                hw = waitbar(0,['Conducting segmentation, please wait']);
                for x=1:size(u,1)
                    if ~isempty(hw), waitbar(x/size(u,1),hw); drawnow, end;                            
                    for y=1:size(u,2)
                        H = [uxx1(x,y) uxy1(x,y); uxy1(x,y) uyy1(x,y)];
                        [~,D] = eig(H); 
                        upp = D(1,1); 
                        uqq = D(2,2);
                        %
                        switch mode 
                            case 'Ridge'
                                if upp < 0
                                    str1(x,y) = abs(upp);
                                end
                            case 'Peak'
                            if upp < 0 && uqq < 0
                                str1(x,y) = sqrt(upp*uqq);
                            end             
                        end
                        %
                        H = [uxx2(x,y) uxy2(x,y); uxy2(x,y) uyy2(x,y)];
                        [~,D] = eig(H); 
                        upp = D(1,1); 
                        uqq = D(2,2);
                        %
                        switch mode 
                            case 'Ridge'
                                if upp < 0
                                    str2(x,y) = abs(upp);
                                end
                            case 'Peak'
                            if upp < 0 && uqq < 0
                                str2(x,y) = sqrt(upp*uqq);
                            end             
                        end                                                         
                    end                    
                end
                if ~isempty(hw), delete(hw), drawnow; end;
                %                             
                z = a*nth1.*map((str1),0,1) + (1-a)*nth2.*map((str2),0,1);
                z = imclose(z>t,strel('disk',fix(max(1,min(S1,S2)/2)),0));
                %
                p = z;        
        end
%-------------------------------------------------------------------------%                
        function nucs = segment_nuclei(u,puncta)
        
                u = max(u(:)) - u; % want ditches                
 
    % https://stackoverflow.com/questions/46000633/nlinfit-appears-better-than-fitgmdist-for-fitting-normal-mixture     
    % Increase no. of iterations. default is 100.
%     Data = u(:);
%     [Yhist, x] = histcounts(Data, 100, 'normalization', 'pdf');
%     x = (x(1:end-1) + x(2:end))/2;    
%     opts.MaxIter = 300;
%     % Ensure that it does all iterations.
%     opts.TolFun = 0;
% 
%     GMModel = fitgmdist(Data, 3, 'Options', opts, 'Start', 'plus');
%     wts = GMModel.ComponentProportion;
%     mu = GMModel.mu;
%     sig = sqrt(squeeze(GMModel.Sigma));
%     Ygmfit = wts(1)*normpdf(x(:),mu(1),sig(1)) + ...
%              wts(2)*normpdf(x(:),mu(2),sig(2)) + ...
%              wts(3)*normpdf(x(:),mu(3),sig(3));             
%     figure;
%     hold on
%     plot(x(:), Yhist, 'b');
%     plot(x(:), Ygmfit, 'k');
                                    
    Data = u(:);
    opts.MaxIter = 300;
    % Ensure that it does all iterations.
    opts.TolFun = 0;
    %
    GMModel = fitgmdist(Data, 3, 'Options', opts, 'Start', 'plus');
    mu = GMModel.mu;
    sig = sqrt(squeeze(GMModel.Sigma));
    %
    min_index = find(mu==min(mu));
    max_index = find(mu==max(mu));
    inds = [1 2 3];
    middle_index = inds(inds~=max_index & inds~=min_index);
    %
    % this is "bright" mask
    mask = u > mu(max_index) - 0.5*sig(max_index);
 
    % set puncta to bright value before segmenting - can help a bit to
    % nuclear segmentation?
    u_ = u;
    u_(puncta) = mu(max_index) + sig(max_index);
    
                K = 2.5;
                S1 = 9;
                S2 = 24;
                a = 0.5;
                t = 0.05;
                 
                nth1 = nonlinear_tophat(u_,S1,K)-1;
                nth1(nth1<0)=0;
                nth2 = nonlinear_tophat(u_,S2,K)-1;
                nth2(nth2<0)=0;
                
                str1 = zeros(size(u));
                str2 = zeros(size(u));
                %
                sigma1 = fix(S1/2);
                sigma2 = fix(S2/2);
                [uxx1,uxy1,uyy1] = gsderiv(u_,sigma1,2);
                [uxx2,uxy2,uyy2] = gsderiv(u_,sigma2,2);
                %
                mode = 'Ridge';
                %
                % dirty solution
                hw = waitbar(0,['Conducting segmentation, please wait']);
                for x=1:size(u,1)
                    if ~isempty(hw), waitbar(x/size(u,1),hw); drawnow, end;                            
                    for y=1:size(u,2)
                        H = [uxx1(x,y) uxy1(x,y); uxy1(x,y) uyy1(x,y)];
                        [~,D] = eig(H); 
                        upp = D(1,1); 
                        uqq = D(2,2);
                        %
                        switch mode 
                            case 'Ridge'
                                if upp < 0
                                    str1(x,y) = abs(upp);
                                end
                            case 'Peak'
                            if upp < 0 && uqq < 0
                                str1(x,y) = sqrt(upp*uqq);
                            end             
                        end
                        %
                        H = [uxx2(x,y) uxy2(x,y); uxy2(x,y) uyy2(x,y)];
                        [~,D] = eig(H); 
                        upp = D(1,1); 
                        uqq = D(2,2);
                        %
                        switch mode 
                            case 'Ridge'
                                if upp < 0
                                    str2(x,y) = abs(upp);
                                end
                            case 'Peak'
                            if upp < 0 && uqq < 0
                                str2(x,y) = sqrt(upp*uqq);
                            end             
                        end                                                         
                    end                    
                end
                if ~isempty(hw), delete(hw), drawnow; end;
                %                             
                z = a*nth1.*map((str1),0,1) + (1-a)*nth2.*map((str2),0,1);
                %z = (a*map((str1),0,1) + (1-a)*map((str2),0,1))>0.5;
                z = imclose(z>t,strel('disk',fix(max(1,min(S1,S2)/2)),0));
                %
                nucs = z; 
                                    
    nucs = nucs | mask;
     
    nucs = imclearborder(nucs);
    nucs = bwareaopen(nucs,81);        
                                                                
        end
        
        