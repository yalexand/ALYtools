        function sgm = do_t_dependent_Nuclei_ratio_FRET_Segmentation(obj,send_to_Icy,~) 

            D_channel = obj.t_dependent_Nuclei_ratio_FRET_donor_channel;
            if 1 == D_channel
                A_channel = 2;
            else
                A_channel = 1;
            end
        
            % all mask sizes are tuned for 1.66 microns/pixel
            fac = 1.66/obj.microns_per_pixel;
            
            [sX,sY,sZ,sC,sT] = size(obj.imgdata);
            if sZ~=1
                obj.imgdata = reshape(obj.imgdata,[sX,sY,sT,sC,sZ]);
                [sX,sY,sZ,sC,sT] = size(obj.imgdata);
            end
            % to reduce to "T" dependence
            
            % segment
            % sT = 6; % .. for debugging
            icyvol = zeros(sX,sY,3,1,sT,'single');
            for k=1:sT                
               ud = single(squeeze(obj.imgdata(:,:,1,D_channel,k)));
               ua = single(squeeze(obj.imgdata(:,:,1,A_channel,k)));
               if 0==sum(ud(:)) || 0==sum(ud(:)) , continue, end % ??
                 
               %S = 16;  
               S = ceil(16*fac);
               smooth_scale = ceil(4*fac);
               max_hole_size = smooth_scale^2;
               
               K  = 2.5;
               % t = 0.055; % good for multiphoton data
               t = 0.1; % less generous
               
               % once - donor
               z = imresize(ud,[2*sX 2*sY],'bicubic');               
               z = nonlinear_tophat(z,S,K)-1;
               z(z<t)=0;                 
               z = imopen(z,strel('disk',smooth_scale,0));
               %
                 % fill SMALL holes                             
                 z1 = imfill(z,'holes') - z;                    
                 z2 = bwareaopen(z1,max_hole_size);
                 z1(z2) = 0;
                 z = z | z1;
                 % fill small holes - ends                                
               %               
               nukes = z;               
               nukes = bwareaopen(nukes,100); % safety
               sgm_d = nukes;
               % another one - acceptor
               z = imresize(ua,[2*sX 2*sY],'bicubic');
               z = nonlinear_tophat(z,S,K)-1;
               z(z<t)=0;                 
               z = imopen(z,strel('disk',smooth_scale,0));               
               %
                 % fill SMALL holes                             
                 z1 = imfill(z,'holes') - z;
                 z2 = bwareaopen(z1,max_hole_size);
                 z1(z2) = 0;
                 z = z | z1;
                 % fill small holes - ends                                
               %               
               nukes = z;               
               nukes = bwareaopen(nukes,ceil(100*fac*fac)); % safety
               sgm_a = nukes;
               
               %%%%%%%%%%%%%%%%%
               z = sgm_d | sgm_a;
                 
                    % break nuclear clumps - start
                    z = imdilate(z,strel('disk',max(2,ceil(2*fac)),0));
                    z2 = bwlabel(z);
                    D = bwdist(~z2,'euclidean'); %distance map 
                       nuc_breacking_distmap_smoothing_scale = max(6,ceil(6*fac));
                       D = medfilt2(D,[nuc_breacking_distmap_smoothing_scale nuc_breacking_distmap_smoothing_scale]);
                    D = -D;
                    D(~z2) = -Inf;                                        
                    L = watershed(D);                                                                                
                    % remove background    
                    stats = regionprops(L,'Area');    
                    bckgind = find([stats.Area]==max([stats.Area]));
                    L(L==bckgind) = 0;
                    nukes = (L>0);
                    %
                    % this paragraph corrects for excessive clump-breaking
                    stats = regionprops(z2, 'Area','ConvexArea'); 
                    area_ratio_threshold = 0.9; % looks alright
                    dont_break = ismember(z2,find([stats.Area]./[stats.ConvexArea]>area_ratio_threshold));  
                    % icy_imshow(z+dont_break); % to see what is going on                   
                    nukes = nukes | dont_break; % print the "dont_break" over, undoing excessive breaks
                    % break nuclear clumps - end                                               
                    %
                nukes = bwareaopen(nukes,ceil(100*fac*fac)); % safety               

                nukes = imresize(double(nukes),[sX sY],'bicubic');
                nukes = (nukes>0.8); % binarize this way
                                             
                icyvol(:,:,1,1,k) = ud;
                icyvol(:,:,2,1,k) = ua;               
                icyvol(:,:,3,1,k) = nukes;
               % 
               disp([num2str(k) ' ' num2str(mean(ud(:))) ' ' num2str(mean(ua(:))) ' ' num2str(sum(nukes(:)))]);
            end
            %
            % fix possible missing frames by assigning a previous one
            for k=1:sT
                nukes = squeeze(icyvol(:,:,3,1,k));
                if 0==sum(nukes(:)) && k > 1
                    ud_prev     = squeeze(icyvol(:,:,1,1,k-1));
                    ua_prev     = squeeze(icyvol(:,:,2,1,k-1));                    
                    nukes_prev  = squeeze(icyvol(:,:,3,1,k-1));
                    %
                    icyvol(:,:,1,1,k) = ud_prev;
                    icyvol(:,:,2,1,k) = ua_prev;
                    icyvol(:,:,3,1,k) = nukes_prev;
                    disp(['fixed missing frame ' num2str(k)]);
                end
            end                            
            % fix possible missed frames by assigning a previous one            
            %
            sgm = uint16(icyvol);
                                                           
                if send_to_Icy                
                    try
                        notification = [obj.current_filename ' - Segmentation:SmallNucleiTimeStack_Segmentation'];
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
        