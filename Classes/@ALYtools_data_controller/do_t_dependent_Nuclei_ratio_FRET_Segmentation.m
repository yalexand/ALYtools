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
            
            z = sum(obj.imgdata,5);
            i1 = squeeze(z(:,:,:,1));
            i2 = squeeze(z(:,:,:,2));
            registration_exclusion_mask = (i1==0) | (i2==0);
                        
            % segment
            % sT = 32; % .. for debugging
            if 2==sC
                icyvol = zeros(sX,sY,3,1,sT,'single');
            elseif 3==sC
                icyvol = zeros(sX,sY,5,1,sT,'single');
            end
            
            for k=1:sT                
               ud = single(squeeze(obj.imgdata(:,:,1,D_channel,k)));
               ua = single(squeeze(obj.imgdata(:,:,1,A_channel,k)));
                                             
               if 0==sum(ud(:)) || 0==sum(ud(:)) , continue, end % ??
                 
               S = ceil(16*fac);
               smooth_scale = ceil(4*fac);
               max_hole_size = smooth_scale^2;
                              
               % mask zeros that may ocprv due to usage of warp transform
               if 0~=sum(sum((0==ud)))
                   z = imdilate(ud,strel('disk',S));
                   s = z(0==ud);
                   ud(0==ud) = min(s(:))/1.75;
               end
               if 0~=sum(sum((0==ua)))
                   z = imdilate(ua,strel('disk',S));
                   s = z(0==ua);
                   ua(0==ua) = min(s(:))/1.75;
               end               
               % mask zeros that may ocprv due to usage of warp transform
                                             
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
               nukes = bwareaopen(nukes,ceil(100*fac*fac)); % safety
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
                     
                nukes = nukes &~ registration_exclusion_mask;     
                nukes = bwareaopen(nukes,ceil(100/4*fac*fac)); % 100 divided by 4, because resized back
            
                icyvol(:,:,1,1,k) = ud;
                icyvol(:,:,2,1,k) = ua;               
                icyvol(:,:,3,1,k) = nukes;
                
                disp(['basic nukes ' num2str(k) ' ' num2str(mean(ud(:))) ' ' num2str(mean(ua(:))) ' ' num2str(sum(nukes(:)))]);
            end
                
%             %%%%%%%%%%%%%%%% if k>=2, ensure that "separated stays separated" by (k-1)-th             
            for k=2:sT  
                    prv = bwlabel(icyvol(:,:,3,1,k-1));
                    cur = bwlabel(icyvol(:,:,3,1,k));

                    stats_prv = regionprops(prv,'Centroid');
                      
                    %centroids image of (k-1) nuclei
                    cim_prv = zeros(size(cur));
                     for m=1:max(prv(:))
                        x = fix(stats_prv(m).Centroid(2));
                        y = fix(stats_prv(m).Centroid(1));
                        cim_prv(x,y) = m;
                     end
                     %
                     % main loop 
                     cur_fixed = cur.*0;
                     for m=1:max(cur(:))
                         z=(cur==m).*cim_prv;
                         s=z(z~=0);
                         prevlab = unique(s(:)); % previous centroid labels within current body
                         if numel(prevlab)>1 % then (cur==m) needs to be corrected

                                D = bwdist(z>0,'euclidean'); %distance map 
                                z = D.*(cur==m);
                                D(~z)=-Inf;
                                L = watershed(D);                                                                                
                                % remove background    
                                stats = regionprops(L,'Area');    
                                bckgind = find([stats.Area]==max([stats.Area]));
                                L(L==bckgind) = 0;
                                %result = (L>0) + (cur>0);
                                %icy_imshow(result,['frame ' num2str(k) ', object ' num2str(m)]);
                                cur_fixed = cur_fixed + (L>0);
                                disp(['fixed! #labels ' num2str(numel(prevlab))]);
                         else
                                cur_fixed = cur_fixed + (cur==m);
                         end  
                     end
                     icyvol(:,:,3,1,k) = cur_fixed;
                     disp(['separated stays separated ' num2str(k) ' ' num2str(sum(cur(:)))]);
            end
%             %%%%%%%%%%%%%%%% if k>=2, ensure that "separated stays separated" by (k-1)-th           

            if 3==sC
            for k=1:sT
                       uref = single(squeeze(obj.imgdata(:,:,1,3,k)));
                       %
                       if true % hey, - would it be better to substitute threshold-type segmentation here? 
                           S = ceil(18*fac); % ..say..
                           smooth_scale = ceil(2*fac);
                           max_hole_size = smooth_scale^2;

                           % skip? - mask zeros that may ocprv due to usage of warp transform
                           K  = 2.5;
                           % t = 0.055; 
                           t = 0.1; % less generous
                           z = nonlinear_tophat(uref,S,K)-1;
                           z(z<t)=0;                 
                           z = imopen(z,strel('disk',smooth_scale,0));
                       end
                       z = z | nukes;
                         % fill SMALL holes                             
                         z1 = imfill(z,'holes') - z;                    
                         z2 = bwareaopen(z1,max_hole_size);
                         z1(z2) = 0;
                         z = z | z1;
                         % fill small holes - ends                                
                       %               
                       cell_bodies = z;
                       %break cell body clumps
                            z2 = bwmorph(nukes,'thicken',Inf);
                            sep_lines = ~z2;
                            sep_lines(~cell_bodies)=0;
                            cell_bodies(sep_lines)=0;

                            % remove orphan pieces of cellular stuff..
                            L_n = bwlabel(nukes);
                            stats_n = regionprops(L_n,'Area','Centroid');    
                            L_c = bwlabel(cell_bodies);                                    
                            %
                            for n = 1:numel(stats_n)
                                xcn = fix(stats_n(n).Centroid(2));
                                ycn = fix(stats_n(n).Centroid(1));
                                cell_label = L_c(xcn,ycn);
                                L_c(L_c==cell_label)=0;
                            end            
                            cell_bodies(L_c~=0)=0;          
                            % remove orphan pieces of cellular stuff - ends
                            icyvol(:,:,4,1,k) = uref;
                            icyvol(:,:,5,1,k) = cell_bodies;
                            disp(['cell segmentation ' num2str(k) ' ' num2str(sum(uref(:))) ' ' num2str(sum(cell_bodies(:)))]);
            end
            end
            %
            % fix possible missing frames by assigning a previous one
            for k=1:sT
                nukes = squeeze(icyvol(:,:,3,1,k));
                if 0==sum(nukes(:)) && k > 1
                    ud_prev     = squeeze(icyvol(:,:,1,1,k-1));
                    ua_prev     = squeeze(icyvol(:,:,2,1,k-1));                    
                    nukes_prev  = squeeze(icyvol(:,:,3,1,k-1));
                    if 3==sC
                        uref_prev = icyvol(:,:,4,1,k-1);
                        cell_body_prev = icyvol(:,:,5,1,k-1);
                    end                    
                    %
                    icyvol(:,:,1,1,k) = ud_prev;
                    icyvol(:,:,2,1,k) = ua_prev;
                    icyvol(:,:,3,1,k) = nukes_prev;
                    if 3==sC
                        icyvol(:,:,4,1,k) = uref_prev;
                        icyvol(:,:,5,1,k) = cell_body_prev;
                    end                                        
                    disp(['fixed missing frame ' num2str(k)]);
                end
            end                            
            % fix possible missed frames by assigning a previous one            
            %
            sgm = uint16(icyvol);
                                                           
                if send_to_Icy                
                    try
                        notification = [obj.prvrent_filename ' - Segmentation:SmallNucleiTimeStack_Segmentation'];
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
          