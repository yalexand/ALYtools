        function sgm = do_NucleiTimeStack_Segmentation(obj,send_to_Icy,~) 
                       
            [sX,sY,sZ,sC,sT] = size(obj.imgdata);
            if sZ~=1
                obj.imgdata = reshape(obj.imgdata,[sX,sY,sT,sC,sZ]);
                [sX,sY,sZ,sC,sT] = size(obj.imgdata);
            end
            % to reduce to "T" dependence
            
            % segment
            %sT = 16; % debug
            icyvol = zeros(sX,sY,3,1,sT);
            for k=1:sT                
               ud = single(squeeze(obj.imgdata(:,:,1,1,k)));
               ua = single(squeeze(obj.imgdata(:,:,1,2,k)));
               if 0==sum(ud(:)) || 0==sum(ud(:)) , continue, end % ??
                 
               S = 8;  
               K  = 2.5;
               t = 0.055;
               % once - donor
               z = nonlinear_tophat(ud,S,K)-1;
               z(z<t)=0;                 
               z = imopen(z,strel('disk',2,0));               
               %
                 % fill SMALL holes                             
                 z1 = imfill(z,'holes') - z;
                 z2 = bwareaopen(z1,4);
                 z1(z2) = 0;
                 z = z | z1;
                 % fill small holes - ends                                
               %               
               nukes = z;               
               nukes = bwareaopen(nukes,9); % safety
               sgm_d = nukes;
               % another one - acceptor
               z = nonlinear_tophat(ua,S,K)-1;
               z(z<t)=0;                 
               z = imopen(z,strel('disk',2,0));               
               %
                 % fill SMALL holes                             
                 z1 = imfill(z,'holes') - z;
                 z2 = bwareaopen(z1,4);
                 z1(z2) = 0;
                 z = z | z1;
                 % fill small holes - ends                                
               %               
               nukes = z;               
               nukes = bwareaopen(nukes,9); % safety
               sgm_a = nukes;
               
               %%%%%%%%%%%%%%%%%
               %%%%%%%%%%%%%%%%%
                 z = sgm_d | sgm_a;
                    % break nuclear clumps
                    z = imdilate(z,strel('disk',1,0));
                    z2 = bwlabel(z);
                    D = bwdist(~z2); %distance map 
                       nuc_breacking_distmap_smoothing_scale = 2;
                       D = medfilt2(D,[nuc_breacking_distmap_smoothing_scale nuc_breacking_distmap_smoothing_scale]);
                    D = -D;
                    D(~z2) = -Inf;                                        
                    L = watershed(D);                                                                                
                    % remove background    
                    stats = regionprops(L,'Area');    
                    bckgind = find([stats.Area]==max([stats.Area]));
                    L(L==bckgind) = 0;
                    z = (L>0);
                    %            
                 nukes = z;
                    %
                 nukes = bwareaopen(nukes,25); % safety               
               %%%%%%%%%%%%%%%%%
               %%%%%%%%%%%%%%%%%
                                             
               icyvol(:,:,1,1,k) = ud;
               icyvol(:,:,2,1,k) = ua;               
               icyvol(:,:,3,1,k) = nukes;
               k 
            end                                    
            %                        
            sgm = icyvol;
                                                           
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

        