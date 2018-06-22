        function sgm = do_NucleiTimeStack_Segmentation(obj,send_to_Icy,~) 
                       
            [sX,sY,sZ,sC,sT] = size(obj.imgdata);
            
            % segment
            icyvol = zeros(sX,sY,2,1,sZ);            
            for k=1:sZ                
                u = single(squeeze(obj.imgdata(:,:,k)));
                if 0==sum(u(:)), continue, end;
                scale1 = 5;
                rel_bg_scale1 = 2;
                scale2 = 15;
                rel_bg_scale2 = 2;
                threshold = 0.1;
                smoothing = 1;
                min_size = 16;
                z = two_scale_nth(u,scale1,rel_bg_scale1,scale2,rel_bg_scale2,threshold,smoothing,min_size);
                z(z>0)=1;
                %
                    z2 = bwlabel(z);
                    D = bwdist(~z2); %distance map                
                        nuc_breacking_distmap_smoothing_scale = 3;
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
                    nukes = bwareaopen(nukes,9); % safety
                    icyvol(:,:,1,1,k) = u;
                    icyvol(:,:,2,1,k) = nukes;
                    k 
            end                                    
            %                        
            sgm = icyvol;
            %                                    
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

        