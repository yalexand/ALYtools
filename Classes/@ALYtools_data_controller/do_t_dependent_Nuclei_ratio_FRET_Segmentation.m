        function sgm = do_t_dependent_Nuclei_ratio_FRET_Segmentation(obj,send_to_Icy,~) 
                       
            [sX,sY,sZ,sC,sT] = size(obj.imgdata);
            if sZ~=1
                obj.imgdata = reshape(obj.imgdata,[sX,sY,sT,sC,sZ]);
                [sX,sY,sZ,sC,sT] = size(obj.imgdata);
            end
            % to reduce to "T" dependence
            
            % segment
            %sT = 6; % 256 is record, July 5; % debug
            icyvol = zeros(sX,sY,3,1,sT);
            for k=1:sT                
               ud = single(squeeze(obj.imgdata(:,:,1,1,k)));
               ua = single(squeeze(obj.imgdata(:,:,1,2,k)));
               if 0==sum(ud(:)) || 0==sum(ud(:)) , continue, end % ??
                 
               S = 16;  
               K  = 2.5;
               t = 0.055;

               % once - donor
               z = imresize(ud,[2*sX 2*sY],'bicubic');               
               z = nonlinear_tophat(z,S,K)-1;
               z(z<t)=0;                 
               z = imopen(z,strel('disk',4,0));
               %
                 % fill SMALL holes                             
                 z1 = imfill(z,'holes') - z;
                 z2 = bwareaopen(z1,16);
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
               z = imopen(z,strel('disk',4,0));               
               %
                 % fill SMALL holes                             
                 z1 = imfill(z,'holes') - z;
                 z2 = bwareaopen(z1,16);
                 z1(z2) = 0;
                 z = z | z1;
                 % fill small holes - ends                                
               %               
               nukes = z;               
               nukes = bwareaopen(nukes,100); % safety
               sgm_a = nukes;
               
               %%%%%%%%%%%%%%%%%
               z = sgm_d | sgm_a;
                 
                    % break nuclear clumps
                    z = imdilate(z,strel('disk',2,0));
                    z2 = bwlabel(z);
                    D = bwdist(~z2,'euclidean'); %distance map 
                       nuc_breacking_distmap_smoothing_scale = 4;
                       D = medfilt2(D,[nuc_breacking_distmap_smoothing_scale nuc_breacking_distmap_smoothing_scale]);
                    D = -D;
                    D(~z2) = -Inf;                                        
                    L = watershed(D);                                                                                
                    % remove background    
                    stats = regionprops(L,'Area');    
                    bckgind = find([stats.Area]==max([stats.Area]));
                    L(L==bckgind) = 0;
                    z = (L>0);
                    % break nuclear clumps - end                    
                    %            
                 nukes = z;
                    %
                 nukes = bwareaopen(nukes,100); % safety               

                nukes = imresize(double(nukes),[sX sY],'bicubic');
                nukes = (nukes>0.8); % binarize this way
                                             
                icyvol(:,:,1,1,k) = ud;
                icyvol(:,:,2,1,k) = ua;               
                icyvol(:,:,3,1,k) = nukes;
               k 
            end
            %
            sgm = uint8(icyvol);
                                                           
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

        