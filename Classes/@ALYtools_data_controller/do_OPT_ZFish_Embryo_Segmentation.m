function sgm = do_OPT_ZFish_Embryo_Segmentation(obj,send_to_Icy,~) 
            
        hw = waitbar(0,'Segmenting OPT ZFish_Embryos, please wait');
        if ~isempty(hw), waitbar(0.2,hw); drawnow, end
        
        tic        
            sgm = [];
            %
            if isempty(obj.imgdata), return, end
            %
            if 3==numel(size(obj.imgdata))
                I = single(squeeze(obj.imgdata));
                % [sx,sy,sz]=size(I);
            elseif 4==numel(size(obj.imgdata))
                I = single(squeeze(obj.imgdata(:,:,:,obj.OPT_ZFish_Embryo_channel_body)));
                head = single(squeeze(obj.imgdata(:,:,:,obj.OPT_ZFish_Embryo_channel_rostral)));
                tail = single(squeeze(obj.imgdata(:,:,:,obj.OPT_ZFish_Embryo_channel_posterior)));                
            end
            %
            % the image at imgdata is used to create multiple embryo volumes
            obj.M_imgdata = cell(0);
            obj.M_sgm = cell(0);                                                
            %
            S1 = round(obj.OPT_ZFish_Embryo_sgm_primary_scale/obj.microns_per_pixel);
            S2 = round(obj.OPT_ZFish_Embryo_sgm_K21*S1);
            S3 = round(obj.OPT_ZFish_Embryo_sgm_K31*S1);
            U1 = gauss3filter(I,S1);
            U2 = gauss3filter(I,S2);
            U3 = gauss3filter(I,S3);
            %
            u1 = (U1-U2)./U2;
            u2 = (U2-U3)./U3;
            % pixelwise max?
            % z = reshape(max(u1(:),u2(:)),[sy_m,sx_m,sz_m]);
            % or, try proportionality ?
            a1 = obj.OPT_ZFish_Embryo_sgm_a1;
            t  = obj.OPT_ZFish_Embryo_sgm_t;
            %
            z = u1*a1 + u2*(1-a1);                            
            z = z > t;                           
            % remove too thin objects
            se = strel('sphere',round(S1/2));
            z1 = imerode(z,se);
            z2 = imdilate(z1,se);
            embr_sgm = bwmorph3(z2,'majority');
            embr_sgm = bwmorph3(embr_sgm,'fill');
            
            min_embr_size = round(obj.OPT_ZFish_Embryo_sgm_min_vol/obj.microns_per_pixel^3);
            embr_sgm = bwareaopen(embr_sgm,min_embr_size,26);

            L_embr = bwlabeln(embr_sgm);
            stats_embr = regionprops3(L_embr,'BoundingBox','Centroid','VoxelList');
            n_embr = size(stats_embr,1);
            
            for k=1:n_embr
                    p = stats_embr.VoxelList{k};
                    y = round(p(:,1));
                    x = round(p(:,2));
                    z = round(p(:,3));
                    minx=min(x(:)); maxx=max(x(:));
                    miny=min(y(:)); maxy=max(y(:));
                    minz=min(z(:)); maxz=max(z(:));                
                    rx=minx:maxx;
                    ry=miny:maxy;
                    rz=minz:maxz;                
                    %
                    obj.M_sgm{k} = embr_sgm(rx,ry,rz);
                    %
                    if 4==numel(size(obj.imgdata)) % this is 3-channels case
                        body_k = I(rx,ry,rz);
                        head_k = head(rx,ry,rz);
                        tail_k = tail(rx,ry,rz);
                        v_k = zeros(size(body_k,1),size(body_k,2),size(body_k,3),3);
                        v_k(:,:,:,1) = body_k;
                        v_k(:,:,:,2) = head_k;
                        v_k(:,:,:,3) = tail_k;                        
                        obj.M_imgdata{k} = v_k;                        
                    else                    
                        obj.M_imgdata{k} = I(rx,ry,rz);
                    end                    
                    %
                    if send_to_Icy
                        iv = zeros(maxx-minx+1,maxy-miny+1,2,maxz-minz+1,1);
                        iv(:,:,1,:,:) = obj.M_imgdata{k};
                        iv(:,:,2,:,:) = obj.M_sgm{k};
                        icy_imshow(uint16(iv));                    
                    end
                    if ~isempty(hw), waitbar(k/n_embr,hw); drawnow, end  
            end                                                         
            if ~isempty(hw), delete(hw), drawnow; end
        toc        
end        