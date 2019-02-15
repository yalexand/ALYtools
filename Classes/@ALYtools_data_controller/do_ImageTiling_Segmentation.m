function icyvol = do_ImageTiling_Segmentation(obj,send_to_Icy,~) 
        
        org = [];

        single_image_input = false;
        
        if ~isempty(obj.imgdata) && numel(size(obj.imgdata))>2
            org = single(squeeze(obj.imgdata));
            single_image_input = true;
        elseif ~isempty(obj.M_imgdata)
            z = obj.M_imgdata{1};
            [sX,sY,] = size(z);
            org = zeros(sX,sY,numel(obj.M_imgdata),'single');
            for k=1:numel(obj.M_imgdata)
                org(:,:,k) = single(obj.M_imgdata{k});
                % sometimes it needs this instead !!!???
                %org(:,:,k) = single(obj.M_imgdata{numel(obj.M_imgdata)-k+1});
            end
        else
            return;
        end

        [sX,sY,nImg] = size(org);        
                
        if strcmp('brightfield',obj.ImageTiling_mode)
            
                % gradmod
                grad_scale = 3;
                GMD = zeros(size(org));
                %
                for k=1:nImg
                    u = squeeze(org(:,:,k));
                        [gx,gy] = gsderiv(u,grad_scale,1);
                        g = sqrt(gx.^2+gy.^2);
                        GMD(:,:,k) = sqrt(gx.^2+gy.^2);
                end
                % minimum intensity projection
                mip = squeeze(GMD(:,:,1));
                for k=2:nImg
                        g = squeeze(GMD(:,:,k));
                        mip = pixelwise_min(mip,g);
                end
                %icy_imshow(mip);                
                %subtract mip
                for k=1:nImg
                        g = GMD(:,:,k);
                        GMD(:,:,k) = g - mip;
                end            
        
                icyvol = zeros(sX,sY,2,1,nImg);            
                %
                for k=1:nImg                                          
                    S = 5;
                    K = 2.5;
                    smooth_scale = 2;
                    t = 0.1;                 
                    z = nonlinear_tophat(squeeze(GMD(:,:,k)),S,K)-1;
                    z = imclose(z>t,strel('disk',smooth_scale,0));
                    z = bwareaopen(z,smooth_scale*smooth_scale);
                    %
                    icyvol(:,:,1,1,k) = z;
                    icyvol(:,:,2,1,k) = squeeze(org(:,:,k));
                    [k nImg]
                end                                    
                                    
        elseif strcmp('bleached_fluor',obj.ImageTiling_mode)
                        
                icyvol = zeros(sX,sY,2,1,nImg);            
                obj.PR_ref_channel = 1;
                for k=1:nImg                                          
                    S = 5;
                    K = 2.5;
                    smooth_scale = 2;
                    t = 0.1;
                    z = nonlinear_tophat(squeeze(org(:,:,k)),S,K)-1;
                    z = imclose(z>t,strel('disk',smooth_scale,0));
                    z = bwareaopen(z,smooth_scale*smooth_scale);
                    %                    
                    icyvol(:,:,1,1,k) = z;
                    icyvol(:,:,2,1,k) = squeeze(org(:,:,k));
                    [k nImg]
                end                                    
        end        
        %icy_imshow(icyvol);
        %
        if single_image_input
            obj.imgdata = org;
        end
        %                                    
                if send_to_Icy                
                    try
                        notification = [obj.current_filename ' - Segmentation: ImageTiling_Segmentation'];
                        if isempty(obj.h_Icy_segmentation_adjustment)
                            obj.h_Icy_segmentation_adjustment = icy_imshow(icyvol,notification);                    
                        else
                            icy_imshow(obj.h_Icy_segmentation_adjustment,icyvol,notification);                    
                        end
                    catch
                        errordlg('problem with Icy, - might be not running');
                    end                
                end        
end                

        