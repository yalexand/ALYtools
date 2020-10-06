function sgm = do_SIFNE_Segmentation(obj,send_to_Icy,~) 
                    
        tic        
            sgm = [];
            %
            if isempty(obj.imgdata), return, end
            %
            
            R                   = obj.SIFNE_LFT_OFT_Radius_of_Filter;
            NofOrientations_FT  = obj.SIFNE_LFT_OFT_Number_of_Filter_Orinetations;
                        
            u = single(obj.imgdata);
            s = u(u~=0);
            t = quantile(s(:),0.995);
            u(u>t)=t;
            OriginImg = uint8(map(u,0,255));
                                    
            [H, W] = size(OriginImg);
            OriginImg_Margin = uint8(zeros(H+R+R,W+R+R));
            OriginImg_Margin(R+1:H+R,R+1:W+R) = OriginImg;
            
            % try vulgar segmentation
            ROI_Mask = zeros(size(OriginImg_Margin));
            % segmnentation setups
            S = obj.SIFNE_vulgar_ROI_sgm_scale;
            t = obj.SIFNE_vulgar_ROI_sgm_threshold; 
            % segmentation setups - ends
            u = single(obj.imgdata);
            u = gsderiv(u,S,0);
            u = map(u,0,1);
            ROI_Mask(R+1:H+R,R+1:W+R) = u>t; 
            %
            % remove small patches
            min_patch_size = round(obj.SIFNE_vulgar_ROI_sgm_min_patch_size/((obj.microns_per_pixel).^2)); %remove objects smaller than XX square microns
            ROI_Mask = bwareaopen(ROI_Mask,min_patch_size);
            % should be OK ...
            ROI_Mask = fill_small_holes(ROI_Mask,round(min_patch_size/4));
            % vulgar segmentation - ends

            [OFT_Img, LFT_Img, LFT_Orientations] = LFT_OFT_mex(double(OriginImg_Margin),double(R),double(NofOrientations_FT),double(ROI_Mask));
                                                                    
            sigma = obj.SIFNE_SGM_Junction_Size;
            OFT_Img = gsderiv(OFT_Img,sigma,0); % smoother !!
                                   
            % global, dimensional threshold
            BW = OFT_Img/(2*R)>obj.SIFNE_SGM_Filaments_Threshold;
            
            BW = imclose(BW,strel('disk',1)); % :) smoother!!
            
            RawSke = bwmorph(BW,'thin',Inf);
            
            % remove empty patches
            L = bwlabel(ROI_Mask);
            for k=1:max(L(:))
                if 0==sum(RawSke.*(L==k),'All')
                    ROI_Mask = ROI_Mask &~ (L==k); 
                end
            end
                        
                    if send_to_Icy
                        [sx,sy]=size(ROI_Mask);
                        iv = zeros(H,W,3,1,1,'uint8');                        
                      iv(:,:,1,1,1)=OriginImg; %uint8(map(OriginImg_Margin,0,255));
                      iv(:,:,2,1,1)=uint8(ROI_Mask(R+1:H+R,R+1:W+R)*100);
                      iv(:,:,3,1,1)=uint8(RawSke(R+1:H+R,R+1:W+R)*200);                      
                      try                          
                        icy_imshow(iv);
                      catch
                          disp('cannot send image to Icy');
                      end
                      return;
                    end            
                                                                            
            sgm = zeros(size(RawSke,1),size(RawSke,2),2,1,1);
            sgm(:,:,1,1,1) = RawSke;
            sgm(:,:,2,1,1) = ROI_Mask;
                        
        toc
        
        disp('do_SIFNE_Segmentation');           
end  
