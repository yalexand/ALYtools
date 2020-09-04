function sgm = do_SIFNE_Segmentation(obj,send_to_Icy,~) 
                    
        tic        
            sgm = [];
            %
            if isempty(obj.imgdata), return, end
            %
            
            R                   = obj.SIFNE_LFT_OFT_Radius_of_Filter;
            NofOrientations_FT  = obj.SIFNE_LFT_OFT_Number_of_Filter_Orinetations;
            
            OriginImg = imadjust(im2uint8(obj.imgdata));
            [H, W] = size(OriginImg);
            OriginImg_Margin = uint8(zeros(H+R+R,W+R+R));
            OriginImg_Margin(R+1:H+R,R+1:W+R) = OriginImg;
            
            %ROI_Mask = ones(size(OriginImg_Margin)); % !!!
            % try vulgar segmentation
            ROI_Mask = zeros(size(OriginImg_Margin));
            % segmnentation setups
            S = obj.SIFNE_vulgar_ROI_sgm_scale;
            t = obj.SIFNE_vulgar_ROI_sgm_threshold; 
            % segmnentation setups - ends
            u = single(obj.imgdata);
            u = gsderiv(u,S,0);
            u = map(u,0,1);
            ROI_Mask(R+1:H+R,R+1:W+R) = u>t; 
            %
            % remove small patches
            min_patch_size = round(obj.SIFNE_vulgar_ROI_sgm_min_patch_size/((obj.microns_per_pixel).^2)); %remove objects smaller than 110 square microns
            ROI_Mask = bwareaopen(ROI_Mask,min_patch_size);
            % should be OK ...
            ROI_Mask = fill_small_holes(ROI_Mask,round(min_patch_size/4));
            %
            % vulgar segmentation - ends

                    if send_to_Icy
                        [sx,sy]=size(ROI_Mask);
                        iv = zeros(sx,sy,3,1,1,'uint8');                        
                      iv(:,:,1,1,1)=uint8(OriginImg_Margin); %uint8(map(OriginImg_Margin,0,255));
                      iv(:,:,2,1,1)=uint8(ROI_Mask*100);
                      try                          
                        icy_imshow(iv);
                      catch
                          disp('cannot send image to Icy');
                      end
                      return;
                    end            
                        
            [OFT_Img, LFT_Img, LFT_Orientations] = LFT_OFT_mex(double(OriginImg_Margin),double(R),double(NofOrientations_FT),double(ROI_Mask));
                                                                    
            % Automatic threshold button                    
            DefaultFactor = 1.42;
            I = mat2gray(OFT_Img);
            t = DefaultFactor*graythresh(I);
            if t>=1
                t = graythresh(I);
            end
            BW = im2bw(I,t);
            RawSke = bwmorph(BW,'thin',Inf);
                                                    
            %sgm = RawSke;
            sgm = zeros(size(RawSke,1),size(RawSke,2),2,1,1);
            sgm(:,:,1,1,1) = RawSke;
            sgm(:,:,2,1,1) = ROI_Mask;
                        
        toc
        
        disp('do_SIFNE_Segmentation');           
end  

