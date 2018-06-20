        function sgm = do_DarkNuclei_Segmentation(obj,send_to_Icy,~) 
   
            u = single(obj.imgdata);
                    
            [sX,sY,sC,sZ,sT] = size(u);
            
            u = squeeze(sum(u,5));
            %
            cells = zeros(size(u)); 
            nuclei = zeros(size(u)); 
            %            
            %        
            % segment
            % backup PR setups                                    
            obj_PR_K = obj.PR_K;
            obj_PR_S1 = obj.PR_S1;
            obj_PR_S2 = obj.PR_S2;
            obj_PR_a = obj.PR_a;
            obj_PR_t = obj.PR_t;
            obj_PR_mode = obj.PR_mode;
            obj_PR_min_size = obj.PR_min_size;
            obj_PR_ref_channel = obj.PR_ref_channel;
            obj_imgdata = obj.imgdata;
            %
            
            % apply PR segmentation   
            vol = zeros(sX,sY,1,1,1);
            vol(:,:,1,1,1) = u;
            obj.PR_ref_channel = 1;                        

            % THIS WORKS WELL FOR ALL IMAGES EXCEPT 4.sdt
            % cells
            obj.imgdata = vol;            
            obj.PR_K = 2.5;
            obj.PR_S1 = 3;
            obj.PR_S2 = 7;
            obj.PR_a = 0.5;
            obj.PR_t = 0.01;
            obj.PR_mode = 'Ridge';
            obj.PR_min_size = 8;            
            cells = obj.do_PR_Segmentation(false);
            
            % nuclei;
            obj.imgdata = max(vol(:))-vol;            
            obj.PR_K = 2.5;
            obj.PR_S1 = 4;
            obj.PR_S2 = 8;
            obj.PR_a = 0.65;
            obj.PR_t = 0.002;
            obj.PR_mode = 'Peak';
            obj.PR_min_size = 8;            
            nuclei = obj.do_PR_Segmentation(false);
                        
            % restore PR setups
            obj.PR_K = obj_PR_K;
            obj.PR_S1 = obj_PR_S1;
            obj.PR_S2 = obj_PR_S2;
            obj.PR_a = obj_PR_a;
            obj.PR_t = obj_PR_t;
            obj.PR_mode = obj_PR_mode;
            obj.PR_min_size = obj_PR_min_size;
            obj.PR_ref_channel = obj_PR_ref_channel;
            obj.imgdata = obj_imgdata; 
                        
            %
            cell_smoothing_scale = 4;
            smthd_cells = gsderiv(u,cell_smoothing_scale,0);
            %            
            cell_overpeak_ratio = 2.5; % permissive
            %            
            minval = min(smthd_cells(:));
            %
            [cnt,vls] = hist(smthd_cells(:),100);
            [pks,locs]= findpeaks(cnt,'MINPEAKDISTANCE',10);
            %
            candidate_index = min(locs);
            %
            if cnt(1) > cnt(candidate_index)
                candidate_index = 1;
            end
            %
            t = vls(candidate_index) + cell_overpeak_ratio*(vls(candidate_index)-minval);
            %                        
            stuff1 = smthd_cells > t;
            stuff1 = imerode(stuff1,strel('disk',6,0));
            %
            stuff2 = imclose(cells,strel('disk',6,0));
            %
            
            %all_clumps = stuff1 | stuff2;         
            
            all_clumps = stuff2; % only this if background is bad
            imfill(all_clumps,'holes');
            
            max_cell_size = 430; %pix
            min_cell_size = 49; % ditto
            %
            all_clumps = bwareaopen(all_clumps,min_cell_size);
            %
            clumps = bwareaopen(all_clumps,max_cell_size);
            separate_cells = all_clumps.*(~clumps);
            
            z = nuclei.*clumps;
            z = bwmorph(z,'thicken',Inf);
            broken = z & clumps;
            
            resulting_cells = separate_cells | broken;
            resulting_cells  = bwareaopen(resulting_cells,min_cell_size);
            %
            % average intensity per cell should be high enough
            L = bwlabel(resulting_cells);
            for k=1:max(L(:))
                intensity = u(L==k);
                target = L==k;
                avrint = sum(intensity(:))/sum(target(:));
                if avrint < 1.2*t
                    resulting_cells(L==k) = 0;
                end
            end
            %
            % THIS WORKS WELL FOR ALL IMAGES EXCEPT 4.sdt
            
            
            
            
            % speical for 4.sdt
%             % cells
%             obj.imgdata = vol;            
%             obj.PR_K = 2.5;
%             obj.PR_S1 = 3;
%             obj.PR_S2 = 7;
%             obj.PR_a = 0.5;
%             obj.PR_t = 0.002;
%             obj.PR_mode = 'Ridge';
%             obj.PR_min_size = 8;            
%             cells = obj.do_PR_Segmentation(false);
%             
%             % nuclei;
%             obj.imgdata = max(vol(:))-vol;            
%             obj.PR_K = 2.5;
%             obj.PR_S1 = 6;
%             obj.PR_S2 = 14;
%             obj.PR_a = 0.2;
%             obj.PR_t = 0.0002;
%             obj.PR_mode = 'Peak';
%             obj.PR_min_size = 8;            
%             nuclei = obj.do_PR_Segmentation(false);
%                         
%             % restore PR setups
%             obj.PR_K = obj_PR_K;
%             obj.PR_S1 = obj_PR_S1;
%             obj.PR_S2 = obj_PR_S2;
%             obj.PR_a = obj_PR_a;
%             obj.PR_t = obj_PR_t;
%             obj.PR_mode = obj_PR_mode;
%             obj.PR_min_size = obj_PR_min_size;
%             obj.PR_ref_channel = obj_PR_ref_channel;
%             obj.imgdata = obj_imgdata; 
%                         
%             %
%             cell_smoothing_scale = 4;
%             smthd_cells = gsderiv(u,cell_smoothing_scale,0);
%             %            
%             cell_overpeak_ratio = 1.1; % permissive
%             %            
%             minval = min(smthd_cells(:));
%             %
%             [cnt,vls] = hist(smthd_cells(:),100);
%             [pks,locs]= findpeaks(cnt,'MINPEAKDISTANCE',10);
%             %
%             candidate_index = min(locs);
%             %
%             if cnt(1) > cnt(candidate_index)
%                 candidate_index = 1;
%             end
%             %
%             t = vls(candidate_index) + cell_overpeak_ratio*(vls(candidate_index)-minval);
%             %                        
%             stuff1 = smthd_cells > t;
%             stuff1 = imerode(stuff1,strel('disk',3,0));
%             %
%             stuff2 = imclose(cells,strel('disk',6,0));
%             %
%             all_clumps = stuff1 | stuff2;
%             
%             all_clumps = imfill(all_clumps,'holes');
%             all_clumps = imopen(all_clumps,strel('disk',4,0));
%                         
%             max_cell_size = 430; %pix
%             min_cell_size = 49; % ditto
%             %
%             all_clumps = bwareaopen(all_clumps,min_cell_size);
%             %
%             clumps = bwareaopen(all_clumps,max_cell_size);
%             separate_cells = all_clumps.*(~clumps);
%             
%             z = nuclei.*clumps;
%             z = bwmorph(z,'thicken',Inf);
%             broken = z & clumps;
%             
%             resulting_cells = separate_cells | broken;
%             resulting_cells  = bwareaopen(resulting_cells,min_cell_size);
%             %
%             % average intensity per cell should be high enough
%             L = bwlabel(resulting_cells);
%             for k=1:max(L(:))
%                 intensity = u(L==k);
%                 target = L==k;
%                 avrint = sum(intensity(:))/sum(target(:));
%                 if avrint < t
%                     resulting_cells(L==k) = 0;
%                 end
%             end
            
            % speical for 4.sdt            
            
            
            % speical for G1 7.sdt
            % cells
            obj.imgdata = vol;            
            obj.PR_K = 2.5;
            obj.PR_S1 = 3;
            obj.PR_S2 = 7;
            obj.PR_a = 0.5;
            obj.PR_t = 0.002;
            obj.PR_mode = 'Ridge';
            obj.PR_min_size = 8;            
            cells = obj.do_PR_Segmentation(false);
            
            % nuclei;
            obj.imgdata = max(vol(:))-vol;            
            obj.PR_K = 2;
            obj.PR_S1 = 3;
            obj.PR_S2 = 8;
            obj.PR_a = 0.7;
            obj.PR_t = 0.001;
            obj.PR_mode = 'Peak';
            obj.PR_min_size = 8;            
            nuclei = obj.do_PR_Segmentation(false);
                        
            % restore PR setups
            obj.PR_K = obj_PR_K;
            obj.PR_S1 = obj_PR_S1;
            obj.PR_S2 = obj_PR_S2;
            obj.PR_a = obj_PR_a;
            obj.PR_t = obj_PR_t;
            obj.PR_mode = obj_PR_mode;
            obj.PR_min_size = obj_PR_min_size;
            obj.PR_ref_channel = obj_PR_ref_channel;
            obj.imgdata = obj_imgdata; 
                                    
            %
            cell_smoothing_scale = 4;
            smthd_cells = gsderiv(u,cell_smoothing_scale,0);
            %            
            cell_overpeak_ratio = 1.7; % permissive
            %            
            minval = min(smthd_cells(:));
            %
            [cnt,vls] = hist(smthd_cells(:),100);
            [pks,locs]= findpeaks(cnt,'MINPEAKDISTANCE',10);
            %
            candidate_index = min(locs);
            %
            if cnt(1) > cnt(candidate_index)
                candidate_index = 1;
            end
            %
            t = vls(candidate_index) + cell_overpeak_ratio*(vls(candidate_index)-minval);
            %                        
            stuff1 = smthd_cells > t;
            stuff1 = imerode(stuff1,strel('disk',3,0));
            %
            stuff2 = imclose(cells,strel('disk',2,0));
            %
            all_clumps = stuff1 | stuff2;
                        
     % fill SMALL holes                             
     z = all_clumps;
     z1 = imfill(z,'holes') - z;
     z2 = bwareaopen(z1,256);
     z1(z2) = 0;
     z = z | z1;
     all_clumps = uint16(z);
     % fill small holes - ends                 

%             icyvol = zeros(sX,sY,3,1,1);
%             icyvol(:,:,1,1,1) = u;
%             icyvol(:,:,2,1,1) = nuclei;
%             icyvol(:,:,3,1,1) = all_clumps;
               
            all_clumps = imopen(all_clumps,strel('disk',2,0));
                                                                                  
            max_cell_size = 430; %pix
            min_cell_size = 49; % ditto
            %
            all_clumps = bwareaopen(all_clumps,min_cell_size);
            %
            clumps = bwareaopen(all_clumps,max_cell_size);
            separate_cells = all_clumps.*(~clumps);
            
            z = nuclei.*clumps;
            z = bwmorph(z,'thicken',Inf);
            broken = z & clumps;
            
            resulting_cells = separate_cells | broken;
            resulting_cells  = bwareaopen(resulting_cells,min_cell_size);
            
% % % % %             % speical for G1 7.sdt            

            
            

            icyvol = zeros(sX,sY,3,1,1);
            %
            icyvol(:,:,1,1,1) = u;
            icyvol(:,:,2,1,1) = resulting_cells;
            icyvol(:,:,3,1,1) = cells;
            %                        
            sgm = icyvol;
            %                                    
                if send_to_Icy                
                    try
                        notification = [obj.current_filename ' - Segmentation: DarkNuclei_Segmentation'];
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

        