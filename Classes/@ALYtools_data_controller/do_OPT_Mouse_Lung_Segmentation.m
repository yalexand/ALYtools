        function sgm = do_OPT_Mouse_Lung_Segmentation(obj,send_to_Icy,~) 
            
        tic
        
            sgm = [];
            %
            if isempty(obj.imgdata), return, end
            %
            u = single(obj.imgdata);
            %
            overpeak_ratio  = 4;

            minval = min(u(:));

                    Nbins = 100;
                    [cnt,vls] = hist(u(:),Nbins);                    
                    
                    [pks,locs]= findpeaks(cnt,'MINPEAKDISTANCE',10);
                    %
                    candidate_index = min(locs);
                    %
                    if cnt(1) > cnt(candidate_index)
                        candidate_index = 1;
                    end
                    %
                    if candidate_index > 15; % ?? can happen if very crowded only...
                        t = vls(candidate_index); % or, one can try min on inverted histo..
                    else 
                        t = vls(candidate_index) + overpeak_ratio*(vls(candidate_index)-minval);
                    end                             

            z = (u>t);
            strel = ones(1,1,1);
            sgm = imclose(z,strel);
                            
            %%%%%%%%%%%%%%%%%%%%
%             S1 = 8;
%             S2 = 22;
%             U1 = gauss3filter(u,S1);
%             U2 = gauss3filter(u,S2);
%             z1 = u.*U1./U2./U2;
%             %
%             S1 = 2;
%             S2 = 6;
%             U1 = gauss3filter(u,S1);
%             U2 = gauss3filter(u,S2);
%             z2 = u.*U1./U2./U2;
% 
%             sgm = sgm & (max(z1,z2)>4);

            proc_time=toc;
            
            %%%%%%%%%%%%%%%%%%%%
                      
            % u(u<t)=0; % to ease life for Icy HKMeans
            
            [sX,sY,sZ] = size(u);
            icyvol = zeros(sX,sY,1,sZ,1,'uint16');
            
                if send_to_Icy                
                    try
                        icyvol(:,:,1,:,1) = uint16(u);
                        icyvol(:,:,2,:,1) = uint16(sgm);
                        %
                        notification = [obj.current_filename ... 
                            ' - Segmentation: OPT_Mouse_Lung_Segmentation' ' t = ' num2str(proc_time)];
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