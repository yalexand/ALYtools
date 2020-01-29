function [datas, captions, table_names, fig] = analyze_OPT_ZFish_Embryo(obj,~,~) 
     datas = [];
     captions = [];
     table_names = 'OPT Embryos quantification';
     fig = [];
             
     obj.do_OPT_ZFish_Embryo_Segmentation(false);
          
hw = waitbar(0,'quantifying OPT Embryos, please wait');

for m=1:numel(obj.M_imgdata)
TITLE = strrep([obj.current_filename ' #' num2str(m)],'_',' ');

    if 4==numel(size(obj.imgdata)) % this is 3-channels case
        I = obj.M_imgdata{m};
        embr = squeeze(I(:,:,:,1));
    else
        embr = obj.M_imgdata{m};
    end
            
    embr_sgm = obj.M_sgm{m};
    [sx,sy,sz]=size(embr);        
              
        embr_skel = Skeleton3D(embr_sgm);            
        [w,l,h] = size(embr_skel);

        THR = 50;
        [A,node,link] = Skel2Graph3D(embr_skel,THR);
        % this is needed to erase void branchpoints 
        skel2 = Graph2Skel3D(node,link,sx,sy,sz);
        [A2,node2,link2] = Skel2Graph3D(skel2,0);        
        embr_skel_pruned = Graph2Skel3D(node2,link2,w,l,h);
        %
        LJ = bwlabeln(bwmorph3(embr_skel_pruned,'endpoints'));
        s = regionprops3(LJ,'VoxelList');
        p = s.VoxelList;
        % x1 = p(1,2);
        % y1 = p(1,1);
        % z1 = p(1,3);
        % x2 = p(2,2);
        % y2 = p(2,1);
        % z2 = p(2,3);
        %visuals  
          v = zeros(sx,sy,2,sz,1);
          v(:,:,1,:,1) = embr_sgm;
          v(:,:,2,:,1) = embr_skel_pruned;
        %     v(x1,y1,2,z1,1) = 20;
        %     v(x2,y2,2,z2,1) = 20;    
%icy_imshow(v);

        % "shell"
        embr_shell = embr_sgm - imerode(embr_sgm,strel('sphere',1));
        % v = zeros(sx,sy,1,sz,1);
        % v(:,:,1,:,1) = embr_shell;
% icy_imshow(v);  

         L_embr = bwlabeln(embr_sgm);
         stats_embr = regionprops3(L_embr,'EquivDiameter','Solidity','VoxelList','BoundingBox','Centroid','PrincipalAxisLength');

         L_embr_skel_pruned = bwlabeln(embr_skel_pruned);
         stats_embr_skel_pruned = regionprops3(L_embr_skel_pruned,'EquivDiameter','Solidity','VoxelList','BoundingBox','Centroid','PrincipalAxisLength','EigenVectors','EigenValues');

         n_fragments = size(stats_embr,1);
         assert(n_fragments==size(stats_embr_skel_pruned,1));

         skel_to_fat_lut = zeros(1,n_fragments);
         for k=1:n_fragments
                p = stats_embr_skel_pruned.VoxelList{k};
                y = round(p(1,1));
                x = round(p(1,2));
                z = round(p(1,3));
                skel_to_fat_lut(k) = L_embr(x,y,z);
         end
                  
             for k=1:n_fragments                 
                    N_k = -123456789.;
                    V_k = 123456789.;
                    A_k = 123456789.;
                    sphericity_k = 123456789.;
                    EquivDiameter_k = 123456789.;
                    Solidity_k = 123456789.;
                    Rg_k = 123456789.;
                    PrincipalAxisLengthRatio_1_k = 123456789.;
                    PrincipalAxisLengthRatio_2_k = 123456789.;
                    Nskel_k = 123456789.;
                    tortuosity_k = 123456789.;
                    pruning_excess_k = 123456789.;
                    skel_k_effective_plane_distance = 123456789.;        
                    Rg_skel_k = 123456789.;
                    distance_to_surface_k_mean = 123456789.;
                    distance_to_surface_k_variation = 123456789.;
                    distance_to_surface_k_skewness = 123456789.;
                    distance_to_surface_k_kurtosis = 123456789.;
                    u_surface_k_variation = 123456789.;
                    u_surface_k_skewness = 123456789.;
                    u_surface_k_kurtosis = 123456789.;       
                    u_volume_k_variation = 123456789.;
                    u_volume_k_skewness = 123456789.;
                    u_volume_k_kurtosis = 123456789.;
                    N_branchpoints_k = 123456789.;
                    head_x_tail_k = 123456789.;

                pk = stats_embr.VoxelList{k};
                x = pk(:,1);
                y = pk(:,2);
                z = pk(:,3);
                N_k = length(x);
                shp = alphaShape(x,y,z);
                A_k = surfaceArea(shp);
                V_k = volume(shp);
                K = (4/3*pi)^(2/3)/(4*pi);
                sphericity_k = K*A_k/(V_k^(2/3));
                % sphericity_k = 1/sphericity_k; % :)
                Rc_k = stats_embr.Centroid(k,:);
                Rg_k = GyrationRadius(x,y,z,Rc_k(1),Rc_k(2),Rc_k(3));
                %                 
                EquivDiameter_k = stats_embr.EquivDiameter(k);
                Solidity_k = stats_embr.Solidity(k);        
                %
                L = stats_embr.PrincipalAxisLength;
                PrincipalAxisLengthRatio_1_k = L(2)/L(1);
                PrincipalAxisLengthRatio_2_k = L(3)/L(1);
                %
                % skeleton related - via distance map
                sklab_k = skel_to_fat_lut(k); % label of volume fragment corresponding to k
                skvol_k = L_embr_skel_pruned==k;
                embr_shell_k = embr_shell.*(sklab_k==L_embr);
                %
                % count number of branchpoints
                [~,nodes_k,links_k] = Skel2Graph3D(skvol_k,0);
                N_branchpoints_k = length(nodes_k)-sum(cell2mat({nodes_k(:).ep}));
                show_graph(nodes_k,links_k,sx,sy,sz,TITLE);                
                %
                dmap_k = bwdist(skvol_k);
                s_k = dmap_k.*embr_shell_k;
                s_k = s_k(s_k~=0);
                distance_to_surface_k_mean = mean(s_k(:));
                distance_to_surface_k_variation = std(s_k(:))/distance_to_surface_k_mean;
                distance_to_surface_k_skewness = skewness(s_k(:),1); % not adjusted for bias
                distance_to_surface_k_kurtosis = kurtosis(s_k(:),1);
                %        
                % do the same for inensity image - shell only
                s_k = embr.*embr_shell_k;
                s_k = s_k(s_k~=0);
                s_k = s_k - min(s_k(:));
                u_surface_k_variation = std(s_k(:))/mean(s_k(:));
                u_surface_k_skewness = skewness(s_k(:),1); % not adjusted for bias
                u_surface_k_kurtosis = kurtosis(s_k(:),1);
                %
                % do the same for inensity image - the whole volume
                s_k = embr.*(L_embr==k);
                s_k = s_k(s_k~=0);
                s_k = s_k - min(s_k(:));                
                u_volume_k_variation = std(s_k(:))/mean(s_k(:));
                u_volume_k_skewness = skewness(s_k(:),1); % not adjusted for bias
                u_volume_k_kurtosis = kurtosis(s_k(:),1);
                %                        
                Nskel_k = sum(skvol_k(:));                              % #voxels in pruned k-th skel        
                Nskel_k_non_pruned = sum((L_embr==k).*embr_skel,'All'); % #voxels in non-pruned k-th skel
                pruning_excess_k = (Nskel_k_non_pruned - Nskel_k)/Nskel_k_non_pruned;
                %
                % to do:
                % intensity related - cross-correlation for 2 ref ?

                try        
                        % to do: skeleton related - plane/center deviations via eigen vectors
                        skvol_k_labelled = bwlabeln(skvol_k);
                        s = regionprops3(skvol_k_labelled,'VoxelList','EigenVectors','Centroid');
                        p = cell2mat(s.VoxelList); % must be only 1        
                        [pNorm,~,p0] = bestfitplane(p);
                        skel_k_effective_plane_distance = 0;
                        for t=1:size(p,1)
                            skel_k_effective_plane_distance = skel_k_effective_plane_distance + abs((p(t,:)-p0)*pNorm);
                        end
                        %
                        skel_k_effective_plane_distance = skel_k_effective_plane_distance/size(p,1);        
                        Rg_skel_k = GyrationRadius(p(:,1),p(:,2),p(:,3),p0(1),p0(2),p0(3));
                        %
                        % tortuosity
                            skvol_k_pruned_ends_labelled = bwlabeln(bwmorph3(skvol_k,'endpoints'));
                            s = regionprops3(skvol_k_pruned_ends_labelled,'VoxelList');
                            p = s.VoxelList; 
                            if 2==size(p,1)
                                x1 = p(1,2);
                                y1 = p(1,1);
                                z1 = p(1,3);
                                x2 = p(2,2);
                                y2 = p(2,1);
                                z2 = p(2,3);
                                skel_ends_distance_k = norm([x1 y1 z1] - [x2 y2 z2]);
                            else
                                skel_ends_distance_k = Rg_skel_k;
                            end
                            tortuosity_k = Nskel_k/skel_ends_distance_k;
                        % tortuosity
                catch
                end

if 4==numel(size(obj.imgdata)) % this is 3-channels case
            I = obj.M_imgdata{m};
            head = squeeze(I(:,:,:,2));
            tail = squeeze(I(:,:,:,3));    
            s = 2*round(obj.OPT_ZFish_Embryo_sgm_primary_scale/obj.microns_per_pixel);
            %
            z = gauss3filter(head,s);
            [~,idx]=maxN(z); % first found - but it is OK
            z = zeros(size(head));
            z(idx(1),idx(2),idx(3)) = 1;
            z = bwdist(z);
            dmap_head = max(z(:))-z;            
            %
            z = gauss3filter(tail,s);
            [~,idx]=maxN(z); % first found - but it is OK
            z = zeros(size(tail));
            z(idx(1),idx(2),idx(3)) = 1;
            z = bwdist(z);
            dmap_tail = max(z(:))-z;                        
            %
            s_head = head.*dmap_head.*embr_sgm/sum(head.*embr_sgm,'All');
            s_tail = tail.*dmap_tail.*embr_sgm/sum(tail.*embr_sgm,'All');
            %
            % corr.coeff
            head_x_tail_k = corr2(s_head(:),s_tail(:));
            %
            v = zeros(size(head,1),size(head,2),4,size(head,3),1);            
            v(:,:,1,:,1) = head.*embr_sgm;
            v(:,:,2,:,1) = dmap_head.*embr_sgm;
            v(:,:,3,:,1) = tail.*embr_sgm;
            v(:,:,4,:,1) = dmap_tail.*embr_sgm;                        
            icy_imshow(uint16(v),[TITLE ' corr2 = ' num2str(head_x_tail_k)]);
end                
                rec_k = {obj.current_filename, ...
                    m, ...
                    N_k, ...
                    V_k, ...
                    A_k, ...
                    sphericity_k, ...
                    EquivDiameter_k, ...
                    Solidity_k, ...
                    Rg_k, ...
                    PrincipalAxisLengthRatio_1_k, ...
                    PrincipalAxisLengthRatio_2_k, ...
                    Nskel_k, ...
                    tortuosity_k, ...
                    pruning_excess_k, ...
                    skel_k_effective_plane_distance, ...        
                    Rg_skel_k, ...                    
                    distance_to_surface_k_mean, ...
                    distance_to_surface_k_variation, ...
                    distance_to_surface_k_skewness, ...
                    distance_to_surface_k_kurtosis, ...
                    u_surface_k_variation, ...
                    u_surface_k_skewness, ...
                    u_surface_k_kurtosis, ...       
                    u_volume_k_variation, ...
                    u_volume_k_skewness, ...
                    u_volume_k_kurtosis, ...
                    N_branchpoints_k, ...
                    head_x_tail_k, ...
                    };
                datas = [ datas; rec_k];
             end
             
             if 4==numel(size(obj.imgdata)) % this is 3-channels case
                        iv = zeros(sx,sy,5,sz,1);
                        iv(:,:,1,:,1) = embr;
                        iv(:,:,2,:,1) = head;
                        iv(:,:,3,:,1) = tail;                        
                        iv(:,:,4,:,1) = embr_sgm;
                        iv(:,:,5,:,1) = embr_skel_pruned;                        
             else
                        iv = zeros(sx,sy,3,sz,1);
                        iv(:,:,1,:,1) = embr;
                        iv(:,:,2,:,1) = embr_sgm;
                        iv(:,:,3,:,1) = embr_skel_pruned;
             end
             %
             if isempty(fig), fig = cell(1,numel(obj.M_imgdata)); end
             fig{m} = iv;
                                                                
    if ~isempty(hw), waitbar(m/numel(obj.M_imgdata),hw); drawnow, end                 
end
if ~isempty(hw), delete(hw), drawnow; end

% embr_datas
            captions = { ...
            'filename', ...    
            'index', ...    
            'N', ...
            'V', ...
            'A', ...
            'ShapeFactor', ...
            'EquivDiameter', ...
            'Solidity', ...
            'Rg', ...
            'PrincipalAxisLengthRatio_1', ...
            'PrincipalAxisLengthRatio_2', ...
            'Nskel', ...
            'tortuosity', ...
            'pruning_excess', ...
            'skel_effective_plane_distance', ...        
            'Rg_skel', ...
            'distance_to_surface_mean', ...
            'distance_to_surface_variation', ...
            'distance_to_surface_skewness', ...
            'distance_to_surface_kurtosis', ...
            'u_surface_variation', ...
            'u_surface_skewness', ...
            'u_surface_kurtosis', ...       
            'u_volume_variation', ...
            'u_volume_skewness', ...
            'u_volume_kurtosis', ...     
            'N_branchpoints', ... % in pruned skeleton
            'head_x_tail', ... % in pruned skeleton            
            };                                  
end
