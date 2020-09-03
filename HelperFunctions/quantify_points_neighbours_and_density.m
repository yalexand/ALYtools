function out = quantify_points_neighbours_and_density(points,max_length_pixels)
%
%
%
XC = points(:,1);
YC = points(:,2);

    npoints = numel(XC);
    out = zeros(npoints,2);

            % adjacency matrices (symmetric)
                    % FOR SPEED ONLY
                    AJM_density = adjacency_matrix(XC,YC,'Delaunay');
                    AJM_nnghb = adjacency_matrix(XC,YC,'SOI');            
            % reasonable looking, but but slow
            %AJM_density = adjacency_matrix(XC,YC,'Gabriel');
            % AJM_nnghb = AJM_density;

            % distance matrix
            DM = squareform(pdist([XC YC]));
            % number of neighbours vector
            NNGHB = sum(AJM_density,1);
            iNNGHB = 1./NNGHB; % matrix containing inverse #nnghb
            %
            % correct AJM_density for distances that are too long...
             AJM_density_fixed = AJM_density;
            for kk=1:npoints
                for jj=1:npoints
                    if DM(kk,jj)>max_length_pixels
                        AJM_density_fixed(kk,jj)=0;
                    end
                end
            end
%             figure;            
%             ax=gca;
%             plot(XC,YC,'r.');
%             for kk=1:npoints
%                 for jj=1:npoints
%                     if kk<jj && 1==AJM_density(kk,jj)
%                         v1 = [XC(kk) XC(jj)];
%                         v2 = [YC(kk) YC(jj)];
%                         h=line(v1,v2);
%                         set(h,'Color','red');
%                         set(h,'LineStyle',':');
%                         set(h,'LineWidth',1);
%                         hold(ax,'on');
%                     end
%                     if kk<jj && 1==AJM_density_fixed(kk,jj)
%                         v1 = [XC(kk) XC(jj)];
%                         v2 = [YC(kk) YC(jj)];
%                         h=line(v1,v2);
%                         set(h,'Color','blue');
%                         hold(ax,'on');
%                     end
%                     if kk<jj && 1==AJM_nnghb(kk,jj)
%                         v1 = [XC(kk) XC(jj)];
%                         v2 = [YC(kk) YC(jj)];
%                         h=line(v1,v2);
%                         set(h,'Color','black');
%                         set(h,'LineStyle',':');
%                         set(h,'LineWidth',2);
%                         hold(ax,'on');
%                     end                           
%                  end
%             end
%             plot(ax,XC,YC,'r.');
%             hold(ax,'off');
%             daspect(ax,[1 1 1]);
            % correct AJM_density for distances that are too long...            
            
            % for "natural" # neighbours
            NNGHB_NNGHB = sum(AJM_nnghb,1);
            %
            for n=1:npoints
                % # neighbours
                nnghb = NNGHB_NNGHB(n);                                
                % estimate for cell density
                dstncs = DM(n,:).*AJM_density(n,:);
                dstncs = dstncs(0~=dstncs); % no zeros
                if ~isempty(dstncs)
                    davr = mean(dstncs);
                    density = (1 + iNNGHB(n))/(pi*davr^2);
                else
                    density = 0;
                end
                %
                out(n,1) = nnghb;
                out(n,2) = density;                
            end                        
end
