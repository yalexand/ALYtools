        function sgm = do_OPT_Embryo_Segmentation(obj,send_to_Icy,~) 
            
        hw = waitbar(0,'Segmenting OPT Embryos, please wait');        
        
        tic        
            sgm = [];
            %
            if isempty(obj.imgdata), return, end
            %
            I = single(obj.imgdata);
            [sx,sy,sz]=size(I);
            %
            % the image at imgdata is used to create multiple embryo
            % volumes
            obj.M_imgdata = cell(0);
            obj.M_sgm = cell(0);            
            %
f = 4; % downsampling
                %
                Ids = resizeVolume(I,round([sx,sy,sz]/f));
                [sx,sy,sz]=size(Ids);
                % maximum intensity projection
                maipxy = -ones(sx,sy)*Inf;
                parfor m=1:sz
                    u = squeeze(Ids(:,:,m));
                    maipxy = pixelwise_max(maipxy,u);
                end
                maipxy = map(maipxy,0,1);
%icy_imshow(maipxy,'xy');

                maipxz = -ones(sx,sz)*Inf;
                parfor m=1:sy
                    u = squeeze(Ids(:,m,:));
                    maipxz = pixelwise_max(maipxz,u);                    
                end
                maipxz = map(maipxz,0,1);    
%icy_imshow(maipxz,'xz');

                maipyz = -ones(sy,sz)*Inf;
                parfor m=1:sx
                    u = squeeze(Ids(m,:,:));
                    maipyz = pixelwise_max(maipyz,u);
                    m
                end
                maipyz = map(maipyz,0,1);    
%icy_imshow(maipyz,'yz');

if ~isempty(hw), waitbar(0.5,hw); drawnow, end  

t = 0.1; % segmentation threshold
                z1 = maipxz>t;
                z1 = imclose(z1,strel('disk',1*round(16/f)));
                z1 = bwareaopen(z1,4*round((16/f)^2));
                xzline = sum(z1,1);

                z1 = maipyz>t;
                z1 = imclose(z1,strel('disk',1*round(16/f)));
                z1 = bwareaopen(z1,4*round((16/f)^2));
                yzline = sum(z1,1);
%figure();
%plot(1:sz,xzline,'b.-',1:sz,yzline,'r.-');grid on;

                z1 = maipxy>t;
                z1 = imclose(z1,strel('disk',1*round(16/f)));
                z1 = bwareaopen(z1,4*round((16/f)^2));
                xline = sum(z1,1);
                yline = sum(z1,2);
%figure();
%plot(1:sx,xline,'b.-',1:sy,yline,'r.-');grid on;

lxy_threshold = 0.5;
lxy_padfac = 1.6;

                z2 = yline > lxy_threshold;
                lz2 = bwlabel(z2);
                vz2 = (1:sy).*z2;
                vz2 = vz2(0~=vz2);
                c2 = mean(vz2);
                wy = round(sum(lz2~=0)*lxy_padfac/2);
                y1 = max(1,c2-wy);
                y2 = min(sy,c2+wy);
                %
                z1 = xline > lxy_threshold;
                lz1 = bwlabel(z1);
                vz1 = (1:sx).*z1;
                vz1 = vz1(0~=vz1);
                c1 = mean(vz1);
                wx = round(sum(lz1~=0)*lxy_padfac/2);
                x1 = round(max(1,c1-wx));
                x2 = round(min(sx,c1+wx));

%figure;
%xmask = zeros(1,sx);xmask(x1:x2) = 1;    
%ymask = zeros(1,sy);ymask(y1:y2) = 1;    
%plot(1:sx,xline,'b.-',1:sx,xmask,'b-',1:sy,yline,'r.-',1:sy,ymask,'r-');grid on;

zline_thershold = 4; % 0.5?
                z1 = xzline > zline_thershold;
                z2 = yzline > zline_thershold;

zpadfac = 1.15; % if embryos are close

                z = z1 | z2;
                lz = bwlabel(z);

%figure;
%plot(1:sz,xzline,'b.-',1:sz,yzline,'r.-');grid on;

                zcutsegments = cell(1,0);
                nembryos = max(lz(:));
                for m=1:nembryos
                    lk = (lz==m);
                    s = lk.*(1:sz); 
                    s = s(s~=0);
                    ck = mean(s);
                    wk = round(sum(lk)*zpadfac/2);
                    z1 = round(max(1,ck-wk));
                    z2 = round(min(sz,ck+wk));
                    zmask = zeros(1,sz);zmask(z1:z2) = m;
%hold on;
%plot(1:sz,zmask,'k.-');
%hold off;    
                    zcutsegments{m} = [z1 z2];
                end

                % create output images
                    for m=1:numel(zcutsegments)
                        zrange = zcutsegments{m}; 
                        z1 = zrange(1);
                        z2 = zrange(2);
                        sx_m = (x2-x1)*f+1;
                        sy_m = (y2-y1)*f+1;
                        sz_m = (z2-z1)*f+1;
                        rx = x1*f:x2*f;
                        ry = y1*f:y2*f;
                        rz = z1*f:z2*f;
                        v = zeros(sy_m,sx_m,1,sz_m,1);
                        v(:,:,1,:,1) = I(ry,rx,rz);
                        obj.M_imgdata{m} = v;                        
                        %
                        % segmentation
K = 2.5;
S1 = 5;
S2 = round(K*S1);
                            U1 = gauss3filter(squeeze(v),S1);
                            U2 = gauss3filter(squeeze(v),S2);
                            %
                            z = (U1-U2)./(U1+U2);
t  = 0.235; %?
                            embr_sgm = z > t;
                            embr_sgm = bwmorph3(embr_sgm,'majority');
                            embr_sgm = bwmorph3(embr_sgm,'fill');
min_embr_size = 40*40*40;            
                            embr_sgm = bwareaopen(embr_sgm,min_embr_size,26);
                            
                            obj.M_sgm{m} = embr_sgm;
                            obj.M_imgdata{m} = v;                                                    
                        % segmentation                        
                        %
                        if send_to_Icy
                            iv = zeros(sy_m,sx_m,2,sz_m,1);
                            iv(:,:,1,:,:) = v;
                            iv(:,:,2,:,:) = embr_sgm;
                            icy_imshow(uint16(iv));
                        end

                    end               
        disp(toc);
        if ~isempty(hw), delete(hw), drawnow; end        
        end        