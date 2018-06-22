function [datas, captions, table_names, fig] = analyze_NucleiTimeStack(obj,~,~) 

     datas = [];
     captions = [];
     table_names = 'default';
     fig = [];
     
     sgm = obj.do_NucleiTimeStack_Segmentation(false);
     %     
     [sX,sY,sZ,sC,sT] = size(obj.imgdata);
          
     N = floor(sZ/3);
          fig = zeros(sX,sY,2,1,N);
     data = [];
     for k=1:N
         i1 = 3*k-2;
         i2 = i1+1;
         r1 = sgm(:,:,1,1,i1);
         r2 = sgm(:,:,1,1,i2);
         r = r1./r2;
         sgm1 = sgm(:,:,2,1,i1);
         sgm2 = sgm(:,:,2,1,i2);
         z = sgm1 | sgm2;
                    z2 = bwlabel(z);
                    D = bwdist(~z2); %distance map                
                    nuc_breacking_distmap_smoothing_scale = 3;
                    D = medfilt2(D,[nuc_breacking_distmap_smoothing_scale nuc_breacking_distmap_smoothing_scale]);
                    D = -D;
                    D(~z2) = -Inf;                                        
                    L = watershed(D);                                                                                 
                    % remove background    
                    stats = regionprops(L,'Area');    
                    bckgind = find([stats.Area]==max([stats.Area]));
                    L(L==bckgind) = 0;
                    nukes = (L>0);
                    %
         nukes = bwareaopen(nukes,9); % safety
                    r(~nukes) = 0;
                    fig(:,:,1,1,k) = r;
                    fig(:,:,2,1,k) = nukes;     
        L = bwlabel(nukes);
        nnucs = max(L(:));
        nuc_data = zeros(nnucs,2);
        for n=1:nnucs
            sample=r(L==n);
            nuc_data(n,1)=length(sample(:));
            nuc_data(n,2)=mean(sample(:));
        end
        avr_area = mean(nuc_data(:,1));
        std_area = std(nuc_data(:,1));
        avr_ratio = mean(nuc_data(:,2));
        std_ratio = std(nuc_data(:,2));        
        data = [data; {nnucs avr_area std_area avr_ratio std_ratio}];
        k
     end
     captions = {'#nuclei','avr_area','std_area','avr_ratio','std_ratio'};
     datas = data;
end
