function [datas, captions, table_names, fig] = analyze_MNEPR(obj,~,~) 

     datas = [];
     captions = [];
     table_names = 'default';
     fig = [];
     
     sgm = obj.do_MNEPR_Segmentation(false);
     %
     u = single(obj.imgdata);
     [sX,sY,sC,sZ,sT] = size(u);

     data = [];
     for tind = 1:sT
         img = squeeze(u(:,:,1,1,tind));
         puncta = sgm(:,:,1,1,tind);
         nuclei = sgm(:,:,2,1,tind);
         %
         nuc_labs = bwlabel(nuclei); 
         stats_n = regionprops(nuc_labs,'Area','Centroid');
         punct_labs = bwlabel(puncta); 
         stats_p = regionprops(punct_labs,'Area','Centroid','PixelList','EquivDiameter');
         %         
         nuc_dilated = imdilate(nuclei,strel('disk',4));
         %
         % total number of puncta
         % total intensity of puncta
         % total area of puncta
         % total area of nuclei
         % distance between the nuclei
         % number of puncta at nuclei
         % area of puncta at nuclei
         % total intensity of puncta at nuclei
         tot_punct_number = numel(stats_p);
         tot_nuc_number = numel(stats_n);         
         tot_nuc_area = sum(nuclei(:));
         tot_punct_area = sum(puncta(:));          
         sample = img(0~=puncta);            
         tot_punct_intensity = sum(sample(:));
         %   
         number_of_puncta_at_nuclei = 0;
         area_of_puncta_at_nuclei = 0;
         total_intensity_of_puncta_at_nuclei = 0;
         %         
         for p = 1:tot_punct_number
             p_img = (punct_labs==p);
             if 0~=sum(sum(nuc_dilated.*p_img))
                 number_of_puncta_at_nuclei = number_of_puncta_at_nuclei + 1;
                 area_of_puncta_at_nuclei = area_of_puncta_at_nuclei + sum(sum(p_img));
                 sample = img(0~=p_img);            
                 total_intensity_of_puncta_at_nuclei = total_intensity_of_puncta_at_nuclei + sum(sample(:));                 
             end
         end
         %
         % distance between nuclei - start
         Dmaxnuc = 0;
         for n1 = 1:tot_nuc_number
             for n2 = 1:tot_nuc_number
                 XC1 = stats_n(n1).Centroid(1);
                 YC1 = stats_n(n1).Centroid(2);
                 XC2 = stats_n(n2).Centroid(1);
                 YC2 = stats_n(n2).Centroid(2);                                  
                 d = norm([XC1 YC1] - [XC2 YC2]);
                 if d > Dmaxnuc
                     Dmaxnuc = d;
                 end
             end
         end
         % distance between nuclei - ends
         record = {obj.current_filename, ...
             tind, ...
             tot_punct_number, ...
             tot_nuc_number, ...
             tot_nuc_area, ...
             tot_punct_area, ...
             tot_punct_intensity, ...                
             number_of_puncta_at_nuclei, ...
             area_of_puncta_at_nuclei, ...
             total_intensity_of_puncta_at_nuclei, ...
             Dmaxnuc}; 
             %
             data = [data; record];
     end
     %
     datas = data;
     table_names = {'MNEPR'};
     %
     icyvol = zeros(sX,sY,3,sZ,sT,'uint16');     
     icyvol(:,:,1,1,:) = uint16(u);
     icyvol(:,:,2,1,:) = uint16(sgm(:,:,1,1,:));
     icyvol(:,:,3,1,:) = uint16(sgm(:,:,2,1,:));
     fig = icyvol;
          
     captions = {'filename','time_frame', ...
             'tot_punct_number', ...
             'tot_nuc_number', ...
             'tot_nuc_area', ...
             'tot_punct_area', ...
             'tot_punct_intensity', ...                
             'number_of_puncta_at_nuclei', ...
             'area_of_puncta_at_nuclei', ...
             'total_intensity_of_puncta_at_nuclei', ...
             'Dmaxnuc'};                                      
end
