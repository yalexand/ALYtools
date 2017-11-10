function [datas, captions, table_names, fig] = analyze_OPT_Mouse_Lung(obj,~,~) 
     datas = [];
     captions = [];
     table_names = 'default';
     fig = [];
             
     if ~isempty(obj.OPT_data_controller.proj) 
        % pay attention to SETTINGS 
        obj.imgdata = obj.OPT_data_controller.perform_reconstruction(false);
     end
     
     u = double(obj.imgdata);          
     
     sgm = obj.do_OPT_Mouse_Lung_Segmentation(false);     
     %
     % this is quite general morphometry set - use regionprops3?
     L = bwlabeln(sgm);
     stats = regionprops(L,'Centroid','Area','PixelList');
  
     datas = [];
     
     K = (4/3*pi)^(2/3)/(4*pi);
     hw = waitbar(0,'gathering 3D Morphometry parameters, please wait');
     for k=1:numel(stats)
         pl = stats(k).PixelList;
        x = pl(:,1);
        y = pl(:,2);
        z = pl(:,3);
        shp = alphaShape(x,y,z);
        A = surfaceArea(shp);
        V = volume(shp);
        sphericity = K*A/(V^(2/3));
        sphericity = 1/sphericity; % :)
        N = size(pl,1);
        if isnan(sphericity) || N < 26
            sphericity = 1; % too small, - like a ball
        end
        %
        if 0==A, A=1; end
        if 0==V, V=1; end
        %
        Rc = stats(k).Centroid;
        if size(x,1)==1
            Rg = (0.5*(1+sqrt(2))/2)/sqrt(5/3);
        else
            Rg = GyrationRadius(x,y,z,Rc(1),Rc(2),Rc(3));
        end
        % 
        rec = {N V A sphericity Rg};
        datas = [ datas; rec];
        if ~isempty(hw), waitbar(k/numel(stats),hw); drawnow, end   
     end
     if ~isempty(hw), delete(hw), drawnow; end
     % this is quite general morphometry set - rearrange as separate proc?
     
     %
     captions = {'Npix','Volume','Area','sphericity','Rg'};
     %
     %
     [sX,sY,sZ] = size(u);
     icyvol = zeros(sX,sY,1,sZ,1,'uint16');          
     icyvol(:,:,1,:,1) = uint16(u);
     icyvol(:,:,2,:,1) = uint16(L);
     %
     fig = icyvol;            
          
end