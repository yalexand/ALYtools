function [datas, captions, table_names, fig] = analyze_DarkNuclei(obj,~,~) 

     datas = [];
     captions = [];
     table_names = 'default';
     fig = [];
     
     sgm = obj.do_DarkNuclei_Segmentation(false);
     %
     % and that is it
     fig = sgm;
              
     fname = obj.current_filename;
     fname = strrep(fname,'.sdt','');
     fname = strrep(fname,'.ome.tiff','');
     fname = strrep(fname,'.OME.tiff','');
     fname = [ obj.RootDirectory filesep fname '_segmentation.tif'];
     mask = squeeze(sgm(:,:,2,1,1));
     imwrite(mask,fname,'tif');
     
end
