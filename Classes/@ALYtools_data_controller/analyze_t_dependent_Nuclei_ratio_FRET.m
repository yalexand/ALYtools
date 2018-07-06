function [datas, captions, table_names, fig] = analyze_t_dependent_Nuclei_ratio_FRET(obj,~,~) 

     datas = [];
     captions = [];
     table_names = 'default';
     fig = [];
     %
     fig = obj.do_t_dependent_Nuclei_ratio_FRET_Segmentation(false);
end
