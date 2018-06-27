function [datas, captions, table_names, fig] = analyze_NucleiTimeStack(obj,~,~) 

     datas = [];
     captions = [];
     table_names = 'default';
     fig = [];
     
     fig = obj.do_NucleiTimeStack_Segmentation(false);
     
end
