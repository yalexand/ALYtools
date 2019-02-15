function perform_image_tiling(in_src_dir,in_dst_dir,in_Extension,in_ncols,in_nrows,ovlp_x,ovlp_y,QT)

% perform_image_tiling('D:\Users\alexany\SMLM_2019\2019_02_04_Sunil_for_stitching\SLICE_2','D:\Users\alexany\SMLM_2019\2019_02_04_Sunil_for_stitching','ome.tif',5,1,0.5,0.5,0)

tic

dc = ALYtools_data_controller(true,[]);
dc.problem = 'Experimental';

Extension = in_Extension; %'ome.tif';

dc.ImageTiling_Ncols = in_ncols;
dc.ImageTiling_Nrows = in_nrows;
dc.ImageTiling_Ovlp_X = ovlp_x;
dc.ImageTiling_Ovlp_Y = ovlp_y; 
dc.ImageTiling_QT = QT; 

dirdata = dir([in_src_dir filesep '*.' Extension]);
fnames = {dirdata.name};
fnames = sort_nat(fnames);
dc.load_multiple(fnames,in_src_dir,false);

[datas, captions, table_names, fig] = dc.analyze_ImageTiling(false);

str = strsplit(in_src_dir,filesep);
xlsname = [in_dst_dir filesep char(str(numel(str))) '.xls'];
xlswrite( xlsname,[captions; datas]);
filesavename = [in_dst_dir filesep char(str(numel(str))) '.ome.tif'];
bfsave(fig,filesavename);

disp(['execution time ' num2str(toc/60) ' min']);

end

