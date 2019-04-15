function perform_image_tiling(in_src_dir,in_dst_dir,in_Extension,in_ncols,in_nrows,ovlp_x,ovlp_y,mode,QT)

if ~isdeployed
    resolve_xlwrite_path_issue();
end

tic

dc = ALYtools_data_controller(true,[]);
dc.problem = 'Image_Tiling';

Extension = in_Extension; %'ome.tif';

dc.ImageTiling_mode = mode;

if isnumeric(in_ncols)
    dc.ImageTiling_Ncols = in_ncols;
else
    dc.ImageTiling_Ncols = str2num(in_ncols);
end
    if isnumeric(in_nrows)
        dc.ImageTiling_Nrows = in_nrows;
    else
        dc.ImageTiling_Nrows = str2num(in_nrows);
    end
        if isnumeric(ovlp_x)
            dc.ImageTiling_Ovlp_X = ovlp_x;
        else
            dc.ImageTiling_Ovlp_X = str2num(ovlp_x);
        end
            if isnumeric(ovlp_y)
                dc.ImageTiling_Ovlp_Y = ovlp_y;
            else
                dc.ImageTiling_Ovlp_Y = str2num(ovlp_y);
            end
                if isnumeric(QT)
                    dc.ImageTiling_QT = QT;
                else
                    dc.ImageTiling_QT = str2num(QT);
                end
       
dirdata = dir([in_src_dir filesep '*.' Extension]);
fnames = {dirdata.name};
fnames = sort_nat(fnames);
dc.load_multiple(fnames,in_src_dir,false);

[datas, captions, table_names, fig] = dc.analyze_ImageTiling(false);

str = strsplit(in_src_dir,filesep);
xlsname = [in_dst_dir filesep char(str(numel(str))) '.xls'];
xlwrite(xlsname,[captions; datas]);
filesavename = [in_dst_dir filesep char(str(numel(str))) '.ome.tif'];
bfsave(fig,filesavename);

disp(['execution time ' num2str(toc/60) ' min']);

end


%-------------------------------------------------------------%
function resolve_xlwrite_path_issue()
            jPath = javaclasspath('-dynamic');
            % first check it isn't already in the dynamic path
            WriteXLPath = false;
            for i = 1:length(jPath)
                if strfind(jPath{i},'jxl.jar');
                    WriteXLPath = true;
                    break;
                end
            end
                
            if ~WriteXLPath
                path = which('jxl.jar');
                if isempty(path)
                    path = fullfile(fileparts(mfilename('fullpath')), 'jxl.jar');
                end
                if ~isempty(path) && exist(path, 'file') == 2
                    javaaddpath(path);
                else 
                     assert('Cannot automatically locate an jxl JAR file');
                end
            end
            
             % first check it isn't already in the dynamic path
            WriteXLPath = false;
            for i = 1:length(jPath)
                if strfind(jPath{i},'MXL.jar');
                    WriteXLPath = true;
                    break;
                end
            end
                
            if ~WriteXLPath
                path = which('MXL.jar');
                if isempty(path)
                    path = fullfile(fileparts(mfilename('fullpath')), 'MXL.jar');
                end
                if ~isempty(path) && exist(path, 'file') == 2
                    javaaddpath(path);
                else 
                     assert('Cannot automatically locate an MXL JAR file');
                end
            end           
end