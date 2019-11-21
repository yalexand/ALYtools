function perform_image_tiling(in_src_dir,in_dst_dir,in_Extension,in_ncols,in_nrows,ovlp_x,ovlp_y,mode,QT)

if ~isdeployed
    resolve_xlwrite_path_issue();
end

tic

dc = ALYtools_data_controller(true,[]);
dc.problem = 'Image_Tiling';

Extension = in_Extension; %'ome.tif';

dc.ImageTiling_mode = mode;
dc.send_analysis_output_to_Icy = false;
dc.save_analysis_output_as_OMEtiff = false;
dc.save_analysis_output_as_xls = false;

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
if ~isempty({dirdata.name})
    fnames = {dirdata.name};
    fnames = sort_nat(fnames);
    dc.load_multiple(fnames,in_src_dir,false);
else
    % try this arrangement: /output/fovdir1, fovdir2 etc.
    dirdata = dir([in_src_dir filesep 'output']);
    dirnames = setdiff({dirdata.name},{'.','..','.directory'});
    dirnames = sort_nat(dirnames);
     dc.imgdata = [];
     for k=1:numel(dirnames)
         curname = char(dirnames(k)); % FOV file name base is the same as folder name
         fullfname = [in_src_dir filesep 'output' filesep curname filesep curname Extension];
         [~,~,I] = bfopen_v(fullfname);
         dc.M_imgdata{k} = single(I);
         dc.M_filenames{k} = curname;        
     end    
end

[datas, captions, table_names, fig] = dc.analyze_ImageTiling(false);

str = strsplit(in_src_dir,filesep);
xlsname = [in_dst_dir filesep char(str(numel(str))) '.xls'];
xlwrite(xlsname,[captions; datas]);
filesavename = [in_dst_dir filesep char(str(numel(str))) '_stitched.ome.tif'];
bfsave(fig,filesavename,'BigTiff',true,'Compression','LZW');

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