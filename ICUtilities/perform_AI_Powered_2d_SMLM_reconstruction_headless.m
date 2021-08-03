function perform_AI_Powered_2d_SMLM_reconstruction_headless()

    def_filename = 'headless_def.txt';

    try
        thisdir = pwd;
        fulldeffilename = [thisdir filesep def_filename];
    catch
        disp('error while trying to read "headless_def.txt"');
    end
    
    fileID = fopen(fulldeffilename);
    C = textscan(fileID,'%s');
        fclose(fileID);
    strings = C{1};    
    %
    % CONVENTION: FIRST TWO LINE IN THE DEF FILE ARE SRC DST FOLDERS    
    SRC_DIR_NAME = strings{1};
    DST_DIR_NAME = strings{2};    
    
    if ~isfolder(SRC_DIR_NAME) || ~isfolder(DST_DIR_NAME)
        disp('error - wrong src or dst directory spec');
        return
    end

    addpath_ALYtools;
    
    bfCheckJavaMemory;
    bfCheckJavaPath;
    bfUpgradeCheck;

    dc = ALYtools_data_controller(true,[]);
    dc.application = 'AI_Powered_2D_SMLM_Reconstruction';
    % NB - all settings should be properly specified in xml file
                   
    t_start = tic;
    for k=3:numel(strings)
            cur_name = char(strings{k});
            dc.open_image([SRC_DIR_NAME filesep cur_name]);    
            [data, caption, ~, fig] = dc.analyze_AI_Powered_2D_SMLM_Reconstruction;
            cur_dir_name = strrep(cur_name,'.OME','');
            cur_dir_name = strrep(cur_dir_name,'.ome','');
            cur_dir_name = strrep(cur_dir_name,'.tiff','');
            cur_dir_name = strrep(cur_dir_name,'.tif','');
            cur_dir_created = mkdir(DST_DIR_NAME,cur_dir_name);
                if ~cur_dir_created, continue, end
            cur_base_name = [DST_DIR_NAME filesep cur_dir_name filesep cur_dir_name];
            % localisations
            save([cur_base_name '.mat'],'caption','data');
            cell2csv([cur_base_name '.csv'],[caption; data]);
            % image
            imagefilesavename = [cur_base_name '_2D.ome.tiff'];
            bfsave(fig,imagefilesavename,'dimensionOrder','XYCZT','BigTiff',true,'Compression','LZW');
            disp(imagefilesavename);
            disp(toc(t_start));
    end      
    
end

