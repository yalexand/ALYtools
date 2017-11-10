function ic_OPTtools_headless_dir_dir(SRC_DIR_NAME,DST_DIR_NAME)

    if ~isdir(SRC_DIR_NAME) || ~isdir(DST_DIR_NAME)
        disp('input parameters are not valit directory names, can not continue');
    end

    ext = 'OME.tiff';
    old_dir = pwd;
    cd(char(SRC_DIR_NAME));
    source_file_names = ls(['*.' ext]);
    cd(old_dir);

    addpath_OMEkit;
    
    bfCheckJavaMemory;
    bfCheckJavaPath;
    bfUpgradeCheck;

    dc = ic_OPTtools_data_controller([]);
                    
    for k=1:size(source_file_names,1)
            cur_name = source_file_names(k,:);            
            SRC = [SRC_DIR_NAME filesep cur_name];
                 dc.Set_Src_Single(SRC,false);    
                 dc.volm = dc.perform_reconstruction(false);    
             timestamp = datestr(now,'yyyy-mm-dd_HH-MM-SS');
            DST = [DST_DIR_NAME filesep timestamp '_REC_' cur_name];
                 dc.save_volume(DST,false);
            disp(SRC);
    end      
end

