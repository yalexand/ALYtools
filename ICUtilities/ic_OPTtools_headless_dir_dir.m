function ic_OPTtools_headless_dir_dir(SRC_DIR_NAME,DST_DIR_NAME)

    if ~isdir(SRC_DIR_NAME) || ~isdir(DST_DIR_NAME)
        disp('input parameters are not valid directory names, can not continue');
    end

    addpath_ALYtools;
    
    bfCheckJavaMemory;
    bfCheckJavaPath;
    bfUpgradeCheck;

    dc = ic_OPTtools_data_controller([]);
    
    dc.BatchSrcDirectory = SRC_DIR_NAME;
    dc.BatchDstDirectory = DST_DIR_NAME;
    
    dc.run_batch([],false); % not opening progress bar
            
end

