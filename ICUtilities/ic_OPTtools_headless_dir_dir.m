    function ic_OPTtools_headless_dir_dir(SRC_DIR_NAME,DST_DIR_NAME)

    if ~isdir(SRC_DIR_NAME) || ~isdir(DST_DIR_NAME)
        disp('input parameters are not valid directory names, can not continue');
        return;
    end

    addpath_ALYtools;
    
    bfCheckJavaMemory;
    bfCheckJavaPath;
    bfUpgradeCheck;
       
    dc = ic_OPTtools_data_controller([]);
    
    % load settings filename in SRC directory if there is one
    old_dir = pwd;
    cd(SRC_DIR_NAME);
    xmlfiles = dir('*.xml');
    for k=1:numel(xmlfiles)
        fname = [SRC_DIR_NAME filesep xmlfiles(k).name];
        try
            if dc.load_settings(fname)
                disp(fname);
                break;
            end
        catch
        end
    end;
    cd(old_dir);
    %
        
    dc.BatchSrcDirectory = SRC_DIR_NAME;
    dc.BatchDstDirectory = DST_DIR_NAME;
    
    dc.run_batch([],false); % not opening progress bar
    
    if isdeployed
        exit
    end
            
end

