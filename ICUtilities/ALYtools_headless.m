% example function showing how to run analysis headless
% NB - results appear by default in the "pwd" directory
function ALYtools_headless(SRC_DIR,application,ext)

    if ~isdir(SRC_DIR)
        disp('input parameters are not valit directory names, can not continue');
    end
    
    dirdata = dir( fullfile(SRC_DIR,['*.' ext]));

    addpath_ALYtools;
    
    bfCheckJavaMemory;
    bfCheckJavaPath;
    bfUpgradeCheck;

    dc = ALYtools_data_controller([]);
    dc.application = application;
                    
    for k=1:numel(dirdata)
            cur_name = dirdata(k).name;
            dc.open_image([SRC_DIR filesep cur_name]);
            dc.analyze_current;
    end      

end

