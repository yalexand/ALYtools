function simple_reconstruct_headless()

% USAGE:
% matlab -nodisplay -nosplash -nodesktop -r "run('path/to/your/script.m');exit;"
% matlab -nodisplay -nosplash -nodesktop -r "run('path/to/your/script.m');exit;" | tail -n + 11

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
    
    if ~isdir(SRC_DIR_NAME) || ~isdir(DST_DIR_NAME)
        disp('error - wrong src or dst directory spec');
        return
    end

    addpath_OMEkit;
    
    bfCheckJavaMemory;
    bfCheckJavaPath;
    bfUpgradeCheck;

    dc = ic_OPTtools_data_controller([]);

    % SETTINGS - this block is not necessary, 
    % if one simply relies on what is in the settings file
    dc.Prefiltering_Size = 'None';
    dc.swap_XY_dimensions = 'AUTO';
    dc.registration_method = 'None';
    dc.imstack_filename_convention_for_angle = 'C1';
    %
    dc.downsampling = 2;
    dc.angle_downsampling = 1;
    %
    dc.FBP_interp = 'linear';
    dc.FBP_filter = 'Ram-Lak';
    dc.FBP_fscaling = 1;
    %
    dc.Reconstruction_Method = 'FBP';
    dc.Reconstruction_GPU = 'OFF';
    dc.Reconstruction_Largo = 'OFF';
    %
    dc.TwIST_TAU = 0.0008;
    dc.TwIST_LAMBDA = 0.0001;
    dc.TwIST_ALPHA = 0;
    dc.TwIST_BETA = 0;
    dc.TwIST_STOPCRITERION = 1;
    dc.TwIST_TOLERANCEA = 0.0001;
    dc.TwIST_TOLERANCED = 0.0001;
    dc.TwIST_DEBIAS = 0;
    dc.TwIST_MAXITERA = 10000;
    dc.TwIST_MAXITERD = 200;
    dc.TwIST_MINITERA = 5;
    dc.TwIST_MINITERD = 5;
    dc.TwIST_INITIALIZATION = 0;
    dc.TwIST_MONOTONE = 1;
    dc.TwIST_SPARSE = 1;
    dc.TwIST_VERBOSE = 0;
    % SETTINGS    
                    
    for k=3:numel(strings)
            cur_name = char(strings{k});            
            SRC = [SRC_DIR_NAME filesep cur_name];
                dc.Set_Src_Single(SRC,false);    
                dc.volm = dc.perform_reconstruction(false);    
            timestamp = datestr(now,'yyyy-mm-dd_HH-MM-SS');
            DST = [DST_DIR_NAME filesep timestamp '_REC_' cur_name];
                dc.save_volume(DST,false);
            disp(SRC);
    end      
    
end

