function ic_OPTtools_reconstruct(SRC_DIR_NAME,DST_DIR_NAME,swap)

    if ~( (isdir(SRC_DIR_NAME) || exist(SRC_DIR_NAME,'file')) && isdir(DST_DIR_NAME) )
        disp('input parameters are not valid directory names, can not continue');
        return;
    end

    addpath_ALYtools;
    
    bfCheckJavaMemory;
    bfCheckJavaPath;
    bfUpgradeCheck;
       
    dc = ic_OPTtools_data_controller([]);
    if strcmpi(swap,'y') % then swap dimensions
        dc.swap_XY_dimensions = 'Y';
    else
        dc.swap_XY_dimensions = 'N';
    end    
    
    % load settings filename in SRC directory if there is one
    if isdir(SRC_DIR_NAME)
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
    end
    %            
    if 2==exist(SRC_DIR_NAME) % SRC_DIR_NAME is single OME.tiff file not directory
        try
                        fname = SRC_DIR_NAME;
                        infostring = dc.Set_Src_Single(fname,false);
                        if isempty(infostring)
                            infostring = dc.Set_Src_FLIM(fname,'sum',false);
                        end
                        if ~isempty(infostring) 
                            dc.Z_range = []; % no selection
                            if ~isempty(dc.delays) %FLIM
                                dc.perform_reconstruction_FLIM;
                            else
                                if strcmp(dc.Reconstruction_Largo,'ON')
                                    dc.perform_reconstruction_Largo;
                                else
                                    dc.volm = dc.perform_reconstruction(false);
                                end
                            end
                            %
                            % save volume on disk
                            s = strsplit(fname,filesep);
                            iName = char(s(numel(s)));
                            L = length(iName);
                            savefilename = [iName(1:L-9) '_VOLUME.OME.tiff'];
                            if isempty(dc.delays) % non-FLIM
                                dc.save_volume([DST_DIR_NAME filesep savefilename],false); % silent
                            else
                                dc.save_volm_FLIM([DST_DIR_NAME filesep savefilename],false); % silent
                            end
                            %
                            clear('dc');
                            if isdeployed
                                exit;
                            else
                                return;                                    
                            end                               
                        end
        catch
            disp(['error while trying to reconstruct file ' SRC_DIR_NAME]);
        end
    elseif 7==exist(SRC_DIR_NAME) % presume that SRC_DIR_NAME is single projections stack 
        try
                            pth = SRC_DIR_NAME
                            if isempty(dc.imstack_get_delays(pth))                    
                                infostring = dc.imstack_Set_Src_Single(pth,false);
                            else
                                infostring = dc.imstack_Set_Src_Single_FLIM(pth,'sum',false);
                            end
                            if ~isempty(infostring)
                                dc.Z_range = []; % no selection
                                if ~isempty(dc.delays) %FLIM
                                    dc.perform_reconstruction_FLIM;
                                else
                                    if strcmp(dc.Reconstruction_Largo,'ON')
                                        dc.perform_reconstruction_Largo;
                                    else
                                        dc.volm = dc.perform_reconstruction(false);
                                    end
                                end
                                %
                                % save volume on disk
                                s = strsplit(pth,filesep);
                                iName = char(s(numel(s)));
                                savefilename = [iName '_VOLUME.OME.tiff'];
                                if isempty(dc.delays) % non-FLIM
                                    dc.save_volume([DST_DIR_NAME filesep savefilename],false); % silent
                                else
                                    dc.save_volm_FLIM([DST_DIR_NAME filesep savefilename],false); % silent
                                end
                                %
                                clear('dc');
                                if isdeployed
                                    exit;
                                else
                                    return;                                    
                                end                               
                            end
        catch
            disp(['error while trying to reconstruct stack ' SRC_DIR_NAME]);
        end
    end
    %    
    % SRC_DIR_NAME is the directory containing or OME.tiff images or projections-directories
    try
        dc.BatchSrcDirectory = SRC_DIR_NAME;
        dc.BatchDstDirectory = DST_DIR_NAME;    
        dc.run_batch([],false); % not opening progress bar
        clear('dc');
        %
        if isdeployed
            exit;
        end
    catch
         disp(['error while trying to reconstruct stack ' SRC_DIR_NAME]);        
    end
end

