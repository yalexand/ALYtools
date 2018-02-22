function ic_OPTtools_split(SRC_DIR_NAME,DST_DIR_NAME,N_str,swap)

    if ~isdir(DST_DIR_NAME)
        disp('error in DST directory spec');
        return;
    end
    %
    N = fix(str2num(N_str));
    if N>40 && N<2
        disp('bad split number spec, can not continue');
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

    infostring = [];
    if 2==exist(SRC_DIR_NAME) % SRC_DIR_NAME is single OME.tiff file not directory
        try
                        fname = SRC_DIR_NAME;
                            infostring = dc.Set_Src_Single(fname,false);
                        if isempty(infostring)
                            infostring = dc.Set_Src_FLIM(fname,'sum',false);
                        end
                        if isempty(infostring) 
                            if isdeployed
                                exit;
                            else
                                return;                                    
                            end                               
                        end
        catch
            disp(['error while trying to split file ' SRC_DIR_NAME]);
        end
    elseif 7==exist(SRC_DIR_NAME) % presume that SRC_DIR_NAME is single projections stack 
        try
                            pth = SRC_DIR_NAME
                            if isempty(dc.imstack_get_delays(pth))                    
                                infostring = dc.imstack_Set_Src_Single(pth,false);
                            else
                                infostring = dc.imstack_Set_Src_Single_FLIM(pth,'sum',false);
                            end
                            if  isempty(infostring)
                                if isdeployed
                                    exit;
                                else
                                    return;                                    
                                end                               
                            end
        catch
            disp(['error while trying to split stack ' SRC_DIR_NAME]);
        end
    end
    %    
    s = strsplit(infostring,filesep);
    name = char(s(numel(s)));
    name = strrep(name,'.OME.tiff','');    
    %
    try
    dc.split_original(DST_DIR_NAME,name,N); % current angles will be used    
    clear('dc');
        if isdeployed
            exit;
        end
    catch
         disp(['error while trying to split data ' SRC_DIR_NAME]);        
    end            
end

