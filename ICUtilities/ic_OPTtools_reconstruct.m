function ic_OPTtools_reconstruct(SRC_DIR_NAME,DST_DIR_NAME,swap)

    if ~( (isdir(SRC_DIR_NAME) || 2==exist(SRC_DIR_NAME)) && isdir(DST_DIR_NAME) )
        disp('input parameters are not valid directory names, can not continue');
        return;
    end

            addpath_ALYtools;
                
            % verify that enough memory is allocated
            bfCheckJavaMemory();

            % load bioformats 
            autoloadBioFormats = 1;
            % load the Bio-Formats library into the MATLAB environment
            status = bfCheckJavaPath(autoloadBioFormats);
            assert(status, ['Missing Bio-Formats library. Either add loci_tools.jar '...
                    'to the static Java path or add it to the Matlab path.']);

            bfCheckJavaPath;
            bfUpgradeCheck;    

            % initialize logging
            loci.common.DebugTools.enableLogging('INFO');
            java.lang.System.setProperty('javax.xml.transform.TransformerFactory', 'com.sun.org.apache.xalan.internal.xsltc.trax.TransformerFactoryImpl');   
       
    dc = ic_OPTtools_data_controller(true,[]); %headless
    dc.run_headless = true;
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
                        tic
                        infostring = dc.Set_Src_Single(fname,false);
                        disp(infostring);
                        disp(['loading time = ' num2str(toc)]);
                        if isempty(infostring)
                            tic
                            infostring = dc.Set_Src_FLIM(fname,'sum',false);
                            disp(infostring);
                            disp(['loading time = ' num2str(toc)]);                            
                        end
                        if ~isempty(infostring) 
                            tic
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
                            disp(['reconstruction time = ' num2str(toc)]);
                            %
                            % save volume on disk
                            tic
                            s = strsplit(fname,filesep);
                            iName = char(s(numel(s)));
                            L = length(iName);
                            savefilename = [iName(1:L-9) '_VOLUME.OME.tiff'];
                            if isempty(dc.delays) % non-FLIM
                                dc.save_volume([DST_DIR_NAME filesep savefilename],false); % silent
                            else
                                dc.save_volm_FLIM([DST_DIR_NAME filesep savefilename],false); % silent
                            end
                            disp([DST_DIR_NAME filesep savefilename]);
                            disp(['saving time = ' num2str(toc)]);
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
                            pth = SRC_DIR_NAME;
                            if isempty(dc.imstack_get_delays(pth))     
                                tic
                                infostring = dc.imstack_Set_Src_Single(pth,false);
                                disp(infostring);
                                disp(['loading time = ' num2str(toc)]);                                
                            else
                                tic
                                infostring = dc.imstack_Set_Src_Single_FLIM(pth,'sum',false);
                                disp(infostring);
                                disp(['loading time = ' num2str(toc)]);                                
                            end
                            if ~isempty(infostring)
                                tic
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
                                disp(['reconstruction time = ' num2str(toc)]);
                                %
                                % save volume on disk
                                tic
                                s = strsplit(pth,filesep);
                                iName = char(s(numel(s)));
                                savefilename = [iName '_VOLUME.OME.tiff'];
                                if isempty(dc.delays) % non-FLIM
                                    dc.save_volume([DST_DIR_NAME filesep savefilename],false); % silent
                                else
                                    dc.save_volm_FLIM([DST_DIR_NAME filesep savefilename],false); % silent
                                end
                                disp([DST_DIR_NAME filesep savefilename]);
                                disp(['saving time = ' num2str(toc)]);
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

