function AI_Powered_2d_SMLM_reconstruct(SRC_FILE_NAME,DST_DIR_NAME,settings_fname,src_filename_must_have_token)

    if ~( isfile(SRC_FILE_NAME) && isfolder(DST_DIR_NAME) && isfile(settings_fname) && ... 
            contains(SRC_FILE_NAME,src_filename_must_have_token) )
        disp('input parameters are not valid directory names, can not continue');
        return;
    end
            if ~isdeployed, addpath_ALYtools, end
                
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
       
    dc = ALYtools_data_controller(true,[]);
    dc.application = 'AI_Powered_2D_SMLM_Reconstruction';
    % NB - all settings should be properly specified in xml file
    
    try
        dc.load_settings(settings_fname)
    catch
        disp('improper settings file spec, can not continue');
        return
    end
    %            
        try
            t_start = tic;
            dc.open_image(SRC_FILE_NAME);
            [data, caption, ~, fig] = dc.analyze_AI_Powered_2D_SMLM_Reconstruction;            
                s = strsplit(SRC_FILE_NAME,filesep);
                cur_file_name = char(s(length(s)));            
            cur_dir_name = strrep(cur_file_name,'.OME','');
            cur_dir_name = strrep(cur_dir_name,'.ome','');
            cur_dir_name = strrep(cur_dir_name,'.tiff','');
            cur_dir_name = strrep(cur_dir_name,'.tif','');
            cur_dir_created = mkdir(DST_DIR_NAME,cur_dir_name);
                if ~cur_dir_created
                    disp(['can not create output directory for a file ' cur_file_name ', leaving..']);
                end
            cur_base_name = [DST_DIR_NAME filesep cur_dir_name filesep cur_dir_name];
            % localisations
            save([cur_base_name '.mat'],'caption','data');
            cell2csv([cur_base_name '.csv'],[caption; data]);
            % image
            imagefilesavename = [cur_base_name '_2D.ome.tiff'];
            bfsave(squeeze(fig(:,:,2,1,1)),imagefilesavename,'dimensionOrder','XYCZT','BigTiff',true,'Compression','LZW');
            disp(imagefilesavename);
            disp(toc(t_start));
            
        catch
            disp(['error while trying to reconstruct file ' SRC_FILE_NAME]);
        end
end
