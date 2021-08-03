function perform_t_dependent_Nuclei_ratio_FRET_headless(src_dir,dst_dir,extension,full_settings_file_name)

%{
IMPORTANT
this application operates FIJI TrackMate plugin
the xml settings file should include correct specification for FIJI scripts directory, for example:
 <FijiScriptsDirectory>C:\Users\yalexand\Fiji.app\scripts</FijiScriptsDirectory>
OTHERWISE IT WON'T WORK
%}

if ~isdeployed
    resolve_xlwrite_path_issue();
end

tic

dc = ALYtools_data_controller(true,[]);
dc.application = 't_dependent_Nuclei_ratio_FRET';

if isfile(full_settings_file_name)
    try
    dc.load_settings(full_settings_file_name);
    catch ex
        disp(ex.message());
        disp('wrong or corrupt settings file, cannot continue');
        return;
    end
else
    disp('wrong settings file path, cannot continue');
    return;
end

dirdata = dir([src_dir filesep '*.' extension]);
if ~isempty({dirdata.name})
    fnames = sort_nat({dirdata.name});
else
    disp('no data, cannot continue');
    return;    
end

%%%%%%%%
            % verify that enough memory is allocated
            bfCheckJavaMemory();
                                   
            % load both bioformats & OMERO
            autoloadBioFormats = 1;
            % load the Bio-Formats library into the MATLAB environment
            status = bfCheckJavaPath(autoloadBioFormats);
            assert(status, ['Missing Bio-Formats library. Either add loci_tools.jar '...
                'to the static Java path or add it to the Matlab path.']);
                        
            % initialize logging
            loci.common.DebugTools.enableLogging('INFO');
            java.lang.System.setProperty('javax.xml.transform.TransformerFactory', 'com.sun.org.apache.xalan.internal.xsltc.trax.TransformerFactoryImpl');            
            
            FijiScriptsDirectory = dc.FijiScriptsDirectory;
            if contains(version,'2018') || contains(version,'2019')
                if isfolder(FijiScriptsDirectory)
                    addpath(FijiScriptsDirectory);
                    Miji(false);
                end
            else
                if isdir(FijiScriptsDirectory)
                    addpath(FijiScriptsDirectory);
                    Miji(false);
                end                
            end           
%%%%%%%%

tic

for k = 1:numel(fnames)
    dc.current_filename = fnames{k};
    dc.open_image([src_dir filesep dc.current_filename]);
    %
    [~,~,~,fig] = dc.analyze_t_dependent_Nuclei_ratio_FRET;
    fig = dc.t_dependent_Nuclei_ratio_FRET_postprocess(fig,dst_dir);
    %
    ometiffsavename = [dst_dir filesep dc.current_filename '_analysis_output.OME.tiff'];
    bfsave(fig,ometiffsavename,'Compression','LZW','BigTiff', true,'dimensionOrder','XYCTZ');
end

disp(['execution time ' num2str(toc/60) ' min']);

end


%-------------------------------------------------------------%
function resolve_xlwrite_path_issue()
            jPath = javaclasspath('-dynamic');
            % first check it isn't already in the dynamic path
            WriteXLPath = false;
            for i = 1:length(jPath)
                if strfind(jPath{i},'jxl.jar');
                    WriteXLPath = true;
                    break;
                end
            end
                
            if ~WriteXLPath
                path = which('jxl.jar');
                if isempty(path)
                    path = fullfile(fileparts(mfilename('fullpath')), 'jxl.jar');
                end
                if ~isempty(path) && exist(path, 'file') == 2
                    javaaddpath(path);
                else 
                     assert('Cannot automatically locate an jxl JAR file');
                end
            end
            
             % first check it isn't already in the dynamic path
            WriteXLPath = false;
            for i = 1:length(jPath)
                if strfind(jPath{i},'MXL.jar');
                    WriteXLPath = true;
                    break;
                end
            end
                
            if ~WriteXLPath
                path = which('MXL.jar');
                if isempty(path)
                    path = fullfile(fileparts(mfilename('fullpath')), 'MXL.jar');
                end
                if ~isempty(path) && exist(path, 'file') == 2
                    javaaddpath(path);
                else 
                     assert('Cannot automatically locate an MXL JAR file');
                end
            end               
end