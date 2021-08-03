function load_segmentation(obj, ~, ~)
       directoryname = uigetdir(obj.DefaultDirectory);
    if isnumeric(directoryname) || isempty(directoryname), return, end;
    %
    Extension = 'tif';
    files = dir([directoryname filesep '*.' Extension]);
    num_seg_files = numel(files);                        
    seg_names_list = {files.name}; 
    
            switch obj.application

                case 'Fungus Dependent Granule Release'                    
                case 'CIDR'                    
                case 'TTO'                    
                case 'PR'                    
                case 'HL1'
                case 'NucCyt'
                case 'MPHG'
                case 'Sparks'
                case 'Experimental'
                case {'per_image_TCSPC_FLIM','per_image_TCSPC_FLIM_PHASOR'}
                    if isempty(obj.M_imgdata), return, end;
                    try
                        for k=1:numel(obj.M_imgdata)
                            fname = char(obj.M_filenames{k});
                            fname = strrep(fname,'.sdt','');
                            fname = strrep(fname,'.ome.tiff','');
                            fname = strrep(fname,'.OME.tiff','');
                            
                            for m=1:num_seg_files
                                segname_vanilla = seg_names_list{m};
                                segname = segname_vanilla;
                                segname = strrep(segname,'.tif','');
                                segname = strrep(segname,'.tiff','');
                                segname = strrep(segname,'.bmp','');
                                segname = strrep(segname,'_segmentation','');
                                segname = strrep(segname,' segmentation','');
                                if strcmp(fname,segname)
                                    sgm = imread([directoryname filesep segname_vanilla]);
                                    obj.M_sgm{k} = sgm;
                                    [directoryname filesep segname_vanilla]
                                end            
                            end
                        end
                    catch
                    end
            end
end    