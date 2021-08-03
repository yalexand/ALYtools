function save_segmentation(obj, ~, ~)
    %
    directoryname = uigetdir(obj.DefaultDirectory);
    if isnumeric(directoryname) || isempty(directoryname), return, end;
    %
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
                case 'per_image_TCSPC_FLIM'
                    if isempty(obj.M_imgdata), return, end;
                    try
                        sgm = obj.do_per_image_TCSPC_FLIM_Segmentation(false);
                        hw = waitbar(0,'..saving segmentations..');
                        for k=1:numel(obj.M_imgdata)
                            fname = char(obj.M_filenames{k});
                            fname = strrep(fname,'.sdt','');
                            fname = strrep(fname,'.ome.tiff','');
                            fname = strrep(fname,'.OME.tiff','');
                            fname = [ directoryname filesep fname '_segmentation.tif'];
                            mask = sgm{k};
                            imwrite(mask,fname,'tif');
                            if ~isempty(hw), waitbar(k/numel(obj.M_imgdata),hw); drawnow, end;
                        end
                        if ~isempty(hw), delete(hw), drawnow; end;
                    catch
                    end
            end
end    