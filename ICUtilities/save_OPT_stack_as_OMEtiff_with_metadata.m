function speed = save_OPT_stack_as_OMEtiff_with_metadata(folder, ometiffilename, varargin )

addpath_OMEkit;

            extension  = 'tif';
            
            files = dir([folder filesep '*.' extension]);
            num_files = length(files);
            if 0 ~= num_files
                file_names = cell(1,num_files);
                for k = 1:num_files
                    file_names{k} = char(files(k).name);
                end
            else
                return, 
            end;               

            if isempty(file_names) || 0 == numel(file_names), return, end;
            
            num_files = numel(file_names);
            %
            acc = zeros(1,num_files);
            for i = 1 : num_files                                   
                try I = imread([folder filesep file_names{i}],extension); catch err, msgbox(err.mesasge), return, end;                
                % telapsed = add_plane_to_OMEtiff_with_metadata(I, i, num_files, folder, ometiffilename);                
                
                % telapsed = add_plane_to_OMEtiff_with_metadata(I, i, num_files, folder, ometiffilename, varargin{:} );
                
                z = i;
                c = 1;
                t = 1;
                sizeZ = num_files;
                sizeC = 1;
                sizeT = 1;
                telapsed = add_plane_to_OMEtiff_with_metadata_ZCT(I, [z c t], [sizeZ sizeC sizeT], folder, ometiffilename, varargin{:} );
                
                acc(i)=telapsed;                                                   
            end
            speed = mean(acc);

