function [ mapfile_name,memmap ] = initialize_memmap( size_arr, repeats, var_name, datatype, varargin )

            ip = inputParser;
            ip.addOptional('ini_data', [], @isnumeric);
            ip.parse(varargin{:});

            mapfile_name = global_tempname;                        
            mapfile = fopen(mapfile_name,'w');
            ini_data = zeros(1,prod(size_arr)*repeats,datatype);
            fwrite(mapfile,ini_data,datatype);
            fclose(mapfile);
            memmap = memmapfile(mapfile_name,'Writable',true,'Repeat',repeats,'Format',{datatype, size_arr, var_name},'Offset',0);
            
            if ~isempty(ip.Results.ini_data)
                memmap.Data.(var_name) = ip.Results.ini_data;
            end

end

