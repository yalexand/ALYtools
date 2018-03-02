   function out = parse_string_for_attribute_value(string,specs)

    out = [];

    NN = length(string);

        for kk = 1 : numel(specs)

            tok = specs{kk};
            toklen = length(tok); 

            startind = strfind(string,tok) + toklen;

            %fr000Rot325_0000.tif
            vallen = 0;
            for mm = startind : NN    
                if ~isempty(regexp(string(mm),'\d','once'))
                    vallen = vallen + 1;
                else
                    break, 
                end;            
            end

            endind = startind + vallen - 1; 

            try
                val = str2num(string(startind:endind));
            catch
            end;

            out{kk} = [];
            if exist('val','var') && ~isempty(val)
                out{kk}.attribute = tok;
                out{kk}.value = val;
                clear('val');
            end

        end

    end