function ic_OPTtools_merge(SRC_DIR,DST_DIR,dim)

    if ~isdir(SRC_DIR) || ~isdir(DST_DIR)
        disp('error in directory spec');
        return;
    end
    %
    addpath_ALYtools;
    
    bfCheckJavaMemory;
    bfCheckJavaPath;
    bfUpgradeCheck;    

    try
        dim = str2num(dim);

        if ~isdir(SRC_DIR)
            disp('not  directory, can not continue');
            return;
        end
        % presume RAM is enough big
        volm = [];
        %
        files = dir([SRC_DIR filesep '*.OME.tiff']);
        num_files = length(files);
        names_list = [];
        if 0 ~= num_files
            names_list = cell(1,num_files);
                for k = 1:num_files
                    names_list{k} = char(files(k).name);
                end
        end
        names_list = sort_nat(names_list);
        %
        for k = 1:num_files
            [~,~,I] = bfopen_v([ SRC_DIR filesep char(names_list{k})]);
            disp(size(I));
            volm = cat(dim,volm,I);
        end
        %
        s = strsplit(names_list{1},'_ic_split_');
        common_name  = char(s{1});
        %
        [szX, szY, szZ] = size(volm);
        full_filename = [DST_DIR filesep common_name '_VOLUME' '.OME.tiff'];
        bfsave(reshape(volm,[szX,szY,1,1,szZ]),full_filename, ...
            'dimensionOrder','XYCTZ','Compression','LZW','BigTiff',true);
    catch
        disp('ic_OPTtools_merge - error occurred');
    end
end

