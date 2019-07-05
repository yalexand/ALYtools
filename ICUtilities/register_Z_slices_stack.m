function register_Z_slices_stack(src_dir,dst_dir,tmp_dir,Extension,n_begin_str,n_end_str,direction_flag)

    if ~( isfolder(src_dir) && isfolder(dst_dir) && isfolder(tmp_dir) ) || ~ischar(Extension) ...
        || ~ischar(n_begin_str) || ~ischar(n_end_str) ...
        || ~ismember(direction_flag,{'refer_to_first_image','refer_to_last_image'})
        disp('bad input parameters, can not continue');
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
               
% Extension = 'ome.tif';

dirdata = dir([src_dir filesep '*.' Extension]);
if ~isempty({dirdata.name})
    fnames = {dirdata.name};
    fnames = sort_nat(fnames);
else
    disp('no suitable images found, nothing to do, exiting..');
    return;
end

n_begin = int64(str2num(n_begin_str));
n_end = int64(str2num(n_end_str));

if isempty(n_begin) || isempty(n_end) || n_begin<1 || n_end>numel(fnames) || n_begin>=n_end
    disp('bad slices range specs, nothing to do, exiting..');
    return;
end

fnames = fnames(n_begin:n_end);

if strcmp(direction_flag,'refer_to_last_image')    
    fnames = fliplr(fnames);
end

fac = 2; % downsampling

transforms = cell(1,numel(fnames)-1);

% THIS CREATRES TRANSFORMS FOR 10X UNDER-SAMPLED IMAGES

[optimizer, metric] = imregconfig('monomodal');
% [optimizer, metric] = imregconfig('multimodal');


% CAN COMMENT THIS AS IMAGES ARE SAVED IN TMP...
parfor k=1:numel(fnames)-1   
    tstart = tic;
    
    fname1 = char(fnames(k))
    [~,~,u1] = bfopen_v([src_dir filesep fname1]);
    u1 = single(u1);
    
    fname2 = char(fnames(k+1))
    [~,~,u2] = bfopen_v([src_dir filesep fname2]);
    u2 = single(u2);                
    
    u1 = imresize(u1,1/fac);
    u2 = imresize(u2,1/fac);
    [gx,gy] = gsderiv(u1,3,1);
        fixed = sqrt(gx.*gx + gy.*gy) + u1/3;        
    [gx,gy] = gsderiv(u2,3,1);
        moving = sqrt(gx.*gx + gy.*gy) + u2/3;
        
    transforms{k} = imregtform(moving,fixed,'rigid',optimizer,metric); 
                           
    toc(tstart)
end

% % perform crude registration of full size images
fname = char(fnames(1))
[~,~,u1] = bfopen_v([src_dir filesep fname]);
    u1 = single(u1);
    v = zeros(size(u1,1),size(u1,2),1,1,1,'uint16');
    v(:,:,1,1,1) = uint16(map(u1,0,65535));        
    fname_reg = [tmp_dir filesep fname]
    bfsave(v,fname_reg,'BigTiff',true,'Compression','LZW','DimensionOrder','XYZCT');

parfor k=2:numel(fnames)   
tstart = tic;
    fname = char(fnames(k))
    [~,~,u_k] = bfopen_v([src_dir filesep fname]);
    u_k = single(u_k);
    %
    T = eye(3);
    for m=k-1:-1:1
        t_m = transforms{m};
        T = T*t_m.T;
        disp([k m]);
    end    
    T(3,1:2) = T(3,1:2)*fac; % !
    tform = affine2d(T);
    u_k_reg = imwarp(u_k,tform,'OutputView',imref2d(size(u1)));    
    
    v = zeros(size(u1,1),size(u1,2),1,1,1,'uint16');
    v(:,:,1,1,1) = uint16(map(u_k_reg,0,65535));
    %
    fname_reg = [tmp_dir filesep fname]
    bfsave(v,fname_reg,'BigTiff',true,'Compression','LZW','DimensionOrder','XYZCT'); 
toc(tstart)
end
% CAN COMMENT THIS AS IMAGES ARE SAVED IN TMP...

DRS = 512; % desired ROI size

N = numel(fnames);
%  GET IMAGES FROM THE ONLY RIGID DIRECTORY AND APPLY DEMONS TO THEM CONSECUTIVELY
for m=1:N-1
    tstart = tic;    
    if m==1
        fname = char(fnames(m));
        [~,~,prev] = bfopen_v([tmp_dir filesep fname]);
        %save
        fname = strrep(fname,['.' Extension],'');
        fullfname = [dst_dir filesep fname '_registered' '.ome.tif'];
        bfsave(prev,fullfname,'BigTiff',true,'Compression','LZW','DimensionOrder','XYZCT');         
    end
    [~,~,next] = bfopen_v([tmp_dir filesep char(fnames(m+1))]);
    prev = single(prev);
    next = single(next);

    %%%%%%%%%%%%%%%%%%%%%%%% "ROI by ROI"
    % moving = next;
    % static = prev;
    [Sx,Sy]=size(prev);

        u = prev + next;
        t=quantile(u(:),0.999);
        u(u>t)=t;
        u=map(u,0,1);

        down_size = round([Sx Sy]/DRS);
        sx = down_size(1);
        sy = down_size(2);
        u = imresize(u,[sx sy]);
        f=Sx/sx;

        t=0.005;
        u=(u>t);
        sum(u(:))

        [x,y] = find(u==1);
        x = round((x-0.5)*f);
        y = round((y-0.5)*f);
        r = round((1+0.08)*DRS/2);
        %
        next_reg = zeros(Sx,Sy);
        %
        store = cell(1,numel(x));
        % loop over ROIs
        for k=1:numel(x)    
            roix=max(1,x(k)-r):min(Sx,x(k)+r);
            roiy=max(1,y(k)-r):min(Sy,y(k)+r);
            moving = next(roix,roiy);
            static = prev(roix,roiy);    

              [gx,gy] = gsderiv(moving,3,1);
                moving_comp = sqrt(gx.*gx + gy.*gy) + moving/3;      
              [gx,gy] = gsderiv(static,3,1);
                static_comp = sqrt(gx.*gx + gy.*gy) + static/3;
              %  
              [dfield_gpu,~] = imregdemons(gpuArray(moving_comp),gpuArray(static_comp),[500 400 200],...
                  'AccumulatedFieldSmoothing',1.3);
              dfield = gather(dfield_gpu);            
            
            store{k} = imwarp(moving,dfield,'SmoothEdges', true);
            disp([k toc(tstart)]);
        end

        for k=1:numel(x)
            roix=max(1,x(k)-r):min(Sx,x(k)+r);
            roiy=max(1,y(k)-r):min(Sy,y(k)+r);
            next_reg(roix,roiy) = store{k};
        end        
    %%%%%%%%%%%%%%%%%%%%%%%% "ROI by ROI"
    
      % save next_reg on disk        
      fname = char(fnames(m+1));
      fname = strrep(fname,['.' Extension],'');
      fullfname = [dst_dir filesep fname '_registered' '.ome.tif'];
      bfsave(uint16(map(next_reg,0,65535)),fullfname,'BigTiff',true,'Compression','LZW','DimensionOrder','XYZCT');
            
      prev = next_reg; % corrected!
      toc(tstart)      
end

