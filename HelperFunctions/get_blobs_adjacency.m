function A = get_blobs_adjacency(flab,patch_mask,verbose)

% TEST CODE
%
% N = 20;
% L = zeros(N);
% for k=1:3
%     for m=1:3
%         L(4*k,4*m)=1;
%     end
% end
% L = bwlabel(L);

% tic
% 
% %A = get_blobs_adjacency(L,[],true);
% 
% mask = zeros(N);
% mask(2:18,2:18)=1;
% A = get_blobs_adjacency(L,mask,true);
% 
% full(A)
%
% toc/60
%
% TEST CODE


MultiCore = feature('numCores')>=4;

if isempty(patch_mask)     

    D = bwdist(flab>0);        
    N = ceil(max(D(:)));
    Lg = flab;
    [sx,sy]=size(D);
        
    if verbose 
        h = waitbar(0,'Scanning Blobs');
    end    
    
    for i=1:N
        [xp,yp] = find(D>=i & D<(i+1));
        for k=1:numel(xp)
            xc=xp(k);
            yc=yp(k);
            if xc>1 && xc<sx && yc>1 && yc<sy && 0==Lg(xc,yc)
                rx = (xc-1):(xc+1);
                ry = (yc-1):(yc+1);
                vic = Lg(rx,ry);
                label = max(vic(:));
                Lg(xc,yc) = label;
            end
        end
            if verbose, try waitbar(i/N,h); catch, end, end        
    end
    
    if verbose, try close(h); catch, end, end
    
    %     v = zeros(size(D,1),size(D,2),2,1,1);
    %     v(:,:,1,1,1)=flab;
    %     v(:,:,2,1,1)=Lg;
    %     icy_imshow(v);        

    %https://stackoverflow.com/questions/22810540/matlab-identify-adjacient-regions-in-3d-image    
    
    Im = Lg(2:sx-1,2:sy-1);
    labels = unique(Im);
    A = zeros(numel(labels));
    N=numel(labels);

    if ~MultiCore
        if verbose 
            h = waitbar(0,'Analysing Adjacency');
        end
        for k=1:N    
            z = Im==k;
            z = imdilate(z,ones(3))-z;
            s = Im(z~=0);
            adjlabs=unique(s(:));    
            for m=1:numel(adjlabs)
                if k~=adjlabs(m)
                A(k,adjlabs(m))=1;
                end
            end 
            if verbose, try waitbar(k/N,h); catch, end, end
        end
        if verbose, try close(h); catch, end, end
    else % use parfor   
        adjlabs = cell(1,N);
        parfor k=1:N    
            z = Im==k;
            z = imdilate(z,ones(3))-z;
            s = Im(z~=0);
            adjlabs{k}=unique(s(:));
            k
        end
        for k=1:N
            adjlabs_k = adjlabs{k};
            for m=1:numel(adjlabs_k)
                if k~=adjlabs_k(m)
                A(k,adjlabs_k(m))=1;
                end
            end 
        end
    end % MultiCore

else % use mask
    try

    D = bwdist(flab>0);
    D = D.*patch_mask;
    flab = flab.*patch_mask;    
    N = ceil(max(D(:)));
    Lg = flab;
    [sx,sy]=size(D);        
        
    if verbose 
        h = waitbar(0,'Scanning Blobs');
    end    
        
    for i=1:N
        [xp,yp] = find(D>=i & D<(i+1));
        for k=1:numel(xp)
            xc=xp(k);
            yc=yp(k);
            if xc>1 && xc<sx && yc>1 && yc<sy && 0==Lg(xc,yc) && 0~=patch_mask(xc,yc)
                rx = (xc-1):(xc+1);
                ry = (yc-1):(yc+1);
                vic = Lg(rx,ry);
                label = max(vic(:));
                Lg(xc,yc) = label;
            end
        end
            if verbose, try waitbar(i/N,h); catch, end, end        
    end
    
    if verbose, try close(h); catch, end, end
        
%         v = zeros(size(D,1),size(D,2),2,1,1);
%         v(:,:,1,1,1)=flab;
%         v(:,:,2,1,1)=Lg;
%         icy_imshow(v);     
    
    Im = Lg + 1;
    labels = unique(Im);
    A = zeros(numel(labels));
    N=numel(labels);

    if ~MultiCore
        if verbose 
            h = waitbar(0,'Analysing Adjacency');
        end
        for k=1:N    
            z = Im==k;
            z = imdilate(z,ones(3))-z;
            s = Im(z~=0);
            adjlabs=unique(s(:));    
            for m=1:numel(adjlabs)
                if k~=adjlabs(m)
                A(k,adjlabs(m))=1;
                end
            end 
            if verbose, try waitbar(k/N,h); catch, end, end
        end
        if verbose, try close(h); catch, end, end
    else % use parfor   
        adjlabs = cell(1,N);
        parfor k=1:N    
            z = Im==k;
            z = imdilate(z,ones(3))-z;
            s = Im(z~=0);
            adjlabs{k}=unique(s(:));
            k
        end
        for k=1:N
            adjlabs_k = adjlabs{k};
            for m=1:numel(adjlabs_k)
                if k~=adjlabs_k(m)
                A(k,adjlabs_k(m))=1;
                end
            end 
        end
    end % MultiCore
        
    A = A(2:size(A,1),2:size(A,1));
    
    catch
    end
end
    
end
