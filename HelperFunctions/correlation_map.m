function corrmap = correlation_map( u1,u2,W )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %gs = (1+sqrt(5))/2;
    %H = round(W/gs);
    H = W;
    corrmap = zeros(W,H);
        
    u1_scaled = map(u1,1,W);
    u2_scaled = map(u2,1,H);
    
    if 0~=sum(isnan(u1_scaled))
        u1_scaled = ones(size(u1_scaled));
    end
    if 0~=sum(isnan(u2_scaled))
        u2_scaled = ones(size(u2_scaled));
    end
    
    sizeX=size(u1);
    for x=1:sizeX
        u1_coord = fix(u1_scaled(x));
        u2_coord = fix(u2_scaled(x));
        corrmap(u1_coord,u2_coord) = corrmap(u1_coord,u2_coord) + 1;        
    end
    
    corrmap = flipud(corrmap');
    
    corrmap = log10(corrmap*10);
    corrmap(isinf(corrmap))=0;
    
end

