function corrmap = correlation_map( u1,u2,W )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %gs = (1+sqrt(5))/2;
    %H = round(W/gs);
    H = W;
    corrmap = zeros(W,H);
        
    u1_ = u1(~isnan(u1));
    u2_ = u2(~isnan(u2));

    u1_scaled = map(u1_,1,W);
    u2_scaled = map(u2_,1,H);
        
    for k=1:size(u1_scaled)
        u1_coord = fix(u1_scaled(k));
        u2_coord = fix(u2_scaled(k));
        corrmap(u1_coord,u2_coord) = corrmap(u1_coord,u2_coord) + 1;        
    end
    
    corrmap = flipud(corrmap');
    
    corrmap = log10(corrmap*10);
    corrmap(isinf(corrmap))=0;




end

