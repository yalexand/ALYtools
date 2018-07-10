function corrmap = correlation_map( u1,u2,W )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %gs = (1+sqrt(5))/2;
    %H = round(W/gs);
    H = W;
    corrmap = zeros(W,H);

    u1_ = map(u1,1,W);
    u2_ = map(u2,1,H);
        
    mask = ~isnan(u1_) & ~isnan(u2_);
    
    u1_ = u1_(mask);
    u2_ = u2_(mask);
        
    for k=1:size(u1_)
            u1_coord = fix(u1_(k));
            u2_coord = fix(u2_(k));
            corrmap(u1_coord,u2_coord) = corrmap(u1_coord,u2_coord) + 1;                
    end
    
    corrmap = flipud(corrmap');
    
    corrmap = log10(corrmap*10);
    corrmap(isinf(corrmap))=0;

end

