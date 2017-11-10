function ret_mask = retain_biggest_object(mask)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
        m2 = bwlabel(mask);
        lab_max = max(m2(:));
        vol_max=0;
        ind_max = 1;
        for k=1:lab_max
            vol_k = sum(m2(:)==k);
            if vol_k>vol_max
                vol_max = vol_k;
                ind_max = k;
            end        
            %disp([vol_k k ind_max]);
        end
        m2(m2~=ind_max)=0;    
        m2(m2==ind_max)=1;       
    ret_mask = m2;
end

