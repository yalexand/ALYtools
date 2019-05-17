%-------------------------------------------------------------------------%     
function [XY_corr,F_corr] = AI_Powered_2D_SMLM_Reconstruction_merge_localisations_in_frames(XY,F,block_size)

    N = size(XY,1);

    N_blocks = ceil(N/block_size);

    XY_corr = [];
    F_corr = [];
    for k=1:N_blocks
        i1 = block_size*(k-1)+1;
        i2 = i1+block_size-1;
            if i2 > N, i2 = N; end
        [XY_corr_k,F_corr_k] = process_block(XY(i1:i2,:),F(i1:i2));
        XY_corr = [XY_corr; XY_corr_k];
        F_corr = [F_corr; F_corr_k];
    end

end
%-------------------------------------------------------------------------%     
function [XY_corr,F_corr] = process_block(XY,F)
    [XY_corr,index,~] = unique(XY,'rows');    
    F_corr = F(index);
end