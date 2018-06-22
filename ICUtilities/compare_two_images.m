function ret = compare_two_images(ffn1_,ffn2_)

    ret = [];

    ffn1 = ffn1_;
    ffn2 = ffn2_;

    %ffn1 = 'W:\SMLM_data\tubulin_red_blue_2d\A2_AF647_ox2_5K_1\ian_script_output\Sigma_PSF_Averaged_shifted_histograms.tif';
    %ffn2 = 'W:\SMLM_data\tubulin_red_blue_2d\A2_AF647_ox2_5K_1\phasor\A2_AF647_ox2_5K_1\Phasor_Averaged_shifted_ histograms.tif';

    [~,~,u1] = bfopen_v(ffn1);
    [~,~,u2] = bfopen_v(ffn2);

    mask = (u1>0) & (u2>0);

    u1 = u1(mask);
    u2 = u2(mask);
    u1 = u1(:);
    u2 = u2(:);

    ret.similarity = 1 -  mean(abs(u1-u2)./(u1+u2));

    % https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
    ret.pearsoncoef = corr(u1,u2);
    ret.reflective_correlation =  sum(u1.*u2)/sqrt(sum(u1.^2)*sum(u2.^2));

end

