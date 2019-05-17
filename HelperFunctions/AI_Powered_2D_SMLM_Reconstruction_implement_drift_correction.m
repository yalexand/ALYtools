%-------------------------------------------------------------------------%     
function XY_UPS_corr = AI_Powered_2D_SMLM_Reconstruction_implement_drift_correction(XY_UPS,F_UPS,DX_DY_DRIFT)
    N = size(XY_UPS,1);
    XY_UPS_corr = zeros(size(XY_UPS));
    for k=1:N
        XY_UPS_corr(k,1) = XY_UPS(k,1) - DX_DY_DRIFT(F_UPS(k),1);
        XY_UPS_corr(k,2) = XY_UPS(k,2) - DX_DY_DRIFT(F_UPS(k),2);
    end
end