%
function [range_out,psd_out,fpeak,confidence] = analyze_profile_periodicity_via_Fourier(profile,fL,fH,df) % frequency in 1/pix
% fL - spectrum start

fft_p = fft(profile);
psd = abs(fft_p).^2;

L = length(psd);
psd = psd(L-fL-1:L);
L = length(psd);

% intensity of high-frequency noise
I_hfn = sum(psd(L-fH:L));
% "filter out" hf noise
psdf = psd;
psdf(L-fH:L) = 0;
% find peak & estimate peak inensity
fpeak = find(psdf==max(psdf));
I_peak = sum(psd(fpeak-df:fpeak+df));

confidence = I_peak/I_hfn;

range_out = (1:L);
psd_out = psd(range_out);


% THIS CODE IS FOR DEBUGGING - LIKELY, IT WOULD BE BETTER TO SWITCH TO "dspdata.psd"
% figure(22);
% pxx = log(periodogram(profile));
% hpsd = dspdata.psd(profile);
% subplot(2,1,1);
% plot(pxx);
% subplot(2,1,2);
% plot(hpsd);




