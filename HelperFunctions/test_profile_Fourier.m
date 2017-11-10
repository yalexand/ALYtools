clear all;
close all;

data = xlsread('nonperiodic4.xls','Intensity c0 t0 z0');
profile = data(:,2);

fL = 100; % counting from high-freq end
fH = 8;
df = 1;
[range_out,psd_out,fpeak,confidence] = analyze_profile_periodicity_via_Fourier(profile,fL,fH,df); % frequency in 1/pix

figure()
subplot(2,1,1);
plot(1:length(profile),profile,'b.-');
grid on;
subplot(2,1,2);
plot(range_out,psd_out,'k.-');
title(['fpeak = ' num2str(fpeak), ', confidence = ' num2str(confidence)]);
grid on;