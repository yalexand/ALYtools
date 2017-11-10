function MFout=conv3Dfreq(cprocl,vker)

cprocl=double(cprocl);
smoothcell=zeros(size(cprocl));
centerpix=floor(size(cprocl)/2);
centerker=floor(size(vker)/2);

embed1i=centerpix(1)-centerker(1)+[1:size(vker,1)];
embed2i=centerpix(2)-centerker(2)+[1:size(vker,2)];
embed3i=centerpix(3)-centerker(3)+[1:size(vker,3)];

smoothcell(embed1i,embed2i,embed3i)=vker;

% disp('matched filter correlating observed volume');
% Correlation through FFT needs:
% a. Finding FFT(S) and FFT(C) 
FFT1=fftn(cprocl);
FFT2=fftn(smoothcell);
% b. ComplexArray=FFT(S) * Conj{FFT(C)} 
CPXARR2=FFT1.*conj(FFT2);
% c. IFFT{ComplexArray}
MFout=2*abs(fftshift(ifftn(CPXARR2)));