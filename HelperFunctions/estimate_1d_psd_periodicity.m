function [ fpeak0,R ] = estimate_1d_psd_periodicity( in_psd,df,hf_cutoff)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% find peak & estimate peak inensity

psd = in_psd;
L = length(psd);
psd(L-hf_cutoff:L)=0;

fpeak0 = max(find(psd==max(psd)));

%left = max(1,fpeak0-df);
%right = min(L,fpeak0+df);
%I_periodic = sum(psd(left:right));

I_periodic = 0;
delta_f = length(psd) - fpeak0;
for m=1:4,
    fpeak = length(psd) - m*delta_f;
    if fpeak < 2*df, break, end;
    I_periodic = I_periodic + sum(psd(max(fpeak-df,1):min(fpeak+df,length(psd))));
end

R = I_periodic/sum(psd);

if R > 1, R=0; end; % wow

end

