function [z,smoothed] = TD_high_pass_filter( s, avr_size )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mask = ones(1,avr_size)/avr_size;
z = padarray(s,avr_size,'symmetric');
z = conv(z,mask,'same');
smoothed = z(avr_size+1:avr_size+length(s));
z = s - smoothed ;
end

