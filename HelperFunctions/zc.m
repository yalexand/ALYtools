% ZC number of zero crossings in x
% [n] = zc(x) calculates the number of zero crossings in x

function [n] = zc(x)

s = sign(x);

t = filter([1 1],1,s);

n = length(find(t==0)) + length(find(x==0));