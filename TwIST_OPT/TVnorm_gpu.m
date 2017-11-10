function y = TVnorm_gpu(x)
[xi,yi]=gradient(x);
y = sum(sum(sqrt(xi.^2+yi.^2)));
