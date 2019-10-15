function z=map(x,ymin,ymax)

xmin = minN(x);
xmax = maxN(x);

z = ymin+(x-xmin)*(ymax-ymin)/(xmax-xmin);







