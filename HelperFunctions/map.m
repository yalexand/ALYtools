function z=map(x,ymin,ymax)

xmin = min(x(:));
xmax = max(x(:));

z = ymin+(x-xmin)*(ymax-ymin)/(xmax-xmin);







