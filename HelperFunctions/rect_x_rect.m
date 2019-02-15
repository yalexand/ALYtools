function r3 = rect_x_rect(x1,y1,w1,h1,x2,y2,w2,h2)
    %r1 = [x1 y1 w1 h1];
    %r2 = [x2 y2 w2 h2];
    l = max(x1,x2);
    r = min(x1+w1,x2+w2);
    b = max(y1,y2);
    t = min(y1+h1,y2+h2);
    r3 =  [l b r-l t-b];
end