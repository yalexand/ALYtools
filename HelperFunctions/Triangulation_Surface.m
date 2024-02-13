function totalArea = Triangulation_Surface(dt)

    c = dt.Points;    
    fb = freeBoundary(dt); 
    
    totalArea = 0;
    for i = 1:size(fb, 1)    
        % Get the vertices of the facet
        v1 = c(fb(i,1),:); 
        v2 = c(fb(i,2),:); 
        v3 = c(fb(i,3),:); 
        % facet area
        area = norm(cross(v2 - v1, v3 - v1)) / 2;
        %
        totalArea = totalArea + area;
    end

