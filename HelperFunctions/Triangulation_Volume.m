function totalVolume = Triangulation_Volume(dt)

c = dt.Points;
tri = dt.ConnectivityList;

totalVolume = 0;
for i = 1:size(tri,1)
    % Get the vertices of the tetrahedron
    vertices = [ c(tri(i, 1),:); c(tri(i, 2),:); c(tri(i, 3),:); c(tri(i, 4),:) ];
    % Calculate the volume of the tetrahedron using the absolute value of the determinant
    volume = abs(det([vertices, ones(4,1)])) / 6;
    %
    totalVolume = totalVolume + volume;    
end