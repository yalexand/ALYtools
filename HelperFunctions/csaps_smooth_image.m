%--------------------------------------------------------------------------
function smoothedImage = csaps_smooth_image(u,p)

    imageData = double(u);  % Convert to double precision if not already
    
    % Get the size of the image
    [rows, cols] = size(imageData);
    
    % Create meshgrid for x and y coordinates
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % Flatten the grids and the image data
    x = 1:cols;
    y = 1:rows;
    z = imageData';
       
    % Apply cubic smoothing spline
    sp = csaps({x, y}, z, p);
    
    % Evaluate the smoothing spline over a grid
    smoothedData = fnval(sp, {x, y});
    
    % Transpose back the smoothed data
    smoothedImage = smoothedData';

end

