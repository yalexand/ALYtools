%--------------------------------------------------------------------------
function smoothedImage = csaps_smooth_image(u, p)
    % CSAPS_SMOOTH_IMAGE Smooth a 2D image using cubic smoothing splines
    %   smoothedImage = csaps_smooth_image(u, p) smooths the 2D image 'u' using
    %   cubic smoothing splines with a smoothing parameter 'p'.
    %
    %   Inputs:
    %     u - A 2D array representing the image to be smoothed.
    %     p - Smoothing parameter (0 <= p <= 1).
    %
    %   Outputs:
    %     smoothedImage - The smoothed 2D image.

    % Validate inputs
    if ~ismatrix(u)
        error('Input must be a 2D array.');
    end
    
    if ~isscalar(p) || p < 0 || p > 1
        error('Smoothing parameter p must be a scalar between 0 and 1.');
    end
    
    % Convert to double precision if not already
    imageData = double(u);  
    
    % Get the size of the image
    [rows, cols] = size(imageData);
    
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