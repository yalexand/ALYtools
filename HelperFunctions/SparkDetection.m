function Parameters = SparkDetection(Parameters),
% SPARKDETECTION    Spark detection based upon a generalized second-derivative test.
%  
%   INPUT: A structure, PARAMETERS, with the following fields:
%  
%   PreProcessedData: Represents {\Delta}F/F_0. All pre-processing of the data, 
%       including background subtraction and smoothing (of whatever type) should
%       be accomplished prior to input into the program.
%       Dimensions: The data matrix can be 2- or 3-D, i.e., (x,t) linescan or
%           (x,y,t) full-frame data.
%  
%   RegionOfInterestMask: a logical matrix, representing the region of interest
%       (ROI) in the field of view, containing ones for pixels that are within 
%       the ROI, zeros for those that are not. 
%       Dimensions: If the data matrix is m x n x p, with the spatial 
%           information contained in dimensions m and n, RegionOfInterestMask
%           is of size m x n. 
%  
%   ThresholdHigh: The higher threshold for spark detection. Set by the user to 
%       an appropriate value. Cheng et al. (1999) use a value of 3.5*sigma,
%       where sigma is the standard deviation of the noise.
%  
%   OUPUT: A structure, PARAMETERS, with the above fields intact, plus two 
%       additional fields:
%  
%   CurrentSparkMap: A matrix containing a set of connected objects ("regions")
%       with each region representing a spark receving a label. 
%       The elements of CurrentSparkMap are non-negative integer values. The 
%       pixels labeled 0 are the background. The pixels labeled 1 make up the 
%       first spark, the pixels labeled 2 make up the second spark, etc.
%       Dimensions: Same size as PreProcessedData.
%  
%   DetectedPeaks: A logical matrix with the locations of the spark maxima
%       receiving a 1, and 0's elsewhere.
%       Dimensions: Same size as PreProcessedData.
%  
%   SUB-FUNCTIONS:
%   Determinant: Computes the determinant of a cell matrix, treating each
%       individual n-D cell as a matrix element.
%     
%   BestStorageClass: Attempts to compress storage of an input matrix of
%       integers by casting it to the class with lowest neccesary precision to
%       hold the maximum value.
%  
%   ComputeBoundingBox: Computes the smallest rectangle containing a region, for
%       all the regions in the input matrix. Output is of the form [corner 
%       width], where 'corner' is in the form [x y z ...] and specifies the 
%       upper left corner of the bounding box, and 'width' is in the form 
%       [x_width y_width ...] and specifies the width of the bounding box along 
%       each dimension. Adapted from MATLAB's 'regionprops' function
%  
%   TOOLBOXES REQUIRED (all R14 or higher): 
%     ImageProcessing: bwlabeln, watershed, bwdist, bwulterode (all essential)
%  
%   Spline: csaps (optional: for derivative calculation, can be substituted with
%       'gradient' function, but runs the risk of lower accuracy)
%  
%   Statistics: regress (optional: can probably be substituted with appropriate
%       matrix division	for same result)
%  
%   Copyright (C) 2006 by Mark-Anthony Bray, Disease Biophysics Group (DBG), School of
% 	Engineering and Applied Sciences. 
%   Email contact: b r a y m p {at} s e a s {dot} h a r v a r d {dot} e d u.
%   DBG director: Prof. Kevin Kit Parker
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.


% Extract the pre-proceeded data (i.e., {\Delta}F/F_0) for ease
pre_processed_data = Parameters.PreProcessedData;
num_dims = ndims(pre_processed_data);

% If there is a region-of-interest (ROI) defined, create the mask for the data
mask = true(size(pre_processed_data));
subscripts = cell(1,num_dims); 
for i = 1:num_dims-1, 
    subscripts{i} = 1:size(pre_processed_data,i); 
end
% Assignment to the mask matrix using a cell subscript takes too long for some reason, so I use a loop
for i = 1:size(mask,num_dims), 
    subscripts{num_dims} = i;
    mask(subscripts{:}) = Parameters.RegionOfInterestMask;
end

% Threshold the data, apply the mask, and label the regions 
values_hi = pre_processed_data > Parameters.ThresholdHigh; 
values_hi(~mask) = 0;
[values_hi,numberofregions_hi] = bwlabeln(values_hi);
values_hi = cast(values_hi,BestStorageClass(values_hi));

stats_hi = ComputeBoundingBox(values_hi);
for i = 1:length(stats_hi), 
    stats_hi(i).BoundingBox = ceil(stats_hi(i).BoundingBox); 
    % Place the BoundingBox stats in a row-column rather than x-y order
    stats_hi(i).BoundingBox([1 2 num_dims+[1 2]]) = stats_hi(i).BoundingBox([2 1 num_dims+[2 1]]);
end

% Detected_peaks: A logical matrix for which the elements at which a spark peak is located are set to TRUE
%   Can be useful as an easy way to keep track of where the maximum is for each spark
%
%   Since detected_peaks is binary, using a sparse matrix is surprisingly more efficient than storing the indices
%   Sparse(m,n) abbreviates sparse([],[],[],m,n,0) for which the first index (m) is the linear index in the first 2
%       dimensions, and the second index (n) contains the third dimension (if any), and the whole matrix is set to 0
detected_peaks = logical(sparse(size(values_hi,1)*size(values_hi,2), ...
                                prod(size(values_hi))/(size(values_hi,1)*size(values_hi,2))));

hw = waitbar(0,'Running SparkDetection, please wait');

% Loop over the number of regions
for i = 1:numberofregions_hi,

    if ~isempty(hw), waitbar(i/numberofregions_hi,hw); drawnow, end;  
    
    % Initialize indices to extract the sub-matrix of data from the full data set
    indices = cell(1,num_dims);
    for j = 1:num_dims,
        indices{j} = stats_hi(i).BoundingBox(j)+(0:stats_hi(i).BoundingBox(j+num_dims)-1);
    end

    subdata = pre_processed_data(indices{:});

    if all(stats_hi(i).BoundingBox((num_dims+1):end) >= 3),     % Minimum size requirement for csaps
        
        % Calculate first derivatives (adapted from MATLAB gradient algorithm)
        first_derivatives = cell(1,num_dims);
        xyz = cell(1,num_dims);
        for k = 1:num_dims, 
            xyz{k} = 1:size(subdata,k); 
        end
        subscripts = cell(1,num_dims);
        order = [2:num_dims 1];  % Cyclic permutation

        for j = 1:num_dims,
            first_derivatives{j} = zeros(size(subdata),'single');
            % Unfortunately, my attempts to vectorize the spline calculation have ended up with out-of-memory 
            %   errors for large subdata sets. Therefore, I have to step through each dimension individually, at 
            %   the cost of some generality
            subscripts = cell(1,3);
            for k = 1:3, subscripts{k} = 1:size(subdata,k); end
            for k = 1:size(subdata,3),
                subscripts{end} = k;
                
                % Create spline representation for better spatial derivatives (MATLAB's gradient function 
                %   takes centered differences, which is bad for small maxima)
                first_derivatives{j}(subscripts{:}) = cast(fnval(fndir(csaps(xyz(1:2),double(subdata(subscripts{:})),1),[1;0]), xyz(1:2)),class(first_derivatives{j}));
            end
            first_derivatives{j} = ipermute(first_derivatives{j},[j:num_dims 1:j-1]);

            % Set up for next pass through the loop
            subdata = permute(subdata,order); xyz = xyz(order);
        end

        % Find the actual maxima by computing zero crossings of derivatives to obtain extrema        
        % Since computation of level-sets is different for 2-D and 3-D in MATLAB, I'm forced to be specific here
        actual_maxima = cell(1,num_dims);
        switch num_dims,
        case 2,
            for j = 1:num_dims,
                actual_maxima{j} = contourc(double(first_derivatives{j}),[0 0]); 
                actual_maxima{j}(:,find(actual_maxima{j}(1,:) == 0)) = [];
                % Contourc outputs the coordinates row-wise, so make it column-wise here
                actual_maxima{j} = actual_maxima{j}';
            end
        case 3,
            for j = 1:num_dims,
                [junk, actual_maxima{j}] = isosurface(first_derivatives{j},0);
            end
        end
        
        % levelset_intersection: the intersection of the level sets in all dimensions
        %   Since I perform this calculation by forcing to a regular grid, I may lose some resolution
        levelset_intersection = intersect(round(actual_maxima{1}),round(actual_maxima{2}),'rows');
        for j = 2:num_dims-1,
            levelset_intersection = intersect(levelset_intersection,round(actual_maxima{j+1}),'rows');
        end

        % Convert levelset_intersection matrix into an index vector for the next step
        if ~isempty(levelset_intersection),
            idx = cell(1,num_dims);
            for j = 1:num_dims,
                idx{j} = levelset_intersection(:,j);
            end
            % Switch to row-column form (contourc and isosurface does x-y)
            idx([1 2]) = idx([2 1]);
            idx = sub2ind(size(subdata),idx{:});
        else
            idx = [];
        end
        
        % Calculate second derivatives
        % The for-loop plus permutation is adapted from MATLAB's 'gradient' function, but
        %   the code uses the function 'csaps' from the Spline toolbox, for better precision
        %   However, if the scalar field isn't too rough, the 'gradient' function can be used.
        second_derivatives = cell(num_dims,num_dims);
        subscripts = cell(1,3);
        order = [2:num_dims 1];  % Cyclic permutation

        for j = 1:num_dims,
            xyz = cell(1,num_dims);
            for k = 1:num_dims, xyz{k} = 1:size(first_derivatives{j},k); end

            for k = 1:num_dims,
                if k < j,
                    second_derivatives{j,k} = second_derivatives{k,j};  % Cross products are identical
                else
                    second_derivatives{j,k} = zeros(size(first_derivatives{j}),'single');

                    for m = 1:3, subscripts{m} = 1:size(first_derivatives{j},m); end
                    for m = 1:size(first_derivatives{j},3),
                        subscripts{end} = m;
                        second_derivatives{j,k}(subscripts{:}) = cast(fnval(fndir(csaps(xyz(1:2),double(first_derivatives{j}(subscripts{:})),1),[1;0]), xyz(1:2)),class(first_derivatives{j})); 
                    end
                    second_derivatives{j,k} = ipermute(second_derivatives{j,k},[k:num_dims 1:k-1]);
                end
                first_derivatives{j} = permute(first_derivatives{j},order); xyz = xyz(order);
            end
        end
        % Clear some of the larger variables to free up memory
        clear first_derivatives;

        % Define the kth leading principal minor M_k as the determinant of the k x k submatrix 
        %   of the n x n matrix H obtained by deleting the last (n-k) rows and columns from H
        %
        % Define the Hessian matrix H as the Jacobian of dg/dx_1, dg/dx_2, ..., dg/dx_n of a 
        %   multivariate function g(x_1,x_2,...x_n)
        %
        % The function g has a local maxima at (x1_0,x2_0,x2_0) if its Hessian matrix H 
        %   at (x1_0,x2_0,x2_0) is negative definite, i.e., (-1)^k*M_k > 0, for each k = 1,...,n

        relative_maxima = true;
        for k = 1:num_dims,
            % Create the principal minor submatrix from the second derivative matrix
            principal_minor = second_derivatives; 
            order = num_dims-k+1;
            principal_minor = principal_minor(1:order,1:order);
            
            % Rather than compute the determinant for every point (x1_0,x2_0,x2_0), I've written 
            %   my own determinant function to vectorize the operation over the n-D cell submatrix (see below)
            relative_maxima = relative_maxima & ( ((-1)^order)*Determinant(principal_minor) > 0 );
        end
        
         % 'Clean' to remove spurious maxima (those with area = 1)
        relative_maxima = bwareaopen(relative_maxima,2);

        % Clear some of the larger variables to free up memory
        clear principal_minor second_derivatives;

        % Mask relative maxima according to upper threshold limit and label them
        relative_maxima = relative_maxima & (values_hi(indices{:}) == i);
        [relative_maxima,number_of_relative_maxima] = bwlabeln(relative_maxima); 
        relative_maxima = cast(relative_maxima,BestStorageClass(relative_maxima));

        if ~isempty(idx),   % This condition should almost always be satisfied, but occasionally 
            %                   not for really small subdata
            % Initialize the extrema matrix
            extrema = false(size(subdata)); 
            extrema(idx) = true; 
            
            % Mask extrema using the relative maxima, erode them to a point, label, and then find the indices
            actual_maxima = bwlabeln(bwulterode(extrema & (relative_maxima > 0)));
            actual_maxima_idx = find(actual_maxima);

            % In the following portion, relative maxima which do not have a corresponding actual maxima are eliminated.
            %   However, until a better method of detecting the actual maxima is found, there may be a chance that
            %   a relevant relative maxima is inadvertently removed. It may be a good idea to tag these maxima for
            %   later inspection rather than eliminate them
            [relative_maxima,number_of_relative_maxima] = bwlabeln(ismember(relative_maxima,relative_maxima(actual_maxima_idx)));
            relative_maxima = cast(relative_maxima,BestStorageClass(relative_maxima));
            
            % Prepare subscript matrix for detected_peaks by collpsing the 1st and 2nd dimensions into the first, 
            %   moving the third dimension into the second, and eliminating the rest
            subscripts = cell(1,num_dims);
            [subscripts{:}] = ind2sub(size(subdata),idx);
            for j = 1:num_dims,
                subscripts{j} = subscripts{j} + stats_hi(i).BoundingBox(j)-1;
            end
            subscripts{1} = sub2ind(size(values_hi),subscripts{1:max(2,end-1)}); 
            if size(subdata,3) == 1,
                subscripts{2} = 1; 
            else
                subscripts(2) = subscripts(end);
            end
            subscripts(3:end) = [];

            % The new detected_peaks is found simply by ORing it with the previous result
            % Use notation sparse(i,j,s,m,n) where i,j are the subscripts that are to be set to the value s (i.e, 1)
            %   and m and n are the number of rows (which are the collapsed 1st and 2nd dimensions) and columns 
            %   (3rd dimension, if any), respectively
            detected_peaks = detected_peaks | sparse(double(subscripts{1}), double(subscripts{2}), 1, ...
                                                     size(values_hi,1)*size(values_hi,2), ...
                                                     prod(size(values_hi))/(size(values_hi,1)*size(values_hi,2)));

            if number_of_relative_maxima > 1,    
                % Check if relative maxima are joined. If so, separate them using the actual maxima 
                %   within the relative maxima blobs
                changed_number_of_obj = 0;
                
                for j = 1:number_of_relative_maxima,
                    idx = find(relative_maxima == j);
                    idx = find(ismember(actual_maxima_idx,idx));
                    if length(idx) > 1,
                        submaxima = false(size(subdata)); 
                        submaxima(actual_maxima_idx(idx)) = true;
                        
                        % Compute the watershed lines and use them to cut apart the relative maxima
                        dividing_surface = (watershed(bwdist(submaxima)) == 0) & (relative_maxima == j);
                        relative_maxima(dividing_surface) = 0;
                        changed_number_of_obj = changed_number_of_obj + 1;
                    end
                end
                
                % Re-label the relative maxima
                relative_maxima = relative_maxima > 0;
                [relative_maxima,number_of_relative_maxima] = bwlabeln(relative_maxima); 
                relative_maxima = cast(relative_maxima,BestStorageClass(relative_maxima));
            end
        end

        % Attempt to fit n-plane to each of the relative maxima; if above 95% confidence fit, remove  
        %   that maxima from consideration since is close to being planar and probably a 'false' positive 
        changed_number_of_obj = 0;
        for j = 1:number_of_relative_maxima,
            if length(find(relative_maxima == j)) == 1, 
                % Degenerate case; remove relative maxima immediately
                relative_maxima(relative_maxima == j) = 0; 
                changed_number_of_obj = changed_number_of_obj + 1;
            else
                idx = cell(1,num_dims);
                [idx{:}] = ind2sub(size(relative_maxima),find(relative_maxima == j)); 
                v = subdata(relative_maxima == j); v = v(:); 
                
                % Turn off rank-deficient warning in 'regress' function
                warning('off','stats:regress:RankDefDesignMat'); lastwarn('');
                [junk,junk,junk,junk,regress_stats] = regress(double(v),[[idx{:}] ones(size(v))]); 
                [msgstr,msgid] = lastwarn;
                
                % If regression fails (for whatever reason), automatically eliminate the relative maxima...
                if strcmp(msgid,'stats:regress:RankDefDesignMat'),
                    relative_maxima(relative_maxima == j) = 0; 
                    changed_number_of_obj = changed_number_of_obj + 1; 
                else
                % ... Otherwise, fit the plane and check for goodness of fit
                    R2 = regress_stats(1);
                    if R2 > 0.95, relative_maxima(relative_maxima == j) = 0; 
                        changed_number_of_obj = changed_number_of_obj + 1; end
                end
                warning('on','stats:regress:RankDefDesignMat');
            end
        end

       % If the number of regions has changed, re-label the matrix
       if changed_number_of_obj, 
            [relative_maxima,number_of_relative_maxima] = bwlabeln(relative_maxima); 
            relative_maxima = cast(relative_maxima,BestStorageClass(relative_maxima));
        end 

        % Update the high-threshold values
        if number_of_relative_maxima > max(values_hi(:)),
            values_hi = cast(values_hi,BestStorageClass(relative_maxima));
        end
        subdata = values_hi(indices{:});
        idx = find(subdata == i);
        warning('off','MATLAB:intConvertNonIntVal');
        subdata(idx) = floor(cast(relative_maxima(idx),class(values_hi)));           
        warning('on','MATLAB:intConvertNonIntVal');
        values_hi(indices{:}) = subdata;
    else
        % Update the high-threshold values
        values_hi(indices{:}) = 0;
    end
end

if ~isempty(hw), delete(hw), drawnow; end;

% Re-label the highj-threshold values and cast to optimal datatype
values_hi = bwlabeln(values_hi);
values_hi = cast(values_hi,BestStorageClass(values_hi)); 

% Set up output
Parameters.CurrentSparkMap = values_hi;

% 'Un-fold' the detected_peaks sparse matrix back into a full matrix
Parameters.DetectedPeaks = false(size(values_hi));
subscripts = cell(1,num_dims); 
[subscripts{1},subscripts{num_dims}] = find(detected_peaks); 
[subscripts{1:2}] = ind2sub(size(Parameters.DetectedPeaks),subscripts{1});
index = sub2ind(size(Parameters.DetectedPeaks),subscripts{:});
Parameters.DetectedPeaks(sub2ind(size(Parameters.DetectedPeaks),subscripts{:})) = true;
    
%====================================
%%%
%%%  Action - Determinant
%%%
%
% DETERMINANT: Recursive definition of determinant using expansion by minors.
%
%   Based on code obtained from http://astronomy.swin.edu.au/~pbourke/analysis/determinant/

function value = Determinant(a),

n = size(a,1);

if n < 1,       % Error
elseif n == 1,  % Shouldn't get used
    value = a{1,1};
elseif n == 2, 
    value = a{1,1}.*a{2,2} - a{2,1}.*a{1,2};
else
    value = 0; 
    i = 1;
    for j = 1:n,
        m = a; m(i,:) = []; m(:,j) = [];
        value = value + ((-1)^(i+j) * a{i,j} .* Determinant(m));
    end
end

%====================================
%%%
%%%  Action - BestStorageClass
%%%
%
% BESTSTORAGECLASS: Returns the string of the datatype which is optimal for storing a matrix
%   in memory.
%
% For storing some of the label matrices, it's best to cast the matrix it to a class which is 
%   appropriate for the number of objects represented in order to save memory
function classname = BestStorageClass(input_matrix)

% The notation below is appropriate for full and sparse matrices
max_value = full(max(input_matrix(find(input_matrix))));
min_value = full(min(input_matrix(find(input_matrix))));

if min_value >= 0 & max_value <= 1,
    classname = 'logical';
    return;
end

possible_classes = {'uint8','int8','uint16','int16','uint32','int32','uint64','int64'};
for i = possible_classes,
    if min_value >= intmin(char(i)) & max_value <= intmax(char(i)),
        classname = char(i);
        return;
    end
end

classname = 'single';

%====================================
%%%
%%% ComputeBoundingBox
%%%
%
% COMPUTEBOUNDINGBOX: A stripped-down version of REGIONPROPS(A,'BoundingBox')
%
% This function is a compilation of the code from the ComputeBoundingBox, 
%   ComputePixelList and ComputePixelIdxList sub-functions in the function 
%   REGIONPROPS.
%
% It was implemented because the sub2ind call in ComputePixelList uses
%   double precision, which can give out-of-memory errors if the input is large.

function stats = ComputeBoundingBox(L),

num_dims = ndims(L);

if isempty(L),
	numObjs = 0;
else
    numObjs = round(double(max(L(:))));
end

% Initialize the stats structure array.
stats = cell2struct(cell(1, numObjs), 'BoundingBox', 1);

% Create a sparse matrix containing one column per region.
S = sparse(find(L), double(L(find(L))), 1);

% Loop over each column of the sparse matrix.   
for k = 1:numObjs,
    % Finding the row indices of the nonzero entries in S(:,P) is equivalent to finding the 
    %   linear indices of pixels in L that equal P.
    % PixelIdxList: A P-by-1 matrix, where P is the number of pixels belonging to the region.  
    %   Each element contains the linear index of the corresponding pixel.
    PixelIdxList = single(find(S(:,k)));
    
    if ~isempty(PixelIdxList),
        % Code ripped from ind2sub to replace [In{:}] = ind2sub(size(L), stats(k).PixelIdxList) 
        %   with single precision calculation to save memory
        
        % Convert the linear indices to subscripts and store the results in the pixel list.
        % PixelList: A P-by-2 matrix, where P is the number of pixels belonging to the region. 
        %   Each row contains the row and column coordinates of a pixel.
        PixelList = zeros(size(PixelIdxList,1),num_dims,'single');
        siz = size(L);
        nout = max(ndims(L),1);
        if length(siz) <= nout,
            siz = [siz ones(1,nout-length(siz))];
        else
            siz = [siz(1:nout-1) prod(siz(nout:end))];
        end
        n = length(siz);
        m = [1 cumprod(siz(1:end-1))];
        PixelIdxList = PixelIdxList - 1;
        for i = n:-1:1,
            PixelList(:,i) = single(floor(PixelIdxList/m(i))+1); 
            PixelIdxList = rem(PixelIdxList,m(i));
        end

        % Reverse the order of the first two subscripts to form x-y order
        PixelList = PixelList(:,[2 1 3:end]);
    else
        PixelList = zeros(0,num_dims,'single');
    end
    
    if isempty(PixelList),
        stats(k).BoundingBox = [0.5*ones(1,num_dims) zeros(1,num_dims)];
    else
        min_corner = min(PixelList,[],1) - 0.5;
        max_corner = max(PixelList,[],1) + 0.5;
        stats(k).BoundingBox = [min_corner (max_corner - min_corner)];
    end
end
