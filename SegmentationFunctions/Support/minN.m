function [v,idx] = minN(A)
%MINN - gets minimum value and its index of n-dimensional matrix
%
% Syntax:  [v,idx] = minN(A)
%
% Inputs:
%    A   - the matrix of which the minimum value shall be found
%
% Outputs:
%    v   - the minimum value of A
%    idx - the index of v in matrix A
%
% Example: 
%    A(:,:,1) = [1,1,1;1,1,1;1,1,1];
%    A(:,:,2) = [1,1,1;1,1,1;1,1,1];
%    A(:,:,3) = [1,1,1;1,1,1;0,1,1];
%    
%    [v,idx] = minN(A);
%
%    % v = 0, and idx = [3,1,3]
%    % it is A(3,1,3) = 0
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: min, max, maxN

% Author: Christopher Haccius, significantly enhanced by Jos (Matlab
% community member)
% Telecommunications Lab, Saarland University, Germany
% email: haccius@nt.uni-saarland.de
% December 2013; Last revision: 19-December-2013

[v, linIdx] = min(A(:));
[idxC{1:ndims(A)}] = ind2sub(size(A),linIdx);
idx = cell2mat(idxC);