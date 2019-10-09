function [v,idx] = maxN(A)
%MAXN - gets maximum value and its index of n-dimensional matrix
%
% Syntax:  [v,idx] = maxN(A)
%
% Inputs:
%    A   - the matrix of which the maximum value shall be found
%
% Outputs:
%    v   - the maximum value of A
%    idx - the index of v in matrix A
%
% Example: 
%    A(:,:,1) = [0,0,0;0,0,0;0,0,0];
%    A(:,:,2) = [0,0,0;0,0,0;0,0,0];
%    A(:,:,3) = [0,0,0;0,0,0;1,0,0];
%    
%    [v,idx] = maxN(A);
%
%    % v = 1, and idx = [3,1,3]
%    % it is A(3,1,3) = 1
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: max, min, minN

% Author: Christopher Haccius, significantly enhanced by Jos (Matlab
% community member)
% Telecommunications Lab, Saarland University, Germany
% email: haccius@nt.uni-saarland.de
% December 2013; Last revision: 19-December-2013

[v, linIdx] = max(A(:));
[idxC{1:ndims(A)}] = ind2sub(size(A),linIdx);
idx = cell2mat(idxC);