function [pNorm,V,xMean] = bestfitplane(x)
% BESTFITPLANE one line description.
% Usage:
%   [pNorm,V,xMean] = bestfitplane(x)
% Where:
%   x is n*3 matrix of points
%   pNorm is the normal to the plane
%   V has the orthonormal basis in its columns (pNorm = V(:,1);)
%   xMean is on the plane and is mean(x,1)
%
% BESTFITPLANE: Computes the plane that fits best (lest square of the
% normal distance to the plane) a set of sample points.

% Author: Nick Linton (2013) based on code by Adrien Leygue
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------

% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------
    persistent isWarned
    if isempty(isWarned)
        warning('this should be rewritten using single value decomposition');
        isWarned = true;
    end
    %the mean of the samples belongs to the plane
    xMean = mean(x,1);
    
    %The samples are reduced:
    R = bsxfun(@minus,x,xMean);
    %Computation of the principal directions if the samples cloud
    [V,~] = eig(R'*R);
    %Extract the output from the eigenvectors
    pNorm = V(:,1);
end      
            
            
            