function Rg = GyrationRadius(X,Y,Z,X0,Y0,Z0)

    N = numel(X);
    r = [X Y Z] - repmat([ X0 Y0 Z0],[N 1]);
    Rg = sqrt(sum(sum(r.*r))/N);
    %
    % for sphere, Rg = R/sqrt(5/3);
    %
end

