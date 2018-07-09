function [directionality,velocity,velocity_sd] = quantify_track(x,y,z,mode)

D = [x;y;z];
N = numel(x);

directionality = 0;
velocity = 0;
velocity_sd = 0;

Ndirectionality = 0;
for k=1:N-2
    
    if strcmp(mode,'3d')
        v1 = D(:,k+2)-D(:,k+1);
        v2 = D(:,k+1)-D(:,k);
    elseif strcmp (mode,'noZ')
        v1 = D(1:2,k+2)-D(1:2,k+1);
        v2 = D(1:2,k+1)-D(1:2,k);        
    end
    
    normv1 = norm(v1);
    normv2 = norm(v2);
    if 0~=normv1 && 0~=normv2
        Ndirectionality = Ndirectionality + 1;
        directionality = directionality + dot(v1,v2)/normv1/normv2;
    end
    velocity = velocity + normv1;
end;

directionality = directionality/Ndirectionality;
velocity = velocity/(N-1);

