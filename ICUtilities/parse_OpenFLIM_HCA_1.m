
function res = parse_OpenFLIM_HCA_1(s)

% Well=C10 X=97219.1621093746 Y=27245.513671875 T=0.0 Filterset=green Z=0.0 ID=00569 Laser intensity=21.21.ome.tiff

res = [];

try
    t = {' Well=',' X=',' Y=',' T=',' Filterset=',' Z=',' ID=',' Laser intensity=','.ome.tiff'};
    begs = cell2mat(regexp(s,t));
    for k=1:numel(t)-1    
        L = length(t{k});
        cur_val = s(begs(k)+L:begs(k+1)-1);
        res = [res; {cur_val}];
    end    
    return;    
catch
end

try
    t = {'Well=','_X=','_Y=','_T=','_Z=','.OME.tiff'};
    begs = cell2mat(regexp(s,t));
    for k=1:numel(t)-1    
        L = length(t{k});
        cur_val = s(begs(k)+L:begs(k+1)-1);
        res = [res; {cur_val}];
    end    
    return;    
catch
end

end