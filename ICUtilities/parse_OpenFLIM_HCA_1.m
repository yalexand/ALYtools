
function res = parse_OpenFLIM_HCA_1(s,mode)

% Well=C10 X=97219.1621093746 Y=27245.513671875 T=0.0 Filterset=green Z=0.0 ID=00569 Laser intensity=21.21.ome.tiff

res = [];

if 1==mode
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
end

if 2==mode
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

%Well=D8 X=77622.03906249983 Y=38745.16796875001 T=0.0 Filterset=Unknown Z=40.0 ID=00016 Laser intensity=3.03.ome.tiff
if 3==mode
    try
        t = {'Well=',' X=',' Y=',' T=',' Filterset=',' Z=',' ID=',' Laser intensity=','.ome.tiff'};
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






end