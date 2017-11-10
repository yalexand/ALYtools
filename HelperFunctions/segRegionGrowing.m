function RG = segRegionGrowing(u,seed,std_factor_int,std_factor_out)

RG = [];

if isempty(seed), return, end;

se1=strel('disk',1,0);

RG = seed;

while(true)
    bnd_out = imdilate(RG,se1) - RG;
    bnd_int = RG - imerode(RG,se1);
    %
    s_out = u(bnd_out~=0);
    s_int = u(bnd_int~=0);
    %
    t_int = mean(s_int(:))-std_factor_int*std(s_int(:));
    t_out = mean(s_out(:))+std_factor_out*std(s_out(:));
    %[t_int t_out]    
    valeursVoisins = ( u > t_int ) & ( u < t_out ) & bnd_out;    
    %
    if 0==sum(valeursVoisins(:)), return, end;    
    %
    RG(valeursVoisins) = 1;
end