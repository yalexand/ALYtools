
function merging_3channel_OpenFLIM_HCA(src,dst)

% 'Z:\User data\Yuiry Alexandrov\Macrophages_test\Macrophages_2_40x_uv_lamp\Sequenced FLIM acquisition 2018-02-21_19.45.38';
% 'Z:\User data\Yuiry Alexandrov\Macrophages_test\Macrophages_2_40x_uv_lamp\Sequenced FLIM acquisition 2018-02-21_19.45.38_merged';

ext = '*.ome.tiff';

D = dir( fullfile(src,ext) );
if isempty(D), return, end;

fnames = {D.name};

channelnames = [];
wxytfz = [];% wxytfz pairs

separator = 'alksdlajdlalnxmkortjfnqoke';
for k=1:numel(fnames)
     res = parse_OpenFLIM_HCA_1(char(fnames(k)));
     w = char(res(1)); 
     x = char(res(2)); 
     y = char(res(3));
     t = char(res(4));
     z = char(res(6));
     %
     wxytfz = [wxytfz; {[w separator x separator y separator t separator z]}];
     %
     channelnames = [channelnames; res(5)];
end

wxytfz=unique(wxytfz);
channelnames = unique(channelnames);

L1 = [];
L2 = [];
L3 = [];

for k=1:numel(wxytfz)
     res = strsplit(char(wxytfz(k)),separator);
     w = char(res(1)); 
     x = char(res(2)); 
     y = char(res(3));
     t = char(res(4));
     z = char(res(5));
     n1 = [];
     n2 = [];
     n3 = [];
     for m=1:numel(fnames)
         str = char(fnames(m));
         %
         res = parse_OpenFLIM_HCA_1(str);
         w_m = char(res(1));
         x_m = char(res(2));
         y_m = char(res(3));
         t_m = char(res(4));
         f_m = char(res(5));
         z_m = char(res(6));
         %    
         if     strcmp(w,w_m) && strcmp(x,x_m) && strcmp(y,y_m) && strcmp(t,t_m) && strcmp(z,z_m) && strcmp(f_m,channelnames(1))
             n1 = str;
         elseif strcmp(w,w_m) && strcmp(x,x_m) && strcmp(y,y_m) && strcmp(t,t_m) && strcmp(z,z_m) && strcmp(f_m,channelnames(2))
             n2 = str;
         elseif strcmp(w,w_m) && strcmp(x,x_m) && strcmp(y,y_m) && strcmp(t,t_m) && strcmp(z,z_m) && strcmp(f_m,channelnames(3))
             n3 = str;
         end
     end
     L1 = [L1; {n1}];
          L2 = [L2; {n2}];
               L3 = [L3; {n3}];
end

I = [];

for k=1:numel(L1)
    %
    f1 = [src filesep char(L1(k))];
    f2 = [src filesep char(L2(k))];
    f3 = [src filesep char(L3(k))];
    %
    res = parse_OpenFLIM_HCA_1(f1);
    basename = ['Well=' char(res(1)) '_X=' char(res(2)) '_Y=' char(res(3)) '_T=' char(res(4)) '_Z=' char(res(6))];
    %
    [uppixX,uppixZ,U1] = bfopen_v(f1);
    [uppixX,uppixZ,U2] = bfopen_v(f2);
    [uppixX,uppixZ,U3] = bfopen_v(f3);
    if isempty(I)
        [szX,szY]=size(U1);
        I = zeros(szX,szY,1,3,1,'like',U1);
    end
        I(:,:,1,1,1)=U1;
        I(:,:,1,2,1)=U2;
        I(:,:,1,3,1)=U3;
        %
        metadata = createMinimalOMEXMLMetadata(I,'XYTCZ');
            if ~isempty(uppixX) && ~isempty(uppixZ)
                toPosFloat = @(x) ome.xml.model.primitives.PositiveFloat(java.lang.Double(x));
                metadata.setPixelsPhysicalSizeX(toPosFloat(uppixX),0);
                metadata.setPixelsPhysicalSizeY(toPosFloat(uppixX),0);
                metadata.setPixelsPhysicalSizeZ(toPosFloat(uppixZ),0);
            end
        bfsave(I,[dst filesep basename '.OME.tiff'],'metadata',metadata);
end

