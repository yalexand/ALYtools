
function merging_3channel_OpenFLIM_HCA(src,dst)

% 'Z:\User data\Yuiry Alexandrov\Macrophages_test\Macrophages_2_40x_uv_lamp\Sequenced FLIM acquisition 2018-02-21_19.45.38';
% 'Z:\User data\Yuiry Alexandrov\Macrophages_test\Macrophages_2_40x_uv_lamp\Sequenced FLIM acquisition 2018-02-21_19.45.38_merged';

ext = '*.ome.tiff';

D = dir( fullfile(src,ext) );
if isempty(D), return, end;

fnames = {D.name};

channelnames = [];
x = [];
y = [];
for k=1:numel(fnames)
     str = char(fnames(k));
     res = parse_OpenFLIM_HCA_1(str);
     channelnames = [channelnames; res(5)];
end
channelnames = unique(channelnames);

L1=[];
L2=[];
L3=[];

xys = [];% xy pairs
separator = 'alksdlajdlalnxmkortjfnqoke';
for k=1:numel(fnames)
     str = char(fnames(k));
     res = parse_OpenFLIM_HCA_1(str);
     x = char(res(2)); % supposed to be unique
     y = char(res(3));
     %chan = char(res(5));
     xys = [xys; {[x separator y]}];
end
xys=unique(xys);

for k=1:numel(xys)
     str = char(xys(k));
     r = strsplit(char(xys(k)),separator);
     x = char(r(1));
     y = char(r(2));
     n1 = [];
     n2 = [];
     n3 = [];
     for k=1:numel(fnames)
         str = char(fnames(k));
         if ~isempty(strfind(str,x)) && ~isempty(strfind(str,y)) && ~isempty(strfind(str,char(channelnames(1))))
             n1 = str;
         elseif ~isempty(strfind(str,x)) && ~isempty(strfind(str,y)) && ~isempty(strfind(str,char(channelnames(2))))
             n2 = str;
         elseif ~isempty(strfind(str,x)) && ~isempty(strfind(str,y)) && ~isempty(strfind(str,char(channelnames(3))))
             n3 = str;
         end
     end
     L1 = [L1; {n1}];
          L2 = [L2; {n2}];
               L3 = [L3; {n3}];
end

I = [];

for k=1:numel(L1)
    % find image file names
    f1 = [src filesep char(L1(k))];
    f2 = [src filesep char(L2(k))];
    f3 = [src filesep char(L3(k))];
    %
    res = parse_OpenFLIM_HCA_1(f1);
    basename = ['Well=' char(res(1)) '_X=' char(res(2)) '_Y=' char(res(3))];    
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

