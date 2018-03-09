
function merging_3channel(src,dst,t1,t2,t3)
% utility to merge image files' triples into OME.tiff depending on wavelength name templates
% src = 'c:\temp\merging_test_data\src';
% dst = 'c:\temp\merging_test_data\dst';
% t1 = '_Wavelength_380nm';
% t2 = '_Wavelength_450nm';
% t3 = '_Wavelength_560nm';
% merging_3channel('c:\temp\merging_test_data\src','c:\temp\merging_test_data\dst','_Wavelength_380nm','_Wavelength_450nm','_Wavelength_560nm')

ext = '*.tif'; D = dir( fullfile(src,ext) );
if isempty(D), ext = '*.tiff';D = dir( fullfile(src,ext) ); end;
if isempty(D), ext = '*.ome.tiff';D = dir( fullfile(src,ext) ); end;
if isempty(D), return, end;

fnames = {D.name};

L1=[];
L2=[];
L3=[];
for k=1:numel(fnames)
    str = char(fnames(k));
    if strfind(str,t1)
        L1 = [L1; cellstr(str)];
    end
    if strfind(str,t2)
        L2 = [L2; cellstr(str)];
    end
    if strfind(str,t3)
        L3 = [L3; cellstr(str)];
    end
end

L1=sort(L1);
L2=sort(L2);
L3=sort(L3);

I = [];

for k=1:numel(L1)
    % find image file names
    f1 = [src filesep char(L1(k))];
    f2 = [src filesep char(L2(k))];
    f3 = [src filesep char(L3(k))];
    basename = strrep(char(L1(k)),t1,'');
    basename = strrep(basename,ext(2:length(ext)),'');  
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

