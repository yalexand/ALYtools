function frame_time = add_plane_to_OMEtiff_with_metadata_ZCT(I, indices, ZCT, folder, ometiffilename, varargin)

global writer;
global getBytes;
global hw;
global ifd;

assert(isnumeric(I), 'First argument must be numeric');
assert(isnumeric(indices), 'Second argument must be numeric');
assert(isnumeric(ZCT), 'Third argument must be numeric');

ip = inputParser;

ip.addOptional('PhysSz_X', [], @isnumeric);
ip.addOptional('PhysSz_Y', [], @isnumeric);
%
ip.addOptional('ModuloZ_Type', [], @ischar);
ip.addOptional('ModuloZ_TypeDescription', [], @ischar);
ip.addOptional('ModuloZ_Unit', [], @ischar);
ip.addOptional('ModuloZ_Start', [], @isnumeric);
ip.addOptional('ModuloZ_Step', [], @isnumeric);
ip.addOptional('ModuloZ_End', [], @isnumeric);
ip.addOptional('ModuloZ_Labels', [], @isnumeric);
%
ip.addOptional('ModuloC_Type', [], @ischar);
ip.addOptional('ModuloC_TypeDescription', [], @ischar);
ip.addOptional('ModuloC_Unit', [], @ischar);
ip.addOptional('ModuloC_Start', [], @isnumeric);
ip.addOptional('ModuloC_Step', [], @isnumeric);
ip.addOptional('ModuloC_End', [], @isnumeric);
ip.addOptional('ModuloC_Labels', [], @isnumeric);
%
ip.addOptional('ModuloT_Type', [], @ischar);
ip.addOptional('ModuloT_TypeDescription', [], @ischar);
ip.addOptional('ModuloT_Unit', [], @ischar);
ip.addOptional('ModuloT_Start', [], @isnumeric);
ip.addOptional('ModuloT_Step', [], @isnumeric);
ip.addOptional('ModuloT_End', [], @isnumeric);
ip.addOptional('ModuloT_Labels', [], @isnumeric);
%
ip.addOptional('Tags', [], @ischar);    
ip.addOptional('Tags_SeparatingSeq', [], @ischar); 
ip.addOptional('verbose', [], @islogical);
ip.addParamValue('Compression', '',  @(x) ismember(x, getCompressionTypes()));
ip.addParamValue('BigTiff', false , @islogical);

ip.parse(varargin{:});

        [sizeX, sizeY] = size(I);
        sizeZ = ZCT(1);
        sizeC = ZCT(2);
        sizeT = ZCT(3);
        final_index = sizeZ*sizeC*sizeT;

z = indices(1);
c = indices(2);
t = indices(3);
index = sub2ind(ZCT,z,c,t);
                
if 1 == index % make all setups

        addpath_OMEkit;
        
        % verify that enough memory is allocated
        bfCheckJavaMemory();
        % Check for required jars in the Java path
        bfCheckJavaPath();
        
        % ini logging
        loci.common.DebugTools.enableLogging('INFO');
        java.lang.System.setProperty('javax.xml.transform.TransformerFactory', 'com.sun.org.apache.xalan.internal.xsltc.trax.TransformerFactoryImpl');

        metadata = createMinimalOMEXMLMetadata(I);
        toInt = @(x) ome.xml.model.primitives.PositiveInteger(java.lang.Integer(x));        
        metadata.setPixelsSizeZ(toInt(sizeZ), 0);
        metadata.setPixelsSizeC(toInt(sizeC), 0);
        metadata.setPixelsSizeT(toInt(sizeT), 0);        

if ~isempty(ip.Results.PhysSz_X)
    metadata.setPixelsPhysicalSizeX(ome.xml.model.primitives.PositiveFloat(java.lang.Double(ip.Results.PhysSz_X)),0);
end    
if ~isempty(ip.Results.PhysSz_Y)
    metadata.setPixelsPhysicalSizeY(ome.xml.model.primitives.PositiveFloat(java.lang.Double(ip.Results.PhysSz_Y)),0);
end

%%%%%%%%%%%%%%%%%%% set up Modulo XML description metadata if present - starts
modlo_cnt = 0;
% ModuloAlongZ
if isfield(ip.Results,'ModuloZ_Type')
    if (~isempty(ip.Results.ModuloZ_Type) && ~isempty(ip.Results.ModuloZ_TypeDescription) && ~isempty(ip.Results.ModuloZ_Unit)) || ~isempty(ip.Results.ModuloZ_Labels)
        modlo = loci.formats.CoreMetadata();        
        %
        modlo.moduloZ.type = ip.Results.ModuloZ_Type;
        modlo.moduloZ.unit = ip.Results.ModuloZ_Unit;
        modlo.moduloZ.typeDescription = ip.Results.ModuloZ_TypeDescription;
    end

    if ~isempty(ip.Results.ModuloZ_Start) && ~isempty(ip.Results.ModuloZ_Step) && ~isempty(ip.Results.ModuloZ_End)
        modlo.moduloZ.start = ip.Results.ModuloZ_Start;
        modlo.moduloZ.end = ip.Results.ModuloZ_End;
        modlo.moduloZ.step = ip.Results.ModuloZ_Step;
        %
    elseif ~isempty(ip.Results.ModuloZ_Labels)
        %
        labels = ip.Results.ModuloZ_Labels;
        modlo.moduloZ.labels = javaArray('java.lang.String',length(labels));
        for i=1:length(labels)
            modlo.moduloZ.labels(i)= java.lang.String(num2str(labels(i)));
        end                                                      
    end
    modlo_cnt = modlo_cnt + 1;
end
% ModuloAlongC
if isfield(ip.Results,'ModuloC_Type')
    if (~isempty(ip.Results.ModuloC_Type) && ~isempty(ip.Results.ModuloC_TypeDescription) && ~isempty(ip.Results.ModuloC_Unit)) || ~isempty(ip.Results.ModuloC_Labels)
        if ~exist('modlo','var')
            modlo = loci.formats.CoreMetadata();        
        end
        %
        modlo.moduloC.type = ip.Results.ModuloC_Type;
        modlo.moduloC.unit = ip.Results.ModuloC_Unit;
        modlo.moduloC.typeDescription = ip.Results.ModuloC_TypeDescription;
    end

    if ~isempty(ip.Results.ModuloC_Start) && ~isempty(ip.Results.ModuloC_Step) && ~isempty(ip.Results.ModuloC_End)
        modlo.moduloC.start = ip.Results.ModuloC_Start;
        modlo.moduloC.end = ip.Results.ModuloC_End;
        modlo.moduloC.step = ip.Results.ModuloC_Step;
        %
    elseif ~isempty(ip.Results.ModuloC_Labels)
        %
        labels = ip.Results.ModuloC_Labels;
        modlo.moduloC.labels = javaArray('java.lang.String',length(labels));
        for i=1:length(labels)
            modlo.moduloC.labels(i)= java.lang.String(num2str(labels(i)));
        end                                                      
    end
    modlo_cnt = modlo_cnt + 1;
end
% ModuloAlongT
if isfield(ip.Results,'ModuloT_Type')
    if (~isempty(ip.Results.ModuloT_Type) && ~isempty(ip.Results.ModuloT_TypeDescription) && ~isempty(ip.Results.ModuloT_Unit)) || ~isempty(ip.Results.ModuloT_Labels)
        if ~exist('modlo','var')
            modlo = loci.formats.CoreMetadata();        
        end
        %
        modlo.moduloT.type = ip.Results.ModuloT_Type;
        modlo.moduloT.unit = ip.Results.ModuloT_Unit;
        modlo.moduloT.typeDescription = ip.Results.ModuloT_TypeDescription;
    end

    if ~isempty(ip.Results.ModuloT_Start) && ~isempty(ip.Results.ModuloT_Step) && ~isempty(ip.Results.ModuloT_End)
        modlo.moduloT.start = ip.Results.ModuloT_Start;
        modlo.moduloT.end = ip.Results.ModuloT_End;
        modlo.moduloT.step = ip.Results.ModuloT_Step;
        %
    elseif ~isempty(ip.Results.ModuloT_Labels)
        %
        labels = ip.Results.ModuloT_Labels;
        modlo.moduloT.labels = javaArray('java.lang.String',length(labels));
        for i=1:length(labels)
            modlo.moduloT.labels(i)= java.lang.String(num2str(labels(i)));
        end                                                      
    end
    modlo_cnt = modlo_cnt + 1;
end
%
if exist('modlo','var') 
    OMEXMLService = loci.formats.services.OMEXMLServiceImpl();
    OMEXMLService.addModuloAlong(metadata, modlo, 0); 
end;
%%%%%%%%%%%%%%%%%%% set up Modulo XML description metadata if present - ends        
                
        % xml annotations
        xml_anno_cnt = 0;
        try
            %
            xmlfilenames = dir([folder filesep '*.xml']);                
            for k = 1 : numel(xmlfilenames)
            xmlfilename = xmlfilenames(k).name;
                fid = fopen([folder filesep xmlfilename],'r');
                fgetl(fid);
                description = fscanf(fid,'%c');
                fclose(fid);
                ind = modlo_cnt + k - 1;
                xml_anno_cnt = xml_anno_cnt + 1;
                metadata.setXMLAnnotationID(['Annotation:' num2str(ind)],ind);
                metadata.setXMLAnnotationValue(description,ind);              
                metadata.setImageAnnotationRef(['Annotation:' num2str(ind)],0,ind); 
            end                        
        catch err
            display(err.message);
        end       
        % xml annotations        
                
        prev_anno_cnt = modlo_cnt + xml_anno_cnt;
        
        % Tags
        tags_cnt = 0;        
        tags = ip.Results.Tags;    
        sepsec = ip.Results.Tags_SeparatingSeq;    
        if ~isempty(tags) && ~isempty(sepsec)
            tags = strsplit(tags,sepsec);
            for k=1:numel(tags)                
                ind = prev_anno_cnt + k - 1;                
                metadata.setTagAnnotationID(['Annotation:' num2str(ind)],k-1);
                metadata.setTagAnnotationValue(tags{k},k-1);              
                metadata.setImageAnnotationRef(['Annotation:' num2str(ind)], 0, ind);
                tags_cnt = tags_cnt + 1;
            end        
        end
        % Tags - end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % backup - use description field for xml annotations? 
        try
            dscr_acc = [];
            %
            xmlfilenames = dir([folder filesep '*.xml']);                
            for k = 1 : numel(xmlfilenames)
            xmlfilename = xmlfilenames(k).name;
                fid = fopen([folder filesep xmlfilename],'r');
                fgetl(fid);
                description = fscanf(fid,'%c');
                fclose(fid);
                %                  
                dscr_acc = [dscr_acc description];
            end                        
            %
        catch err
            display(err.message);
        end
        %
        if ~isempty(dscr_acc) 
            metadata.setImageDescription(sprintf('first line\nsecondline'), 0);
            metadata.setImageDescription(sprintf(dscr_acc),0);
        end 
        % backup - use description field for xml annotations? 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
        % Create ImageWriter
        writer = loci.formats.ImageWriter();
        writer.setWriteSequentially(true);
        writer.setMetadataRetrieve(metadata);        
        
        if ~isempty(ip.Results.Compression)
            writer.setCompression(ip.Results.Compression);
        end
        if ip.Results.BigTiff
            writer.getWriter(ometiffilename).setBigTiff(ip.Results.BigTiff);
        end
                        
        writer.setId(ometiffilename);

        % Load conversion tools for saving planes
        switch class(I)
            case {'int8', 'uint8'}
                getBytes = @(x) x(:);
            case {'uint16','int16'}
                getBytes = @(x) loci.common.DataTools.shortsToBytes(x(:), 0);
            case {'uint32','int32'}
                getBytes = @(x) loci.common.DataTools.intsToBytes(x(:), 0);
            case {'single'}
                getBytes = @(x) loci.common.DataTools.floatsToBytes(x(:), 0);
            case 'double'
                getBytes = @(x) loci.common.DataTools.doublesToBytes(x(:), 0);
        end
                
        ifd = loci.formats.tiff.IFD();
        ifd.putIFDValue(ifd.ROWS_PER_STRIP,sizeY);
        
        hw = [];
        if ~isempty(ip.Results.verbose) && ip.Results.verbose
            hw = waitbar(0, 'Loading planes...');
        end
end        

t0 = tic;

% writer.saveBytes(index-1, getBytes(I')); % slower :)
 
writer.getWriter(ometiffilename).saveBytes(index-1, getBytes(I'), ifd);
ifd.clear();
ifd.putIFDValue(ifd.ROWS_PER_STRIP,sizeY);

frame_time = toc(t0);

if ~isempty(hw)
    waitbar(index/final_index,hw); drawnow;
end
        
if index == final_index
    
    if ~isempty(hw)
        delete(hw); 
        drawnow;
    end
    
    writer.close();
    
    %xmlValidate = loci.formats.tools.XMLValidate();
    %comment = loci.formats.tiff.TiffParser(ometiffilename).getComment()
    %xmlValidate.process(ometiffilename, java.io.BufferedReader(java.io.StringReader(comment)));    
end;

end

function compressionTypes = getCompressionTypes()

% List all values of Compression
writer = loci.formats.ImageWriter();
compressionTypes = arrayfun(@char, writer.getCompressionTypes(),...
    'UniformOutput', false);
end

