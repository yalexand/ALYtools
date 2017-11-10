function frame_time = add_plane_to_OMEtiff_with_metadata(I, index, final_index, folder, ometiffilename, varargin)

global writer;
global getBytes;
global hw;
global ifd;

assert(isnumeric(I), 'First argument must be numeric');
assert(isnumeric(index), 'Second argument must be numeric');
assert(isnumeric(final_index), 'Third argument must be numeric');

ip = inputParser;

ip.addOptional('PhysSz_X', [], @isnumeric);
ip.addOptional('PhysSz_Y', [], @isnumeric);
ip.addOptional('ModuloZ_Type', [], @ischar);
ip.addOptional('ModuloZ_TypeDescription', [], @ischar);
ip.addOptional('ModuloZ_Unit', [], @ischar);
ip.addOptional('ModuloZ_Start', [], @isnumeric);
ip.addOptional('ModuloZ_Step', [], @isnumeric);
ip.addOptional('ModuloZ_End', [], @isnumeric);
ip.addOptional('ModuloZ_Labels', [], @isnumeric);
ip.addOptional('Tags', [], @ischar);    
ip.addOptional('Tags_SeparatingSeq', [], @ischar); 
ip.addOptional('verbose', [], @islogical);
ip.addParamValue('Compression', '',  @(x) ismember(x, getCompressionTypes()));
ip.addParamValue('BigTiff', false , @islogical);

ip.parse(varargin{:});

        [sizeX,sizeY] = size(I);
        sizeZ = final_index;

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

%%%%%%%%%%%%%%%%%%% set up Modulo XML description metadata if present - starts
if ~isempty(ip.Results.PhysSz_X)
    metadata.setPixelsPhysicalSizeX(ome.xml.model.primitives.PositiveFloat(java.lang.Double(ip.Results.PhysSz_X)),0);
end    
if ~isempty(ip.Results.PhysSz_Y)
    metadata.setPixelsPhysicalSizeY(ome.xml.model.primitives.PositiveFloat(java.lang.Double(ip.Results.PhysSz_Y)),0);
end

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

if exist('modlo','var') 
    OMEXMLService = loci.formats.services.OMEXMLServiceImpl();
    OMEXMLService.addModuloAlong(metadata, modlo, 0); 
end;
%%%%%%%%%%%%%%%%%%% set up Modulo XML description metadata if present - ends        

        % DESCRIPTION - one needs to find xml file if there... and so on
        n_anno = 0;  
        dscr_acc = [];
        if exist('modlo','var'), n_anno = n_anno + 1; end;        
        try
            %
            xmlfilenames = dir([folder filesep '*.xml']);                
            for k = 1 : numel(xmlfilenames)
            xmlfilename = xmlfilenames(k).name;
                fid = fopen([folder filesep xmlfilename],'r');
                fgetl(fid);
                description = fscanf(fid,'%c');
                fclose(fid);
                ind = n_anno+k-1;
                metadata.setXMLAnnotationID(['Annotation:' num2str(ind)],ind);
                metadata.setXMLAnnotationValue(description,ind);
                % helped by Jmarie, Sebastien - http://www.openmicroscopy.org/community/posting.php?mode=reply&f=6&t=7659
                metadata.setImageAnnotationRef(['Annotation:' num2str(ind)],0,ind); 
                %                  
                dscr_acc = [dscr_acc description];
            end                        
            %
            n_anno = ind;
            %
        catch err
            display(err.message);
        end
        
        % backup - use description field? 
        if ~isempty(dscr_acc) 
            metadata.setImageDescription(sprintf('first line\nsecondline'), 0);
            metadata.setImageDescription(sprintf(dscr_acc),0);
        end        
        %
        % DESCRIPTION - ends
        
        % Tags
        tags = ip.Results.Tags;    
        sepsec = ip.Results.Tags_SeparatingSeq;    
        if ~isempty(tags) && ~isempty(sepsec)
            tags = strsplit(tags,sepsec);
            for k=1:numel(tags)
                metadata.setTagAnnotationID(['Annotation:' num2str(n_anno+k)],k-1);
                metadata.setTagAnnotationValue(tags{k},k-1);
                metadata.setImageAnnotationRef(['Annotation:' num2str(n_anno+k)], 0, n_anno+k);
            end        
        end
        % Tags - end
                                       
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

