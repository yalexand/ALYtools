function [datas, captions, table_names, fig] = analyze_NucleiTimeStack(obj,~,~) 

     datas = [];
     captions = [];
     table_names = 'default';
     fig = [];
     
     fig = obj.do_NucleiTimeStack_Segmentation(false);

     [sX,sY,sC,sZ,nFovs] = size(fig);
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first define autofluorescence
d_data = [];
a_data = [];

% d_vals = zeros(1,nFovs);
% a_vals = zeros(1,nFovs);
% for k=1:nFovs
%     ud = fig(:,:,1,1,k);
%     ua = fig(:,:,2,1,k);
%     nukes = fig(:,:,3,1,k);
%     d_sample = ud(nukes==1);
%     d_sample = single(d_sample(:));
%     a_sample = ua(nukes==1);
%     a_sample = single(a_sample(:));
%     d_vals(k) = median(d_sample(:));
%     a_vals(k) = median(a_sample(:));
% end
% figure(22);
% plot(1:nFovs,d_vals,'b.-',1:nFovs,a_vals,'r.-')
% xlabel('frame #');
% ylabel('intensity');
% legend({'donor','acceptor'})
% grid on;

for k=1:nFovs
    ud = fig(:,:,1,1,k);
    ua = fig(:,:,2,1,k);
    nukes = fig(:,:,3,1,k);
    d_sample = ud(nukes==1);
    d_sample = single(d_sample(:));
    a_sample = ua(nukes==1);
    a_sample = single(a_sample(:));
    d_data = [d_data d_sample'];
    a_data = [a_data a_sample'];
    k
end
%
IA_min = quantile(a_data,0.01)
IA_max = quantile(a_data,0.99)
%
ID_min = quantile(d_data,0.01)
ID_max = quantile(d_data,0.99)
%
% boundaries for Z
Z_data = [];
for k=1:nFovs
    ud = fig(:,:,1,1,k);
    ua = fig(:,:,2,1,k);
    nukes = fig(:,:,3,1,k);
    d_sample = ud(nukes==1);
    d_sample = single(d_sample(:));
    a_sample = ua(nukes==1);
    a_sample = single(a_sample(:));
    Z = (a_sample - IA_min)./(d_sample - ID_min);
    Z_data = [d_data Z'];
    k    
end

Z = Z(Z>0);
Z = Z(~isinf(Z));

%Zmin = min(Z(:));
%Zmax = max(Z(:));
Zmin = quantile(Z(:),0.01);
Zmax = quantile(Z(:),0.99);

%nFovs = 6;
NUCDATA = cell(1,nFovs);
E = 0.5;
for k=1:nFovs
    ud = fig(:,:,1,1,k);
    ua = fig(:,:,2,1,k);
    nukes = fig(:,:,3,1,k);
    %
        nnucs = 0;
        nuc_data = [];
        %
        Z = (ua - IA_min)./(ud - ID_min);
        %
        try        
            %        
            L = bwlabel(nukes);
            nnucs = max(L(:));
            % should be XY + RATIO + RAW INTENSITIES  = +4 COLUMNS
            nuc_data = zeros(nnucs,3); 
            %                        
            for n=1:nnucs
                sample_a = ua(L==n);
                sample_d = ud(L==n);
                Z_cur = Z(L==n);
                Z_cur = Z_cur(~isnan(Z_cur));
                Z_cur(Z_cur<Zmin)=Zmin;
                Z_cur(Z_cur>Zmax)=Zmax;                
                nuc_data(n,1)=length(sample_a(:)); % area
                z = mean(Z_cur(:));
                x = (z-Zmin)/(Zmax-Zmin-E*(Zmax-z));
                nuc_data(n,2)=x;
                nuc_data(n,3)=corr(sample_a(:),sample_d(:),'type','Pearson');
            end
        catch
            disp(['glitch at index ' num2str(k)]);
        end
        %
        if isempty(nuc_data), continue, end
        %
        k
        NUCDATA{k} = nuc_data;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
     dum = 255*fig(:,:,3,1,:);
     imp = copytoImagePlus(dum,'dimorder','XYCZT');
     imp.show();     
     % test stuff     
     %imp = ij.IJ.openImage('http://fiji.sc/samples/FakeTracks.tif');
     %imp.show();
     %
        % Some of the parameters we configure below need to have
        % a reference to the model at creation. So we create an
        % empty model now.
        model = fiji.plugin.trackmate.Model();

        % Send all messages to ImageJ log window.
        % model.setLogger(fiji.plugin.trackmate.Logger.IJ_LOGGER)

        %------------------------
        % Prepare settings object
        %------------------------

        settings = fiji.plugin.trackmate.Settings();
        settings.setFrom(imp)

        % Configure detector - We use a java map
        settings.detectorFactory = fiji.plugin.trackmate.detection.LogDetectorFactory();
        map = java.util.HashMap();
        map.put('DO_SUBPIXEL_LOCALIZATION', true);
        map.put('RADIUS', 6);
        map.put('TARGET_CHANNEL', 1);
        map.put('THRESHOLD', 0);
        map.put('DO_MEDIAN_FILTERING', false);
        settings.detectorSettings = map;

        % Configure spot filters - Classical filter on quality
        filter1 = fiji.plugin.trackmate.features.FeatureFilter('QUALITY', 1.0, true);
        settings.addSpotFilter(filter1)

        % Configure tracker - We want to allow splits and fusions
        settings.trackerFactory  = fiji.plugin.trackmate.tracking.sparselap.SparseLAPTrackerFactory();
        settings.trackerSettings = fiji.plugin.trackmate.tracking.LAPUtils.getDefaultLAPSettingsMap(); % almost good enough
        settings.trackerSettings.put('ALLOW_TRACK_SPLITTING', true);
        settings.trackerSettings.put('ALLOW_TRACK_MERGING', true);

        % Configure track analyzers - Later on we want to filter out tracks 
        % based on their displacement, so we need to state that we want 
        % track displacement to be calculated. By default, out of the GUI, 
        % not features are calculated. 

        % The displacement feature is provided by the TrackDurationAnalyzer.
        settings.addTrackAnalyzer(fiji.plugin.trackmate.features.track.TrackDurationAnalyzer())

        % Configure track filters - We want to get rid of the two immobile spots at 
        % the bottom right of the image. Track displacement must be above 10 pixels.
        filter2 = fiji.plugin.trackmate.features.FeatureFilter('TRACK_DISPLACEMENT', 2.0, true);
        settings.addTrackFilter(filter2)

        %-------------------
        % Instantiate plugin
        %-------------------
        trackmate = fiji.plugin.trackmate.TrackMate(model, settings);

        %--------
        % Process
        %--------
        ok = trackmate.checkInput();
        if ~ok
            display(trackmate.getErrorMessage())
        end

        ok = trackmate.process();
        if ~ok
            display(trackmate.getErrorMessage())
        else
            % Echo results
            display(model.toString());
        end
        
        %----------------
        % Display results
        %----------------
        selectionModel = fiji.plugin.trackmate.SelectionModel(model);
        displayer =  fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer(model, selectionModel, imp);
        %displayer =  fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer(model, selectionModel, copytoImagePlus(dum,'dimorder','XYZCT'));        
        displayer.render()
        displayer.refresh()
   
 %--------------------------------------------------------------------------------------------------------------------------------
            function imp = copytoImagePlus(I,varargin) % for now keep it here
                    % copytoImagePlus allows you to open an array I with an instance of ImageJ
                    % within MATLAB with a proper data type and hyperstack dimensions.
                    %
                    %
                    % SYNTAX
                    % imp = copytoImagePlus(I)
                    % imp = copytoImagePlus(I,dimorder)
                    % imp = copytoImagePlus(____,'Param',value)
                    % 
                    %
                    % REQUIREMENTS
                    % ImageJ-MATLAB as part of Fiji installation
                    % https://imagej.net/MATLAB_Scripting
                    %
                    % ijmshow assumes a net.imagej.matlab.ImageJMATLABCommands Java object
                    % named 'IJM' is made available in the base Workspace by ImageJ (part of
                    % ImageJ-MATLAB).
                    %
                    %
                    % INPUT ARGUMENTS
                    % I           uint16 | uint8 | double | single
                    %             An array of integers to be opened with ImageJ. This array can
                    %             have from 2 to 5 dimensions.
                    %
                    %
                    % dimorder    char row vector made of 'XYCZT' | 'YXCZT' (default)
                    %
                    %             (Optional) A char row vector composed of 'X', 'Y', 'C' for
                    %             channels, 'Z' for slices, and 'T' for frames. dimorder is
                    %             case insensitive. You cannot repeat any of the five letters
                    %             in dimorder. The first two letters must be either 'X' or 'Y'.
                    %             The length of dimorder must be 5 or match the number of
                    %             dimensions of the array specified by I. The third to the
                    %             fifth letters must be chosen from 'C', 'Z', and 'T'.
                    %
                    %             The default is set 'YXCZT' rather than 'XYZCT', because the X
                    %             and Y axes of an MATLAB array is flipped over in ImageJ by
                    %             IJM.show().
                    %
                    %
                    % OPTIONAL PARAMETER/VALUE PAIRS
                    % NewName     char row vector | 'new' (default)
                    %             The window title of the new image in ImageJ
                    %
                    % FrameInterval
                    %             scalar
                    %             Time frame sampling interval in seconds
                    %
                    % OUTPUT ARGUMENTS
                    % imp         ij.ImagePlus Java object
                    %
                    % EXAMPLES
                    % see https://github.com/kouichi-c-nakamura/copytoImagePlus
                    %
                    %
                    % Written by Kouichi C. Nakamura Ph.D.
                    % MRC Brain Network Dynamics Unit
                    % University of Oxford
                    % kouichi.c.nakamura@gmail.com
                    % 03-May-2018 04:57:24
                    %
                    % See also
                    % ijmshow (this requires a net.imagej.matlab.ImageJMATLABCommands object IJM)
                    % https://github.com/kouichi-c-nakamura/ijmshow (repository for this function)
                    %
                    % ImageJ as part of ImageJ-MATLAB (https://github.com/imagej/imagej-matlab/)
                    %
                    % net.imagej.matlab.ImageJMATLABCommands
                    % evalin, assignin
                    % https://imagej.net/MATLAB_Scripting


                    import ij.process.ShortProcessor
                    import ij.process.ByteProcessor
                    import ij.process.FloatProcessor

                    p = inputParser;
                    p.addRequired('I',@(x) isnumeric(x));
                    p.addOptional('dimorder','YXCZT',@(x) ischar(x) && isrow(x) ...
                        && all(arrayfun(@(y) ismember(y,'XYCZT'),upper(x))) && length(x) >=2 ...
                        && all(arrayfun(@(y) ismember(y,'XY'),upper(x(1:2))))...
                        );
                    p.addParameter('NewName','new',@(x) ischar(x) && isrow(x));
                    p.addParameter('FrameInterval',[],@(x) isreal(x) && x > 0);

                    p.parse(I,varargin{:});

                    dimorder = upper(p.Results.dimorder);
                    newname = p.Results.NewName;
                    frameinterval = p.Results.FrameInterval;



                    switch dimorder(1:2)
                        case 'XY'
                            order1 = [1 2];
                        case 'YX'
                            order1 = [2 1];
                    end


                    switch dimorder(3:ndims(I))
                        case 'CZT'
                            order2 = 3:5;
                        case 'CTZ'
                            order2 = [3 5 4];
                        case 'ZCT'
                            order2 = [4 3 5];
                        case 'ZTC'
                            order2 = [4 5 3];
                        case 'TCZ'
                            order2 = [5 3 4];
                        case 'TZC'
                            order2 = [5 4 3];
                        case 'CZ'
                            order2 = [3 4];
                        case 'CT'
                            order2 = [3 5 4];
                        case 'ZC'
                            order2 = [4 3];
                        case 'ZT'
                            order2 = [5 3 4];        
                        case 'TC'
                            order2 = [4 5 3];
                        case 'TZ'
                            order2 = [5 4 3];
                        case 'C'
                            order2 = [3 4 5];
                        case 'Z'
                            order2 = [4 3 5];
                        case 'T'
                            order2 = [4 5 3];
                        otherwise
                            order2 = 3:5;
                    end


                    I0 = permute(I, [order1, order2]);

                    nX = int32(size(I0,1));
                    nY = int32(size(I0,2));
                    nC = int32(size(I0,3)); 
                    nZ = int32(size(I0,4));
                    nT = int32(size(I0,5));

                    try
                        switch class(I0)
                            case 'uint8'
                                bitdepth = 8;

                            case 'int8'

                                bitdepth = 8;

                            case 'uint16'

                                bitdepth = 16;

                            case 'int16'

                                bitdepth = 8;

                            case 'uint32'

                                bitdepth = 32;

                            case 'int32'

                                bitdepth = 32;

                            case 'uint64'
                                error('MATLAB:copytoImg:UnsupportedType', ...
                                    'uint64 is not supported.');
                            case 'int64'
                                error('MATLAB:copytoImg:UnsupportedType', ...
                                    'uint64 is not supported.');
                            case 'single'

                                bitdepth = 32;

                            case 'double'

                                bitdepth = 32;

                            case 'logical'

                                bitdepth = 8;

                            otherwise
                                error('MATLAB:copytoImg:UnsupportedType', ...
                                    '%s is not supported.', class(I0));

                        end

                    catch merr
                        if strcmp(merr.identifier, 'MATLAB:undefinedVarOrClass')
                                error('MATLAB:copytoImg:undefinedVarOrClass', ...
                                    'Could not find ImgLib2 on the path. Did you forget to run ''Miji(false)'' before calling this function?');
                        else
                            rethrow(merr);
                        end
                    end


                    imp = ij.IJ.createHyperStack(newname,nX,nY,nC,nZ,nT,bitdepth);

                    for t = 1:nT
                        imp.setT(t);
                        for z = 1:nZ
                            imp.setZ(z);
                            for c = 1:nC
                                imp.setC(c);

                                XY = I0(:,:,c,z,t);
                                xy = XY(:)';

                                switch bitdepth
                                    case 16
                                        ip = ShortProcessor(nX,nY);
                                        ip.setPixels(xy);
                                        imp.setProcessor(ip);
                                    case 8
                                        ip = ByteProcessor(nX,nY);
                                        ip.setPixels(xy);
                                        imp.setProcessor(ip);
                                    otherwise
                                        ip = FloatProcessor(nX,nY);
                                        ip.setPixels(single(xy));                  
                                        imp.setProcessor(ip);
                                end
                            end
                        end
                    end

                    imp.setT(1);
                    imp.setZ(1);
                    imp.setC(1);
                    %imp.show();
                    imp.setDisplayMode(ij.IJ.COLOR) %NOTE this is required to enable the next line
                    imp.setDisplayMode(ij.IJ.COMPOSITE)

                    try
                        imp.resetDisplayRanges();
                    catch mexc
                        if strcmpi(mexc.identifier,'MATLAB:UndefinedFunction')
                            warning('resetDisplayRanges did not work')
                        else
                            throw(mexc)
                        end
                    end

                    if ~isempty(frameinterval)

                        fi = imp.getFileInfo();
                        fi.frameInterval = frameinterval;
                        imp.setFileInfo(fi);
                        %TODO Show Info... does not show the frameinterval
                    end
            end

end
