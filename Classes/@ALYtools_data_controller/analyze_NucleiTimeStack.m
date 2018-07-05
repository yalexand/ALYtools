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

% all Z calculateions are bound to the frame
NUCDATA = cell(1,nFovs);
E = 0.5;
for k=1:nFovs
    ud = fig(:,:,1,1,k);
    ua = fig(:,:,2,1,k);
    nukes = fig(:,:,3,1,k);
    %
        a_data = ua(nukes==1);
        d_data = ud(nukes==1);
        IA_min = quantile(a_data(:),0.01);
        ID_min = quantile(d_data(:),0.01);
        Z = (ua - IA_min)./(ud - ID_min);
            sample = Z(Z>0);
            sample = sample(~isinf(sample));
            sample = sample(~isnan(sample));
            Zmin = quantile(sample(:),0.01);
            Zmax = quantile(sample(:),0.99);        
        %            
        try        
            %        
            L = bwlabel(nukes);
            stats = regionprops(L,'Centroid');
            nnucs = max(L(:));
            % should be XY + RATIO + RAW INTENSITIES  = +4 COLUMNS
            nuc_data = zeros(nnucs,8); 
            for n=1:nnucs
                sample_a = single(ua(L==n));
                sample_d = single(ud(L==n));
                %
%                 Z_cur = Z(L==n);
%                 Z_cur = Z_cur(~isnan(Z_cur));
%                 Z_cur(Z_cur<Zmin)=Zmin;
%                 Z_cur(Z_cur>Zmax)=Zmax;                
%                 nuc_data(n,1)=length(sample_a(:)); % area
%                 z_cur = mean(Z_cur(:));                
%                 nuc_data(n,2)=(z_cur-Zmin)/(Zmax-Zmin-E*(Zmax-z_cur));
                nuc_data(n,2)=0;
                %
                nuc_data(n,3)=corr(sample_a(:),sample_d(:),'type','Pearson');
                    mean_sample_a = mean(sample_a(:));
                    mean_sample_d = mean(sample_d(:));
                nuc_data(n,4)=mean_sample_a/mean_sample_d;
                nuc_data(n,5)=mean_sample_a;
                nuc_data(n,6)=mean_sample_d;
                nuc_data(n,7) = stats(n).Centroid(2);
                nuc_data(n,8) = stats(n).Centroid(1);
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
%save('NUCDATA','NUCDATA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
     dum = 255*fig(:,:,3,1,:);
     imp = copytoImagePlus(dum,'dimorder','XYCZT');
        % test image
        % imp = ij.IJ.openImage('http://fiji.sc/samples/FakeTracks.tif');
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
        map.put('RADIUS', 4);
        map.put('TARGET_CHANNEL', 1);
        map.put('THRESHOLD', 0);
        map.put('DO_MEDIAN_FILTERING', false);
        settings.detectorSettings = map;

        % Configure spot filters - Classical filter on quality
        filter1 = fiji.plugin.trackmate.features.FeatureFilter('QUALITY', 0.5, true);
        settings.addSpotFilter(filter1)

        % Configure tracker - We want to allow splits and fusions
        settings.trackerFactory  = fiji.plugin.trackmate.tracking.sparselap.SparseLAPTrackerFactory();
        settings.trackerSettings = fiji.plugin.trackmate.tracking.LAPUtils.getDefaultLAPSettingsMap(); % almost good enough
        settings.trackerSettings.put('ALLOW_TRACK_SPLITTING', false);
        settings.trackerSettings.put('ALLOW_TRACK_MERGING', false);

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
%         selectionModel = fiji.plugin.trackmate.SelectionModel(model);
%         displayer =  fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer(model, selectionModel, imp);
%         %displayer =  fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer(model, selectionModel, copytoImagePlus(dum,'dimorder','XYZCT'));        
%         displayer.render()
%         displayer.refresh()

%         fname = obj.current_filename;
%         fname = strrep(fname,'.OME.tiff','');
%         fname = strrep(fname,'.OME.tif','');
%         fname = strrep(fname,'.tif','');
%         fname = [obj.RootDirectory filesep fname '.xml'];
        %
        fname = [tempname,'.xml'];
        file = java.io.File(fname);
        fiji.plugin.trackmate.action.ExportTracksToXML.export(model, settings, file);        
        tracks = importTrackMateTracks(file,true);

        % DEBUG
        icyvol = fig;                
        for k=1:numel(tracks)
            track = tracks{k};
            for m=1:size(track,1)
                frame = track(m,1)+1;
                x = round(track(m,2))+1;
                y = round(track(m,3))+1;
                try
                icyvol(x,y,3,1,frame) = 255;
                catch
                    disp([x y frame]);
                end
            end
        end        
        icy_imshow(uint16(icyvol));    
end
