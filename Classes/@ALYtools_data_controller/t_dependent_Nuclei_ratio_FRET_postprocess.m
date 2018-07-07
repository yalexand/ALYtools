function fig = t_dependent_Nuclei_ratio_FRET_postprocess(obj,in_fig,output_directory,~) 

     if ~isdir(output_directory), disp('wrong output directory, expect no output'), end
     
     fig = single(in_fig);   

     [sX,sY,sC,sZ,nFovs] = size(fig);

% for saving - don't override!!!     
        fname = obj.current_filename;
        fname = strrep(fname,'.OME.tiff','');
        fname = strrep(fname,'.OME.tif','');
        fname = strrep(fname,'.tif','');
          
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot intensity
d_vals = zeros(1,nFovs);
a_vals = zeros(1,nFovs);
for k=1:nFovs
    ud = fig(:,:,1,1,k);
    ua = fig(:,:,2,1,k);
    nukes = fig(:,:,3,1,k);
    d_sample = ud(nukes==1);
    d_sample = single(d_sample(:));
    a_sample = ua(nukes==1);
    a_sample = single(a_sample(:));
    d_vals(k) = median(d_sample(:));
    a_vals(k) = median(a_sample(:));
end
h = figure;
plot(1:nFovs,d_vals,'b.-',1:nFovs,a_vals,'r.-')
xlabel('frame #');
ylabel('intensity');
legend({'donor','acceptor'})
grid on;
ax_new=gca;
set(ax_new,'Position','default');
legend(ax_new,'-DynamicLegend');
saveName = [output_directory filesep fname '_segmented_intensities'];
saveas(h,saveName,'fig');
close(h);


% all Z calculateions are bound to the frame
NUCDATA = cell(1,nFovs);
E = 0.5;
for k=1:nFovs
    ud = single(fig(:,:,1,1,k));
    ua = single(fig(:,:,2,1,k));
    nukes = fig(:,:,3,1,k);
    %
%         a_data = ua(nukes==1);
%         d_data = ud(nukes==1);
%         IA_min = quantile(a_data(:),0.01);
%         ID_min = quantile(d_data(:),0.01);
%         Z = (ua - IA_min)./(ud - ID_min);
%             sample = Z(Z>0);
%             sample = sample(~isinf(sample));
%             sample = sample(~isnan(sample));
%             Zmin = quantile(sample(:),0.01);
%             Zmax = quantile(sample(:),0.99);        
        %            
        try        
            %
            ratrefXY = zeros(sX,sY);
            %
            L = bwlabel(nukes);
            stats = regionprops(L,'Centroid');
            nnucs = max(L(:));
            %
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
                %
                ratrefXY(L==n)=100*nuc_data(n,4); % A/D ratio multiplied by 100
            end
            % replace binary segmentation in output for the A/D ratio multiplied by 100
            fig(:,:,3,1,k)=ratrefXY;
        catch
            disp(['glitch at index ' num2str(k)]);
        end
        %
        if isempty(nuc_data), continue, end
        %
        k
        NUCDATA{k} = nuc_data;
end
%save('NUCDATA','NUCDATA'); %??

% MOVIE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fix NUCDATA :)
%%%%%%%%%%%%%%%%%%%%%%%%%%%% that is legacy; should't happen
Z = NUCDATA;
cell_nums = zeros(1,numel(Z));
for k=1:numel(Z)
    if isempty(Z{k})
        Z{k}=Z{k-1};        
    end
    cell_nums(k)=length(Z{k});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = (0:numel(NUCDATA)-1)*obj.t_dependent_Nuclei_ratio_FRET_TIMESTEP;

% save N(t) curve
xlswrite([output_directory filesep fname ' cell numbers curve'],[t' cell_nums']);

h=figure; 
movegui(h, 'onscreen');

vidObj = VideoWriter([output_directory filesep fname]);

vidObj.Quality = 100;
open(vidObj);

% get the figure and axes handles
 hFig = gcf;
 hAx  = gca;
 % set the figure to full screen
 set(hFig,'units','normalized','outerposition',[0 0 1 1]);
 % set the axes to full screen
 set(hAx,'Unit','normalized','Position',[0 0 1 1]);
 % hide the toolbar
 set(hFig,'menubar','none')
 % to hide the title
 set(hFig,'NumberTitle','off');

set(h,'Units','pixels'); 
rect = get(h,'Position');
rect(1:2) = [0 0]; 

mean_x = [];
std_x = [];
cell_num = [];

fontsize = 18;
markersize = 12;
linewidth = 2;

for k = 1 : numel(NUCDATA)
    
    try
        nuc_data = NUCDATA{k};
        nuc_size = nuc_data(:,1);
        corr_nuc_FRET_ratio = nuc_data(:,3);
        corr_nuc_FRET_ratio = corr_nuc_FRET_ratio(~isnan(corr_nuc_FRET_ratio));        
        nuc_FRET_ratio = nuc_data(:,4);
        %        
        subplot(2,2,1);    
        histogram(nuc_FRET_ratio,200,'Normalization','probability'); % should display pdf
        xlabel('FRET ratio','fontsize',fontsize);
        ylabel('PDF value','fontsize',fontsize);
        axis([0 2 0 0.1]);        
        set(gca,'FontSize',fontsize);
        title(['frame ' num2str(k)]);
            grid on;

        subplot(2,2,3);
        cell_num = [cell_num numel(nuc_size)];
        semilogy(t(1:numel(cell_num)),cell_num,'k.-','markersize',markersize,'linewidth',linewidth); %
        xlabel('time [h]','fontsize',fontsize);
        ylabel('cell number','fontsize',fontsize);
        axis([t(1) t(numel(cell_num)) min(cell_nums) max(cell_nums)]);
        set(gca,'FontSize',fontsize);        
            grid on;    

        subplot(2,2,4);
        mean_x = [mean_x mean(nuc_FRET_ratio)];
        std_x = [std_x std(nuc_FRET_ratio)];
        errorbar(t(1:numel(mean_x)),mean_x,std_x,'Color','red','Marker','o','MarkerFaceColor','green'); %
        %mseb(t(1:numel(mean_x)),mean_x,std_x); %
        xlabel('time [h]','fontsize',fontsize);
        ylabel('FRET ratio','fontsize',fontsize);
        set(gca,'FontSize',fontsize);        
            grid on;            
            
         subplot(2,2,2);
         nuc_FRET_ratio = nuc_FRET_ratio(~isnan(corr_nuc_FRET_ratio)); % because may be different size
         histogram2(nuc_FRET_ratio,corr_nuc_FRET_ratio,0:0.1:2,-1:0.1:1,'Normalization','probability','DisplayStyle','tile');
         xlabel('FRET ratio');
         ylabel('Donor-Acceptor pixel correlation (Pearson)','fontsize',fontsize);
         axis([0 2 -1 1]);
         set(gca,'FontSize',fontsize);
         title(strrep(fname,'_',' '));
             grid on;            
    catch
        disp(['glitch at ' num2str(k)]);
    end              
    movegui(h, 'onscreen');
    hold all;
    drawnow;
    writeVideo(vidObj,getframe(gcf,rect));
end
close(vidObj);  
close(h);
% MOVIE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TRACKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
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
        map.put('RADIUS', obj.TrackMate_RADIUS);
        map.put('TARGET_CHANNEL', obj.TrackMate_TARGET_CHANNEL);
        map.put('THRESHOLD', obj.TrackMate_THRESHOLD);
        map.put('DO_MEDIAN_FILTERING', obj.TrackMate_DO_MEDIAN_FILTERING);
        settings.detectorSettings = map;

        % Configure spot filters - Classical filter on quality
        filter1 = fiji.plugin.trackmate.features.FeatureFilter('QUALITY', obj.TrackMate_QUALITY, true);
        settings.addSpotFilter(filter1)

        % Configure tracker - We want to allow splits and fusions
        settings.trackerFactory  = fiji.plugin.trackmate.tracking.sparselap.SparseLAPTrackerFactory();
        settings.trackerSettings = fiji.plugin.trackmate.tracking.LAPUtils.getDefaultLAPSettingsMap(); % almost good enough
        settings.trackerSettings.put('ALLOW_TRACK_SPLITTING', obj.TrackMate_ALLOW_TRACK_SPLITTING);
        settings.trackerSettings.put('ALLOW_TRACK_MERGING', obj.TrackMate_ALLOW_TRACK_MERGING);

        % Configure track analyzers - Later on we want to filter out tracks 
        % based on their displacement, so we need to state that we want 
        % track displacement to be calculated. By default, out of the GUI, 
        % not features are calculated. 

        % The displacement feature is provided by the TrackDurationAnalyzer.
        settings.addTrackAnalyzer(fiji.plugin.trackmate.features.track.TrackDurationAnalyzer())

        % Configure track filters - We want to get rid of the two immobile spots at 
        % the bottom right of the image. Track displacement must be above 10 pixels.
        filter2 = fiji.plugin.trackmate.features.FeatureFilter('TRACK_DISPLACEMENT', obj.TrackMate_TRACK_DISPLACEMENT, true);
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

        % fname = [tempname,'.xml'];
        fullfname = [output_directory filesep fname '_RAW_TRACKMATE_OUTPUT' '.xml'];                        
        file = java.io.File(fullfname);
        fiji.plugin.trackmate.action.ExportTracksToXML.export(model, settings, file);        
        tracks = importTrackMateTracks(file,false); % use Z slot to store FRET ratio ?
                        
        % DEBUG
%         icyvol = fig;                
%         for k=1:numel(tracks)
%             track = tracks{k};
%             for m=1:size(track,1)
%                 frame = track(m,1)+1;
%                 x = round(track(m,2))+1;
%                 y = round(track(m,3))+1;
%                 try
%                 icyvol(x,y,3,1,frame) = 255;
%                 catch
%                     disp([x y frame]);
%                 end
%             end
%         end                
%         icy_imshow(uint16(icyvol));    

% matching
        for k=1:numel(tracks)
            track = tracks{k};
            for m=1:size(track,1)
                frame = track(m,1)+1;
                x = round(track(m,2))+1;
                y = round(track(m,3))+1;
                FRET_ratio = fig(x,y,3,1,frame);
                try
                if 0==FRET_ratio % safely use 5x5 vicinity of the orphan point
                    xl = max(x-2,1);
                    xr = min(x+2,sX);
                    yd = max(y-2,1);
                    yu = min(y+2,sY);
                    d_sample = fig(xl:xr,yd:yu,1,1,frame);
                    a_sample = fig(xl:xr,yd:yu,2,1,frame);
                    FRET_ratio = mean(a_sample(:))/mean(d_sample(:))*100;
                    %
                    disp(['fixing orphan TrackMate point at ' num2str(x) ' ' num2str(y) ' ' num2str(frame) ' ' num2str(FRET_ratio)]);
                    %                    
                end
                    track(m,4) = FRET_ratio/100;
                catch
                    disp([x y frame FRET_ratio/100]);
                end
            end
            tracks{k} = track;
        end 
fullfname = [output_directory filesep fname 'FRET_ratio_featured_TRACKMATE_OUTPUT'];
save(fullfname,'tracks');
% matching

% TRACKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fig = uint16(fig);                        
end
