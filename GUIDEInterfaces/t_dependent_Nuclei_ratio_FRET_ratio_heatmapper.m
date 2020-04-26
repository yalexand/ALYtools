function varargout = t_dependent_Nuclei_ratio_FRET_ratio_heatmapper(varargin)
% T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER MATLAB code for t_dependent_Nuclei_ratio_FRET_ratio_heatmapper.fig
%      T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER, by itself, creates a new T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER or raises the existing
%      singleton*.
%
%      H = T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER returns the handle to a new T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER or the handle to
%      the existing singleton*.
%
%      T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER.M with the given input arguments.
%
%      T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER('Property','Value',...) creates a new T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before t_dependent_Nuclei_ratio_FRET_ratio_heatmapper_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to t_dependent_Nuclei_ratio_FRET_ratio_heatmapper_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help t_dependent_Nuclei_ratio_FRET_ratio_heatmapper

% Last Modified by GUIDE v2.5 26-Apr-2020 11:34:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @t_dependent_Nuclei_ratio_FRET_ratio_heatmapper_OpeningFcn, ...
                   'gui_OutputFcn',  @t_dependent_Nuclei_ratio_FRET_ratio_heatmapper_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before t_dependent_Nuclei_ratio_FRET_ratio_heatmapper is made visible.
function t_dependent_Nuclei_ratio_FRET_ratio_heatmapper_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to t_dependent_Nuclei_ratio_FRET_ratio_heatmapper (see VARARGIN)

% Choose default command line output for t_dependent_Nuclei_ratio_FRET_ratio_heatmapper
handles.output = hObject;

handles.TrackPlotter_handles = varargin{1};

[~,NAME,EXT]=fileparts(handles.TrackPlotter_handles.fullfilename);
set(handles.figure1, 'Name', ['FRET ratio heatmap : ' NAME EXT]);

set(handles.features_chooser,'String',handles.TrackPlotter_handles.features);
set(handles.features_chooser,'Value',6);
set(handles.curve_type,'String',{'mean','median'});
set(handles.n_clusters,'String',{'2','3','4','5'});
set(handles.clustering_method,'String',{'kmeans','kmedoids','hierarchical'});

handles.CRV_range = [];     %2 x Ncurv array
handles.CRV_type = [];      % for simplicity: 1(mean), 2(median)
handles.CRV_valarr = [];   % to check bounds
handles.CRV_add = false;    % default - not adding curves

% Update handles structure
guidata(hObject, handles);

visualize(hObject,handles);

% UIWAIT makes t_dependent_Nuclei_ratio_FRET_ratio_heatmapper wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = t_dependent_Nuclei_ratio_FRET_ratio_heatmapper_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in features_chooser.
function features_chooser_Callback(hObject, eventdata, handles)
% hObject    handle to features_chooser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns features_chooser contents as cell array
%        contents{get(hObject,'Value')} returns selected item from features_chooser
handles.CRV_range = [];
handles.CRV_type = [];
handles.CRV_valarr = [];
handles.CRV_add = false;
guidata(hObject, handles);
visualize(hObject,handles);

% --- Executes during object creation, after setting all properties.
function features_chooser_CreateFcn(hObject, eventdata, handles)
% hObject    handle to features_chooser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ------------------------------------------------------
function visualize(hObject,handles)

    % to update mask
    fullfilename = handles.TrackPlotter_handles.fullfilename;
    FIGS = findobj(0, 'type', 'figure');
    for k=1:size(FIGS,1)
        h = FIGS(k,1);
        cur_handles = guidata(h);
        try
        if strcmp(cur_handles.fullfilename,fullfilename) % mine
            handles.TrackPlotter_handles = cur_handles;
            guidata(hObject, handles); 
            break;
        end
        catch
        end
    end
    % to update mask

    tracks = handles.TrackPlotter_handles.MI_tracks;
    display_tracks = handles.TrackPlotter_handles.MI_norm_FRET_ratio;
    
    if 1==get(handles.selection_only,'Value') % if masking option is on
        mask = handles.TrackPlotter_handles.mask;
        tracks = tracks(mask);
        display_tracks = display_tracks(mask,:);
    end
    
    track = tracks{1};
    
    par_ind = get(handles.features_chooser,'Value');        
              index = 0;
                 switch par_ind                    
                    case 4 % #nghbrs                    
                        if ismember(size(track,2),[10 11 15])
                            index = 9;
                        end
                    case 5 % cell density
                        if ismember(size(track,2),[10 11 15])
                            index = 10;
                        end
                    case 6 % FRET ratio
                            index = 4;
                    case 8 % donor intensity                     
                        index = 5;
                    case 9 % acceptor intensity
                        index = 6;
                    case 10 % nucleus size
                        index = 7;
                    case 11 % Pearson                      
                        index = 8;
                    case 14 % FRET molar fraction
                        index = 11;
                    case 15 % cell size
                        if 15== size(track,2)
                        index = 12;
                        end
                    case 16 % 
                       if 15== size(track,2)
                        index = 13;
                        end
                    case 17 % 
                        if 15== size(track,2)
                        index = 14;
                        end
                    case 18 % 
                         if 15== size(track,2)
                        index = 15;
                        end
                 end                  
   
                 if 0==index, return, end
                 
                 valarr = [];
                 for k=1:numel(tracks)
                    track = tracks{k};
                    vals = squeeze(track(:,index));
                    valarr = [valarr; median(vals)];
                 end

     handles.CRV_valarr = valarr; % not sorted, essential
     guidata(hObject, handles);
                 
     [~,I] = sort(valarr,'descend');
     M = display_tracks(I,:);
        
%     import bioma.data.*
%     M=map(M,0,1);
%     DMobj = DataMatrix(M);
%     hmo = HeatMap(DMobj,'Colormap',redgreencmap);
    imagesc(handles.heatmap_axes,M);
    cmap = jet(256);       
    cmap(1,:)=[0,0,0];
    colormap(handles.heatmap_axes,cmap);
    %   
    qntvals = [0 .2 .4 .6 .8 1];
    quantiles = quantile(handles.CRV_valarr,qntvals);
    qntsteps = ceil(length(handles.CRV_valarr)*qntvals);
    qntsteps(1)=1;
    if numel(unique(qntsteps)) < numel(qntvals) % switch to range
        quantiles = [min(handles.CRV_valarr) max(handles.CRV_valarr)];
        qntsteps = [ 1 length(handles.CRV_valarr)];
    end
    yticks(handles.heatmap_axes,qntsteps);
    yticklabels(handles.heatmap_axes,flip(quantiles));
        a = get(handles.heatmap_axes,'YTickLabel');
        set(handles.heatmap_axes,'YTickLabel',a,'fontsize',8)
    str = get(handles.features_chooser,'String'); 
    ylabel(handles.heatmap_axes,str{par_ind},'fontsize',8);
    %
    t = round((0:size(display_tracks,2)-1)*handles.TrackPlotter_handles.dt*60);
    %
    b1=min(t);
    b2=max(t);
    step = fix(length(t)/4);
    xticks(handles.heatmap_axes,0:step:(length(t)));
    step_val = round((b2 - b1)/4);
    xticklabels(handles.heatmap_axes,b1+step_val*(0:4));
    %
    xlabel(handles.heatmap_axes,'time [min]');
    grid(handles.heatmap_axes,'on'); 
    
    % 2 default curves
    if ~ handles.CRV_add % default
        mode = get(handles.curve_type,'Value'); % either 1(mean) or 2(median)
        if isempty(handles.CRV_range)
            handles.CRV_range = [handles.CRV_range; [min(valarr) median(valarr)]];
            handles.CRV_type = [handles.CRV_type; mode];
            handles.CRV_range = [handles.CRV_range; [median(valarr) max(valarr)]];
            handles.CRV_type = [handles.CRV_type; mode];        
        end
        guidata(hObject, handles);
    end
    %
    plot_curves(handles,display_tracks)
    
%%%% 2 default curves
% % % % %    cla(handles.curves_axes,'reset');
% % % % %    N = size(M,1);
% % % % % if true && N>=2
% % % % %     
% % % % %    m1 = M(1:fix(N/2),:);
% % % % %    m2 = M(fix(N/2):N,:);
% % % % %    
% % % % %    mode_ind = get(handles.curve_type,'Value');
% % % % %    mode_str = get(handles.curve_type,'String');
% % % % %    if 1 == mode_ind % 'mean'
% % % % %        curve1 = mean(m1,1);
% % % % %        curve2 = mean(m2,1);
% % % % %    else
% % % % %        curve1 = median(m1,1);
% % % % %        curve2 = median(m2,1);       
% % % % %    end
% % % % %    
% % % % % t = round((0:length(curve1)-1)*handles.TrackPlotter_handles.dt*60);
% % % % % plot(handles.curves_axes,t,curve1,'k--',t,curve2,'k:','linewidth',2);
% % % % % legend(handles.curves_axes,{['upper half: ' mode_str{mode_ind}],['lower half: ' mode_str{mode_ind}]});
% % % % % xlabel(handles.curves_axes,'time [min]');
% % % % % set(handles.curves_axes,'yticklabel', []);
% % % % % axis(handles.curves_axes,[min(t) max(t) min(min(curve1(:)),min(curve2(:))) max(max(curve1(:)),max(curve2(:)))]);
% % % % % grid(handles.curves_axes,'on'); 
% % % % % 
% % % % % str = get(handles.features_chooser,'String'); 
% % % % % ylabel(handles.heatmap_axes,str{par_ind});
% % % % % %
% % % % % sarr = 0.1*round(handles.CRV_valarr*10);
% % % % % b1=min(sarr(:));
% % % % % b2=max(sarr(:));
% % % % % step = min(length(sarr),fix(length(sarr)/10));
% % % % % yticks(handles.heatmap_axes,0:step:(length(sarr)));
% % % % % step_val = (b2 - b1)/10;
% % % % % yticklabels(handles.heatmap_axes,flip(b1+step_val*(0:10)));
% % % % % %
% % % % % b1=min(t);
% % % % % b2=max(t);
% % % % % step = fix(length(t)/4);
% % % % % xticks(handles.heatmap_axes,0:step:(length(t)));
% % % % % step_val = round((b2 - b1)/4);
% % % % % xticklabels(handles.heatmap_axes,b1+step_val*(0:4));
% % % % % %
% % % % % xlabel(handles.heatmap_axes,'time [min]');
% % % % % grid(handles.heatmap_axes,'on'); 
% % % % % end

function plot_curves(handles,display_tracks)

    if length(handles.CRV_type)> 7, return, end % too many curves, c'mon

    colors = zeros(7,3);
    colors(1,:) = [0 0.4470 0.7410];
    colors(2,:) = [0.8500 0.3250 0.0980];
    colors(3,:) = [0.9290 0.6940 0.1250];
    colors(4,:) = [0.4940 0.1840 0.5560];
    colors(5,:) = [0.4660 0.6740 0.1880];
    colors(6,:) = [0.3010 0.7450 0.9330];
    colors(7,:) = [0.6350 0.0780 0.1840];
    
    styles = {'-',':','--','-.','-',':','--'};
    
    t = round((0:size(display_tracks,2)-1)*handles.TrackPlotter_handles.dt*60);

    cla(handles.curves_axes,'reset'); 
    maxval = -inf;
    minval = inf;
    LEGEND = cell(1,0);
    for k=1:length(handles.CRV_type)
        hold(handles.curves_axes,'on');
        mode = handles.CRV_type(k);
        b1 = handles.CRV_range(k,1);
        b2 = handles.CRV_range(k,2);
        I = handles.CRV_valarr>=b1 & handles.CRV_valarr<=b2;
        s = display_tracks(I,:);
        switch mode
            case 1 
                curve = mean(s,1);
                type_lab = 'mean ';
            case 2 
                curve = median(s,1);
                type_lab = 'median ';
        end
        maxval = max(maxval,max(curve));
        minval = min(minval,min(curve));
        plot(handles.curves_axes,t,curve,'linestyle',styles{k},'Color',colors(k,:),'linewidth',2);
        LEGEND = [LEGEND cellstr([type_lab num2str(b1) ' : ' num2str(b2)])];
    end
    hold(handles.curves_axes,'off');
    xlabel(handles.curves_axes,'time [min]');
    set(handles.curves_axes,'yticklabel', []);
    axis(handles.curves_axes,[min(t) max(t) minval maxval]);
    grid(handles.curves_axes,'on'); 
    ylabel(handles.curves_axes,'normalized FRET ratio');
    legend(handles.curves_axes,LEGEND);
    %
    str = get(handles.features_chooser,'String'); 
    par_ind = get(handles.features_chooser,'Value');
    title(handles.curves_axes,['partitioned by ' str{par_ind}],'fontsize',8);

% --- Executes on button press in selection_only.
function selection_only_Callback(hObject, eventdata, handles)
    handles.CRV_range = [];
    handles.CRV_type = [];
    handles.CRV_valarr = [];
    handles.CRV_add = false; % default
    guidata(hObject, handles);
    visualize(hObject,handles);

% --- Executes on button press in add_curve.
function add_curve_Callback(hObject, eventdata, handles)
    range = get_CRV_range(get(handles.sorting_parameter_range,'String'),handles);
    if isempty(range), return, end
    %
    if ~handles.CRV_add % if default
        handles.CRV_range = [];
        handles.CRV_type = [];
        handles.CRV_add = true;
        guidata(hObject, handles);
        cla(handles.curves_axes,'reset');
    end
    %
    curve_type = get(handles.curve_type,'Value'); % 1(mean) 2(median)
    handles.CRV_range = [handles.CRV_range; range];
    handles.CRV_type = [handles.CRV_type; curve_type];
        set(handles.sorting_parameter_range,'String','');
        guidata(hObject, handles);
    visualize(hObject,handles);
    handles.CRV_range
    handles.CRV_type

function sorting_parameter_range_Callback(hObject, eventdata, handles)
    if isempty(get_CRV_range(get(hObject,'String'),handles))
        set(hObject,'String','');
        guidata(hObject, handles);
    end

function range = get_CRV_range(str,handles)
    range = [];
    s = strsplit(str);
    if 2~=numel(s), return, end
        b1 = str2double(s{1});
        b2 = str2double(s{2});
        if isnan(b1) || isnan(b2), return, end 
            minv = min(handles.CRV_valarr);
            maxv = max(handles.CRV_valarr);
                if b1>maxv || b2<minv, return, end
                    range = [b1 b2];


% --- Executes during object creation, after setting all properties.
function sorting_parameter_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sorting_parameter_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in curve_type.
function curve_type_Callback(hObject, eventdata, handles)
if ~handles.CRV_add % default
    handles.CRV_range = [];
    handles.CRV_type = [];
    handles.CRV_valarr = [];
    guidata(hObject, handles);
    visualize(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function curve_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to curve_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clear_curves_axes.
function clear_curves_axes_Callback(hObject, eventdata, handles)
    if handles.CRV_add
        n = length(handles.CRV_type)-1;
        if 0==n
            handles.CRV_add = false;
            handles.CRV_range = [];
            handles.CRV_type = [];        
        else
            handles.CRV_range = handles.CRV_range(1:n,:);
            handles.CRV_type = handles.CRV_type(1:n);
        end
        guidata(hObject, handles);
        visualize(hObject,handles);
    end




% --- Executes on button press in show_clustering.
function show_clustering_Callback(hObject, eventdata, handles)

n_clusters = get(handles.n_clusters,'Value') + 1 % :)

        MI_tracks = handles.TrackPlotter_handles.MI_tracks;
        MI_norm_FRET_ratio = handles.TrackPlotter_handles.MI_norm_FRET_ratio;
        MI_norm_nuc_size = handles.TrackPlotter_handles.MI_norm_nuc_size; 
        MI_fnames = handles.TrackPlotter_handles.MI_fnames;
        MI_track_indices = handles.TrackPlotter_handles.MI_track_indices;
        MI_peak_shift = handles.TrackPlotter_handles.MI_peak_shift;  
        dt = handles.TrackPlotter_handles.dt;
        pixelsize = handles.TrackPlotter_handles.pixelsize;

[N,mi_length] = size(MI_norm_FRET_ratio);

FRET_ratio = zeros(N,mi_length);
nucleus_size = zeros(N,mi_length);
Pearson = zeros(N,mi_length);
nneighbours = zeros(N,mi_length);
cell_density = zeros(N,mi_length);

for k=1:size(MI_tracks,1)
    track = MI_tracks{k};
    %
    FRET_ratio(k,:) = track(:,4);
    nucleus_size(k,:) = track(:,7);
    Pearson(k,:) = track(:,8);
    nneighbours(k,:) = track(:,9);
    cell_density(k,:) = track(:,10);
end

avr_FRET = median(FRET_ratio,2);
avr_nucleus_size = median(nucleus_size,2)*pixelsize*pixelsize;
avr_Pearson = median(Pearson,2);
avr_nneighbours = median(nneighbours,2);
avr_cell_density = median(cell_density,2);

% data = MI_norm_FRET_ratio;
% data = [data MI_peak_shift avr_FRET avr_nucleus_size avr_Pearson avr_nneighbours avr_cell_density];
% data = [data avr_FRET avr_Pearson avr_nneighbours avr_cell_density];
data = [MI_norm_FRET_ratio avr_FRET avr_Pearson avr_nneighbours avr_cell_density MI_peak_shift];

% dimensionality reduction
truncation = 2;
V = cov(data);
SD = sqrt(diag(V));
R = V./(SD*SD');
COEFF = pcacov(R);
U = data*COEFF;
data = U(:,1:truncation);                                

mode = get(handles.clustering_method,'Value');
switch mode
    case 1
        % KMEANS
        IDX = kmeans(data,n_clusters);
        %IDX = kmeans(data,n_clusters,'OnlinePhase','on','Distance','correlation');
        %IDX = kmeans(data,n_clusters,'Replicates',10);
        %IDX = kmeans(data,n_clusters,'Start','cluster');
        %IDX = kmeans(data,n_clusters,'OnlinePhase','on','Distance','correlation','Start','cluster','Replicates',10);        
    case 2
        % KMEDOIDS
        IDX = kmedoids(data,n_clusters);
        %IDX = kmedoids(data,n_clusters,'Algorithm','clara');        
    case 3
        % HIERARCHICAL CLUSTERING
        Z = linkage(data,'ward');
        IDX = cluster(Z,'Maxclust',n_clusters);        
end

% PCA COORDINATES
% [coeff,score] = pca(data);
% figure;
% colormap(jet)
% scatter(score(:,1),score(:,2),20,IDX,'filled');
% colorbar(gca)
% grid on;

% CANONICAL COORDINATES
[~,~,stats] = manova1(data,IDX);
C1 = stats.canon(:,1);
C2 = stats.canon(:,2);

    colors = zeros(7,3);
    colors(1,:) = [0 0.4470 0.7410];
    colors(2,:) = [0.8500 0.3250 0.0980];
    colors(3,:) = [0.9290 0.6940 0.1250];
    colors(4,:) = [0.4940 0.1840 0.5560];
    colors(5,:) = [0.4660 0.6740 0.1880];
    colors(6,:) = [0.3010 0.7450 0.9330];
    colors(7,:) = [0.6350 0.0780 0.1840];
    markers = {'o','s','^','d','p'};
    styles = {'-',':','--','-.','-',':','--'};
    
str = get(handles.clustering_method,'String');   
figure('units','normalized','outerposition',[0 0 1 1],'name',[' clustering method: ' str{mode}]);
nrows = 3;
ncols = 5;
h_canon = subplot(nrows,ncols,1);
h_curves = subplot(nrows,ncols,2);
maxy=-inf;
miny=inf;
t = dt*(0:(size(MI_norm_FRET_ratio,2)-1));

LEGEND = [];
for k=1:n_clusters
    subplot(h_curves);
    hold(h_curves,'on');
    m = MI_norm_FRET_ratio(IDX==k,:);
    crv = mean(m,1);
    maxy = max(maxy,max(crv));
    miny = min(miny,min(crv));
    plot(h_curves,t,crv,'marker',markers{k},'linestyle',styles{k},'markersize',4,'linewidth',2);
    axis(h_curves,[min(t) max(t) miny maxy]);
    xlabel(h_curves,'time [min]');
    ylabel(h_curves,'normalized FRET ratio');
    grid(h_curves,'on');
    hold(h_curves,'off');
    
    subplot(h_canon);
    hold(h_canon,'on');
    c1 = C1(IDX==k);
    c2 = C2(IDX==k);
    plot(h_canon,c1,c2,'color',colors(k,:),'marker',markers{k},'linestyle','none','markersize',4,'linewidth',2);
    xticks([]);yticks([]);
    set(h_canon,'XColor','none');
    set(h_canon,'YColor','none');
    hold(h_canon,'off');
    
    h = subplot(nrows,ncols,ncols+k);
    subplot(h);
    hold(h,'on');
    imagesc(m);
    xticks([]);yticks([]);   
    axis(h,[1 size(m,2) 1 size(m,1)]);
    title(h,['cl. ' num2str(k) ', ncells = ' num2str(size(m,1))],'color',colors(k,:));    
    hold(h,'off');
    
    LEGEND = [LEGEND {['cl. ' num2str(k)]}];
    
end
hold off
subplot(h_curves);
legend(h_curves,LEGEND);

PD = zeros(5,N);
PD(1,:) = avr_FRET;
PD(2,:) = avr_nucleus_size;
PD(3,:) = avr_Pearson;
PD(4,:) = avr_nneighbours;
PD(5,:) = avr_cell_density;
PD_labels = {'FRET ratio','nucleus size [\mum^2]','Pearson coeff','#neighbours','cell density [\mum^{-2}]'};
for m=1:5
    h = subplot(nrows,ncols,2*ncols+m);
    subplot(h);
    d = PD(m,:);
    
    x = linspace(min(d),max(d),40);
    for k=1:n_clusters
        hold(h,'on');
        d_k = d(IDX==k);
        histogram(h,d_k,x,'facecolor',colors(k,:),'facealpha',.5,'edgecolor','none','normalization','probability');
        % histogram(h,d_k,x,'facecolor','none','edgecolor',colors(k,:),'displaystyle','stairs','linewidth',2);
    end
    xlabel(h,PD_labels{m});
    %yticklabels(h,[]);
    if 1==m
        ylabel(h,'probability');
    end
    grid(h,'on');
    hold(h,'off');   
end

% --- Executes on selection change in clustering_method.
function clustering_method_Callback(hObject, eventdata, handles)
% hObject    handle to clustering_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns clustering_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from clustering_method


% --- Executes during object creation, after setting all properties.
function clustering_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clustering_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in n_clusters.
function n_clusters_Callback(hObject, eventdata, handles)
% hObject    handle to n_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns n_clusters contents as cell array
%        contents{get(hObject,'Value')} returns selected item from n_clusters


% --- Executes during object creation, after setting all properties.
function n_clusters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
