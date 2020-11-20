function varargout = t_dependent_Nuclei_ratio_FRET_HTS_MI_heatmapper(varargin)
% T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_MI_HEATMAPPER MATLAB code for t_dependent_Nuclei_ratio_FRET_HTS_MI_heatmapper.fig
%      T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_MI_HEATMAPPER, by itself, creates a new T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_MI_HEATMAPPER or raises the existing
%      singleton*.
%
%      H = T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_MI_HEATMAPPER returns the handle to a new T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_MI_HEATMAPPER or the handle to
%      the existing singleton*.
%
%      T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_MI_HEATMAPPER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_MI_HEATMAPPER.M with the given input arguments.
%
%      T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_MI_HEATMAPPER('Property','Value',...) creates a new T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_MI_HEATMAPPER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before t_dependent_Nuclei_ratio_FRET_HTS_MI_heatmapper_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to t_dependent_Nuclei_ratio_FRET_HTS_MI_heatmapper_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help t_dependent_Nuclei_ratio_FRET_HTS_MI_heatmapper

% Last Modified by GUIDE v2.5 06-Nov-2020 16:56:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @t_dependent_Nuclei_ratio_FRET_HTS_MI_heatmapper_OpeningFcn, ...
                   'gui_OutputFcn',  @t_dependent_Nuclei_ratio_FRET_HTS_MI_heatmapper_OutputFcn, ...
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


% --- Executes just before t_dependent_Nuclei_ratio_FRET_HTS_MI_heatmapper is made visible.
function t_dependent_Nuclei_ratio_FRET_HTS_MI_heatmapper_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to t_dependent_Nuclei_ratio_FRET_HTS_MI_heatmapper (see VARARGIN)

% Choose default command line output for t_dependent_Nuclei_ratio_FRET_HTS_MI_heatmapper
handles.output = hObject;

    handles.TrackPlotter_handles = varargin{1};
    try
    [~,NAME,EXT]=fileparts(handles.TrackPlotter_handles.fullfilename);
    set(handles.figure1, 'Name', ['FRET ratio heatmap : ' NAME EXT]);
    
    % masking doesn't worrk.. likely because haven't switched to MI view..
                        
    catch
    end

% should start calculating all parameters starting from this index
handles.DF = round(handles.TrackPlotter_handles.MI_LLeft/(handles.TrackPlotter_handles.dt*60));     
handles.DF = handles.DF + 2; % ?????    

%handles.max_frame = floor(handles.TrackPlotter_handles.filter_table.Data(1,2)/handles.TrackPlotter_handles.dt);        
            %
            tracks = handles.TrackPlotter_handles.MI_tracks;
            raw_FRET_vals = [];
            frame_vals = [];            
            for k=1:numel(tracks)
                track = tracks{k};
                    frame_vals = [frame_vals; track(:,1)];
                    raw_FRET_vals = [raw_FRET_vals; track(:,4)];
            end
handles.max_frame = max(frame_vals);
handles.minmax_raw_FRET_vals = minmax(raw_FRET_vals');
z = handles.TrackPlotter_handles.MI_norm_FRET_ratio(:,handles.DF:size(handles.TrackPlotter_handles.MI_norm_FRET_ratio,2));
handles.minmax_norm_FRET_vals = minmax(z(:)');
handles.max_MI_number = numel(tracks);
                        
handles.letters = {'A','B','C','D','E','F','G','H'};
handles.numbers = {'1','2','3','4','5','6','7','8','9','10','11','12'};
    
set(handles.letter_start,'String',handles.letters);
set(handles.letter_end,'String',handles.letters);
set(handles.number_start,'String',handles.numbers);
set(handles.number_end,'String',handles.numbers);
set(handles.number_start,'Value',1);
set(handles.number_end,'Value',12);
set(handles.letter_start,'Value',1);
set(handles.letter_end,'Value',8);
set(handles.highlight_selection,'Value',1);

set(handles.clustering_method,'String',{'kmeans','kmedoids','hierarchical'});

set(handles.plate_map_display,'String',{ ...
    'MI number', ...        
    '<MI time>', ...        
    '|(type 1| / |type1 U type2)|', ... 
    '<FRET ratio>', ...    
    '<FRET ratio after-peak amplitude>', ...         
});
set(handles.plot_type,'String',{ ...
    'mean norm FRET ratio profiles', ...    
    'raw FRET ratio histograms', ...        
    'mean raw FRET ratio profiles', ...
    'MI number (t)', ...
    'sorted raw FRET ratio heatmaps', ...
    });

set(handles.show_curves,'String',{'all','total only','type1,type2 only'});
set(handles.show_curves,'Value',1);

set(handles.platemap_visualization_mode,'String',{'platemap','chart'});
set(handles.platemap_visualization_mode,'Value',1);

set(handles.statistic,'String',{'p-value: KS','p-value: t-test','p-value: Wilcoxon','Cohen"s d','|median diff|'});
set(handles.statistic,'visible','off');

% discrimination chart
handles.chart_panel = uipanel(handles.figure1,'visible','off');
handles.chart_panel.Position = [0.08 0.02 .3*1.15 .5*1.15];     

% platemap 
handles.platemap_panel = uipanel(handles.figure1,'visible','off');
handles.platemap_panel.Position = [0.02 0.2 .45 .4];

    raw_filenames = handles.TrackPlotter_handles.MI_fnames;
    %
    all_tokens_10_12 = [];
    all_tokens = [];    
    for k=1:8
        for m=1:12
            if m>=10
                all_tokens_10_12 = [all_tokens_10_12 {[handles.letters{k} '-' handles.numbers{m}]} ];
            end
            all_tokens = [all_tokens {[handles.letters{k} '-' handles.numbers{m}]} ];
        end
    end
    raw_data_tokens = num2cell(zeros(size(raw_filenames)));
    % first pass
    for k=1:numel(raw_data_tokens)
        for m=1:numel(all_tokens)
        if contains(raw_filenames{k},all_tokens{m})
            raw_data_tokens(k) = all_tokens(m);
        end
        end
    end
    % second pass
    for k=1:numel(raw_data_tokens)
        for m=1:numel(all_tokens_10_12)
        if contains(raw_filenames{k},all_tokens_10_12{m})
            raw_data_tokens(k) = all_tokens_10_12(m);
        end
        end
    end
            
    handles.raw_data_tokens = raw_data_tokens;

handles.h = [];
handles.h_chart = [];
% D = rand(8,12);
% update_platemap(hObject,handles,D);

handles.color1 = [0.8500 0.3250 0.0980];
handles.color2 = [0 0.4470 0.7410];
set(handles.clustering_method,'Value',1);
handles.IDX = do_clustering(handles,false);

handles.vis_panel = uipanel(handles.figure1,'visible','off');
handles.vis_panel.Position = [0.5 0.01 .49 .92];        

% Update handles structure
guidata(hObject, handles);

set(handles.plate_map_display,'Value',4);
visualize_platemap(hObject,handles);
set(handles.plot_type,'Value',4);
update_visuals_axis(hObject,handles);

% UIWAIT makes t_dependent_Nuclei_ratio_FRET_HTS_MI_heatmapper wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = t_dependent_Nuclei_ratio_FRET_HTS_MI_heatmapper_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in letter_start.
function letter_start_Callback(hObject, eventdata, handles)
end_letter = get(handles.letter_end,'Value');
my_letter = get(hObject,'Value');
if my_letter>end_letter
    set(hObject,'Value',end_letter);
end
guidata(hObject, handles);
visualize_platemap(hObject,handles);
update_visuals_axis(hObject,handles);


% --- Executes during object creation, after setting all properties.
function letter_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to letter_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in letter_end.
function letter_end_Callback(hObject, eventdata, handles)
start_letter = get(handles.letter_start,'Value');
my_letter = get(hObject,'Value');
if my_letter<start_letter
    set(hObject,'Value',start_letter);
end
guidata(hObject, handles);
visualize_platemap(hObject,handles);
update_visuals_axis(hObject,handles);


% --- Executes during object creation, after setting all properties.
function letter_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to letter_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in number_start.
function number_start_Callback(hObject, eventdata, handles)
end_number = get(handles.number_end,'Value');
my_number = get(hObject,'Value');
if my_number>end_number
    set(hObject,'Value',end_number);
end
guidata(hObject, handles);
visualize_platemap(hObject,handles);
update_visuals_axis(hObject,handles);

% --- Executes during object creation, after setting all properties.
function number_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in number_end.
function number_end_Callback(hObject, eventdata, handles)
start_number = get(handles.number_start,'Value');
my_number = get(hObject,'Value');
if my_number<start_number
    set(hObject,'Value',start_number);
end
guidata(hObject, handles);
visualize_platemap(hObject,handles);
update_visuals_axis(hObject,handles);


% --- Executes during object creation, after setting all properties.
function number_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function IDX = do_clustering(handles,verbose)
    n_clusters = 2; %get(handles.n_clusters,'Value') + 1; % :)

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
    
    MI_norm_FRET_ratio = MI_norm_FRET_ratio(:,handles.DF+1:size(MI_norm_FRET_ratio,2));
    
    % data = MI_norm_FRET_ratio;
    % data = [data MI_peak_shift avr_FRET avr_nucleus_size avr_Pearson avr_nneighbours avr_cell_density];
    % data = [data avr_FRET avr_Pearson avr_nneighbours avr_cell_density];
    data = [MI_norm_FRET_ratio avr_FRET avr_Pearson avr_nneighbours avr_cell_density];

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

    % order type1 and type2 properly
         m1 = mean(MI_norm_FRET_ratio(IDX==1,:),1);
         m2 = mean(MI_norm_FRET_ratio(IDX==2,:),1);
         if max(m1)<max(m2)
             IDX(IDX==1)=nan;
             IDX(IDX==2)=1;
             IDX(isnan(IDX))=2;                                       
         end         
     
    if verbose

            colors = zeros(7,3);
            colors(2,:) = [0 0.4470 0.7410];
            colors(1,:) = [0.8500 0.3250 0.0980];
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
        t = 60*dt*(0:(size(MI_norm_FRET_ratio,2)-1));

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
            title(h_canon,str{mode}); 
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
    end


% --- Executes on selection change in clustering_method.
function clustering_method_Callback(hObject, eventdata, handles)
% hObject    handle to clustering_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.IDX = do_clustering(handles,false);
guidata(hObject,handles);
visualize_platemap(hObject,handles);
update_visuals_axis(hObject,handles);

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


% --- Executes on selection change in plate_map_display.
function plate_map_display_Callback(hObject, eventdata, handles)
    guidata(hObject, handles);
    visualize_platemap(hObject,handles);

% --- Executes during object creation, after setting all properties.
function plate_map_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plate_map_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in open_as_separate_figure.
function open_as_separate_figure_Callback(hObject, eventdata, handles)
    % turns out to be .. too cumbersome to implement
    plot_ind = get(handles.plot_type,'Value');
    s = get(handles.plot_type,'String');
    plot_name = s{plot_ind};
       
    if strcmp('on',get(handles.vis_panel,'visible'))
        % ??
    elseif strcmp('on',get(handles.visuals_axis,'visible'))
        axes(handles.visuals_axis);
        src_ax = gca;        
        h = figure('Name',plot_name);
        dst_ax = gca;
        copyobj(src_ax.Children,dst_ax);
        grid(dst_ax,'on');          
    end

% --- Executes on selection change in plot_type.
function plot_type_Callback(hObject, eventdata, handles)
    update_visuals_axis(hObject,handles);

% --- Executes during object creation, after setting all properties.
function plot_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in highlight_selection.
function highlight_selection_Callback(hObject, eventdata, handles)
    visualize_platemap(hObject,handles);
    update_visuals_axis(hObject,handles);

function visualize_platemap(hObject,handles)   
[~,indices] = select_wells(handles);
    D = create_platemap_values(handles);
    if get(handles.highlight_selection,'Value')        
        D_ = NaN(8,12);
        for k=1:numel(indices)
            D_(indices{k}(1),indices{k}(2)) = D(indices{k}(1),indices{k}(2));
        end
        update_platemap(hObject,handles,D_);
    else
        update_platemap(hObject,handles,D);
    end
    
function [selected_wells, selected_indices] = select_wells(handles)
selected_wells = cell(0);
selected_indices = cell(0);
for k=get(handles.letter_start,'Value'):get(handles.letter_end,'Value')
    for m=get(handles.number_start,'Value'):get(handles.number_end,'Value')
        selected_wells = [selected_wells [handles.letters{k} '-' handles.numbers{m}] ];
        selected_indices = [selected_indices [k m]];
    end
end

function [letter_ind,number_ind] = get_letter_number_indices(handles,well_token)
    s = strsplit(well_token,'-');
    letter_ind = find(ismember(handles.letters,s{1}));
    number_ind = find(ismember(handles.numbers,s{2}));

function D = create_platemap_values(handles)
    D = NaN(8,12);
    %D = rand(8,12);    

    tracks = handles.TrackPlotter_handles.MI_tracks;
    wells = handles.raw_data_tokens;
    dt = handles.TrackPlotter_handles.dt;
        
    C = cell(8,12);
    
    switch get(handles.plate_map_display,'Value')        
        
        case 1 % total number of Mitotic Intervals
          D = zeros(8,12);
          for w=1:numel(wells)
              [k,m] = get_letter_number_indices(handles,wells{w});
              D(k,m) = D(k,m)+1;
          end            
          
        case 2 % average time where mitosis happens
          for w=1:numel(wells)
             [k,m] = get_letter_number_indices(handles,wells{w});
             mitosis_time = (tracks{w}(1,1)+handles.DF)*dt;
             C{k,m} = [C{k,m}; mitosis_time];
          end
          D = nan(8,12);
          for k=1:8
              for m=1:12
                  if~isempty(C{k,m})
                    D(k,m) = mean(C{k,m});
                  end
              end
          end
            
        case 3 % percentage in type1
             D1 = zeros(8,12);
             D2 = zeros(8,12);
             for w=1:numel(wells)
                [k,m] = get_letter_number_indices(handles,wells{w});
                          if 1==handles.IDX(w)
                              D1(k,m)=D1(k,m)+1;
                          else
                              D2(k,m)=D2(k,m)+1;
                          end                          
             end
             D = D1./(D1+D2);
                        
        case 4 % average FRET ratio
         for w=1:numel(wells)
             [k,m] = get_letter_number_indices(handles,wells{w});
             track = tracks{w};
             FRET_values = track(handles.DF:size(track,1),4);
             C{k,m} = [C{k,m}; FRET_values];
          end
          D = nan(8,12);
          for k=1:8
              for m=1:12
                  if~isempty(C{k,m})
                    D(k,m) = mean(C{k,m});
                  end
              end
          end            
                        
        case 5 % average FRET ratio amplitude
         for w=1:numel(wells)
             [k,m] = get_letter_number_indices(handles,wells{w});
             track = tracks{w};
             FRET_values = track(handles.DF+2:size(track,1),4); % to ensure it is well down
             C{k,m} = [C{k,m}; max(FRET_values)-min(FRET_values)];
          end
          D = nan(8,12);
          for k=1:8
              for m=1:12
                  if~isempty(C{k,m})
                    D(k,m) = mean(C{k,m});
                  end
              end
          end                                                     
    end

function update_platemap(hObject,handles,D)

    bckg_color = get(handles.figure1,'Color');

    switch get(handles.platemap_visualization_mode,'Value')
        
        case 1
            set(handles.statistic,'visible','off');
            set(handles.platemap_panel,'visible','on');
            set(handles.chart_panel,'visible','off');
            
            handles.h = heatmap(handles.platemap_panel,handles.numbers,handles.letters,D);
            handles.h.Colormap = jet;

        %     s = get(handles.parameter,'String');
        %     param_name = s{get(handles.parameter,'Value')};
        %     s = get(handles.statistic,'String');
        %     statistic_name = s{get(handles.statistic,'Value')};
        %     
        %     handles.h.Title = ['t = ' get(handles.evaluation_time,'String'), ...    
        %                        'h , ' param_name,' , ',statistic_name];

            handles.h.MissingDataColor = bckg_color;    
            handles.h.GridVisible = 'off';
            
                switch get(handles.plate_map_display,'Value')
                    case 1 % ?
                    case 3
                    caxis(gca,[0 1]);                    
                end            
            
        case 2

            s = get(handles.plate_map_display,'String');
            param = get(handles.plate_map_display,'Value');
            if ismember(param,[1 3]), return, end            
            
            param_name = s{get(handles.plate_map_display,'Value')};            
            
            set(handles.statistic,'visible','on');
            set(handles.platemap_panel,'visible','off');
            set(handles.chart_panel,'visible','on');
                
            [selected_wells,D] = create_discrimination_chart_data(handles);
            if 96*96 ~= sum(isnan(D),'all')
                handles.h_chart = heatmap(handles.chart_panel,selected_wells,selected_wells,D);
                handles.h_chart.Colormap = jet;

                s = get(handles.statistic,'String');
                statistic_name = s{get(handles.statistic,'Value')};        
                handles.h_chart.Title = [param_name,' , ',statistic_name];

                handles.h_chart.FontSize = 8;
                handles.h_chart.MissingDataColor = bckg_color;    
                handles.h_chart.GridVisible = 'off';
                
                if ismember(get(handles.statistic,'Value'),[1 2 3])
                    caxis(gca,[0 1]);
                end
            end
    end
    
    guidata(hObject, handles);

function update_visuals_axis(hObject,handles)

    tracks = handles.TrackPlotter_handles.MI_tracks;
    wells = handles.raw_data_tokens;
    
% handles.max_frame
% handles.minmax_raw_FRET_vals
% handles.minmax_norm_FRET_vals
% handles.max_MI_number

    dt = handles.TrackPlotter_handles.dt;

    [selected_wells,~] = select_wells(handles);
    %indices = find(ismember(wells,selected_wells));

    plot_ind = get(handles.plot_type,'Value');
    s = get(handles.plot_type,'String');
    plot_name = s{plot_ind};
    
    tfaxis = 1:handles.max_frame;
    
    MI_norm_FRET_ratio = handles.TrackPlotter_handles.MI_norm_FRET_ratio;        
    MI_norm_FRET_ratio = MI_norm_FRET_ratio(:,handles.DF:size(MI_norm_FRET_ratio,2));    
     
    set(handles.vis_panel,'visible','off');
    set(handles.visuals_axis,'visible','on');
    
    switch plot_ind        
               
        %%%%%%%%%%
        case 1 % mean FRET ratio profile
                profiles = [];
                profiles1 = [];
                profiles2 = [];
                for k=1:numel(tracks)
                    if ismember(wells{k},selected_wells)
                        cur_profile = MI_norm_FRET_ratio(k,:)';
                        profiles = [profiles; cur_profile'];                        
                          if 1==handles.IDX(k)
                                profiles1 = [profiles1; cur_profile'];
                            else
                                profiles2 = [profiles2; cur_profile'];
                          end                                                
                    end
                end                
                profile = mean(profiles,1);
                taxis = (0:(length(profile)-1))*dt;
                
                N = size(profiles,1);
                N1 = size(profiles1,1);
                N2 = size(profiles2,1);
                L1 = {['total ' num2str(N)],['type1 ' num2str(N1)],['type2 ' num2str(N2)]};
                L2 = {['total ' num2str(N)]};
                L3 = {['type1 ' num2str(N1)],['type2 ' num2str(N2)]};
                                
                % plot(handles.visuals_axis,taxis,profile,'ro-','linewidth',2);
                cla(handles.visuals_axis,'reset');
                axes(handles.visuals_axis);
                
                switch get(handles.show_curves,'Value')
                    case 1 % all
                        shadedErrorBar(taxis,mean(profiles),std(profiles),'lineProps',{'k.-','linewidth',2,'markersize',24});
                        hold(handles.visuals_axis,'on');
                        shadedErrorBar(taxis,mean(profiles1),std(profiles1),'lineProps',{'linewidth',2,'color',handles.color1,'linestyle',':'},'patchSaturation',0.1);
                        hold(handles.visuals_axis,'on');
                        shadedErrorBar(taxis,mean(profiles2),std(profiles2),'lineProps',{'linewidth',2,'color',handles.color2,'linestyle',':'},'patchSaturation',0.1);
                        hold(handles.visuals_axis,'off');   
                        legend(handles.visuals_axis,L1);
                    case 2 % total only
                        shadedErrorBar(taxis,mean(profiles),std(profiles),'lineProps',{'k.-','linewidth',2,'markersize',24});                        
                        legend(handles.visuals_axis,L2);
                    case 3 % types only                        
                        shadedErrorBar(taxis,mean(profiles1),std(profiles1),'lineProps',{'linewidth',2,'color',handles.color1,'linestyle',':'},'patchSaturation',0.1);
                        hold(handles.visuals_axis,'on');
                        shadedErrorBar(taxis,mean(profiles2),std(profiles2),'lineProps',{'linewidth',2,'color',handles.color2,'linestyle',':'},'patchSaturation',0.1);
                        hold(handles.visuals_axis,'off');                                        
                        legend(handles.visuals_axis,L3);
                end
                
                grid(handles.visuals_axis,'on');
                xlabel(handles.visuals_axis,'time [h]');
                ylabel(handles.visuals_axis,plot_name);
                axis(handles.visuals_axis,[ min(taxis) ...
                                            max(taxis) ...
                                            handles.minmax_norm_FRET_vals(1)...
                                            handles.minmax_norm_FRET_vals(2)]);
        
        case 2 % histograms
            %%%%%%%%%%%%%%%
            avr_FRET = [];
            avr_FRET_1 = [];
            avr_FRET_2 = [];
            for k=1:size(tracks,1)                
                 if ismember(wells{k},selected_wells)
                          track = tracks{k};
                          v = median(track(:,4));
                          avr_FRET = [avr_FRET; v];
                          if 1==handles.IDX(k)
                            avr_FRET_1 = [avr_FRET_1; v];
                          else
                            avr_FRET_2 = [avr_FRET_2; v];
                          end                                                
                 end   
            end  
            
                N = numel(avr_FRET);
                N1 = numel(avr_FRET_1);
                N2 = numel(avr_FRET_2);
                L1 = {['total ' num2str(N)],['type1 ' num2str(N1)],['type2 ' num2str(N2)]};
                L2 = {['total ' num2str(N)]};
                L3 = {['type1 ' num2str(N1)],['type2 ' num2str(N2)]};            
            
            x = linspace(handles.minmax_raw_FRET_vals(1),handles.minmax_raw_FRET_vals(2),40);
            h = handles.visuals_axis;
            cla(handles.visuals_axis,'reset');
                switch get(handles.show_curves,'Value')
                    case 1 % all
                        histogram(h,avr_FRET,x,'facecolor','k','facealpha',.2,'edgecolor','none','normalization','probability');
                        hold(h,'on');
                        histogram(h,avr_FRET_1,x,'facecolor',handles.color1,'facealpha',.2,'edgecolor','none','normalization','probability');
                        hold(h,'on');                       
                        histogram(h,avr_FRET_2,x,'facecolor',handles.color2,'facealpha',.2,'edgecolor','none','normalization','probability');                        
                        hold(h,'on');
                        legend(h,L1);
                    case 2 % total
                        histogram(h,avr_FRET,x,'facecolor','k','facealpha',.2,'edgecolor','none','normalization','probability');
                        legend(h,L2);
                    case 3 % types only
                        histogram(h,avr_FRET_1,x,'facecolor',handles.color1,'facealpha',.2,'edgecolor','none','normalization','probability');
                        hold(h,'on');
                        histogram(h,avr_FRET_2,x,'facecolor',handles.color2,'facealpha',.2,'edgecolor','none','normalization','probability');
                        legend(h,L3);
                end
                hold(h,'off');                   
                ylabel(h,'probability');
                xlabel(handles.visuals_axis,'FRET ratio');
                grid(h,'on');
                                                     
        case 3 % mean FRET ratio profile
                profiles = [];
                profiles1 = [];
                profiles2 = [];
                for k=1:numel(tracks)
                    if ismember(wells{k},selected_wells)
                        cur_profile = tracks{k}(handles.DF:size(tracks{k},1),4);
                        %cur_profile = tracks{k}(:,4);
                        profiles = [profiles; cur_profile'];                        
                          if 1==handles.IDX(k)
                                profiles1 = [profiles1; cur_profile'];
                            else
                                profiles2 = [profiles2; cur_profile'];
                          end                                                
                    end
                end                
                profile = mean(profiles,1);
                taxis = (0:(length(profile)-1))*dt;
                
                N = size(profiles,1);
                N1 = size(profiles1,1);
                N2 = size(profiles2,1);
                L1 = {['total ' num2str(N)],['type1 ' num2str(N1)],['type2 ' num2str(N2)]};
                L2 = {['total ' num2str(N)]};
                L3 = {['type1 ' num2str(N1)],['type2 ' num2str(N2)]};
                                
                % plot(handles.visuals_axis,taxis,profile,'ro-','linewidth',2);
                cla(handles.visuals_axis,'reset');
                axes(handles.visuals_axis);
                
                switch get(handles.show_curves,'Value')
                    case 1 % all
                        shadedErrorBar(taxis,mean(profiles),std(profiles),'lineProps',{'k.-','linewidth',2,'markersize',24});
                        hold(handles.visuals_axis,'on');
                        shadedErrorBar(taxis,mean(profiles1),std(profiles1),'lineProps',{'linewidth',2,'color',handles.color1,'linestyle',':'},'patchSaturation',0.1);
                        hold(handles.visuals_axis,'on');
                        shadedErrorBar(taxis,mean(profiles2),std(profiles2),'lineProps',{'linewidth',2,'color',handles.color2,'linestyle',':'},'patchSaturation',0.1);
                        hold(handles.visuals_axis,'off');   
                        legend(handles.visuals_axis,L1);
                    case 2 % total only
                        shadedErrorBar(taxis,mean(profiles),std(profiles),'lineProps',{'k.-','linewidth',2,'markersize',24});                        
                        legend(handles.visuals_axis,L2);
                    case 3 % types only                        
                        shadedErrorBar(taxis,mean(profiles1),std(profiles1),'lineProps',{'linewidth',2,'color',handles.color1,'linestyle',':'},'patchSaturation',0.1);
                        hold(handles.visuals_axis,'on');
                        shadedErrorBar(taxis,mean(profiles2),std(profiles2),'lineProps',{'linewidth',2,'color',handles.color2,'linestyle',':'},'patchSaturation',0.1);
                        hold(handles.visuals_axis,'off');                                        
                        legend(handles.visuals_axis,L3);
                end
                
                grid(handles.visuals_axis,'on');
                xlabel(handles.visuals_axis,'time [h]');
                ylabel(handles.visuals_axis,plot_name);
                axis(handles.visuals_axis,[ min(taxis) ...
                                            max(taxis) ...
                                            handles.minmax_raw_FRET_vals(1)...
                                            handles.minmax_raw_FRET_vals(2)]);
                                                        
        case 4 % number of intervals (t)                
            numbers = zeros(size(tfaxis));
            numbers1 = zeros(size(tfaxis));
            numbers2 = zeros(size(tfaxis));
            for k=1:numel(tracks)
                if ismember(wells{k},selected_wells)
                    cur_frame = tracks{k}(1,1);
                    numbers(cur_frame) = numbers(cur_frame)+1;
                    if 1==handles.IDX(k)
                        numbers1(cur_frame) = numbers1(cur_frame)+1;
                    else
                        numbers2(cur_frame) = numbers2(cur_frame)+1;                        
                    end
                end
            end
            numbers = cumsum(numbers);
            numbers1 = cumsum(numbers1);
            numbers2 = cumsum(numbers2);
            N = numbers(numel(numbers));
            N1 = numbers1(numel(numbers1));
            N2 = numbers2(numel(numbers2));
            L1 = {['total ' num2str(N)],['type1 ' num2str(N1)],['type2 ' num2str(N2)]};
            L2 = {['total ' num2str(N)]};
            L3 = {['type1 ' num2str(N1)],['type2 ' num2str(N2)]};
            
            switch get(handles.show_curves,'Value')
                case 1 % all            
                    plot(handles.visuals_axis,(tfaxis+handles.DF-1)*dt,numbers,'k:','linewidth',3);
                    hold(handles.visuals_axis,'on');
                    plot(handles.visuals_axis,(tfaxis+handles.DF-1)*dt,numbers1,'linestyle',':','color',handles.color1,'linewidth',3);
                    hold(handles.visuals_axis,'on');
                    plot(handles.visuals_axis,(tfaxis+handles.DF-1)*dt,numbers2,'linestyle',':','color',handles.color2,'linewidth',3);
                    hold(handles.visuals_axis,'off');
                    legend(handles.visuals_axis,L1);
                case 2 % total only
                    plot(handles.visuals_axis,(tfaxis+handles.DF-1)*dt,numbers,'k:','linewidth',3);
                    legend(handles.visuals_axis,L2);
                case 3 % types only           
                    plot(handles.visuals_axis,(tfaxis+handles.DF-1)*dt,numbers1,'linestyle',':','color',handles.color1,'linewidth',3);
                    hold(handles.visuals_axis,'on');
                    plot(handles.visuals_axis,(tfaxis+handles.DF-1)*dt,numbers2,'linestyle',':','color',handles.color2,'linewidth',3);
                    hold(handles.visuals_axis,'off');
                    legend(handles.visuals_axis,L3);
            end            
            %
            grid(handles.visuals_axis,'on');
            xlabel(handles.visuals_axis,'time [h]');
            ylabel(handles.visuals_axis,plot_name);            
                fmax_real = min(find(numbers==max(numbers)));
            axis(handles.visuals_axis,[ 0 (fmax_real+handles.DF-1)*dt 0 handles.max_MI_number ]);
            
        case 5 % heatmaps
                profiles = [];
                profiles1 = [];
                profiles2 = [];
                    avr_FRET = [];
                    avr_FRET_1 = [];
                    avr_FRET_2 = [];
                
                for k=1:numel(tracks)
                    if ismember(wells{k},selected_wells)
                        cur_profile = tracks{k}(handles.DF:size(tracks{k},1),4);
                        %cur_profile = tracks{k}(:,4);
                        v = median(tracks{k}(:,4));
                        avr_FRET = [avr_FRET; v];
                        profiles = [profiles; cur_profile'];                        
                          if 1==handles.IDX(k)
                                avr_FRET_1 = [avr_FRET_1; v];
                                profiles1 = [profiles1; cur_profile'];
                          else
                                avr_FRET_2 = [avr_FRET_2; v];
                                profiles2 = [profiles2; cur_profile'];
                          end                                                
                    end
                end

                [~,I] = sort(avr_FRET,'descend');
                profiles = profiles(I,:);
                    [~,I] = sort(avr_FRET_1,'descend');
                    profiles1 = profiles1(I,:);
                        [~,I] = sort(avr_FRET_2,'descend');
                        profiles2 = profiles2(I,:);
                
                N = size(profiles,1);
                N1 = size(profiles1,1);
                N2 = size(profiles2,1);
                Ltot = {['total ' num2str(N)]};
                L1 = {['type1 ' num2str(N1)]};
                L2 = {['type2 ' num2str(N2)]};

                taxis = round((0:size(profiles,2)-1)*handles.TrackPlotter_handles.dt*60);
                b1=min(taxis);
                b2=max(taxis);
                step = fix(length(taxis)/4);
                step_val = round((b2 - b1)/4);
                                
                set(handles.vis_panel,'visible','on');
                set(handles.visuals_axis,'visible','off');
                
                switch get(handles.show_curves,'Value')
                    case 1 % all
                        ax1 = subplot(1,3,1,'Parent',handles.vis_panel);
                        imagesc(ax1,profiles);
                        ax2 = subplot(1,3,2,'Parent',handles.vis_panel);
                        imagesc(ax2,profiles1);
                        ax3 = subplot(1,3,3,'Parent',handles.vis_panel);
                        imagesc(ax3,profiles2);       
                        xticks(ax1,[]);yticks(ax1,[]);
                        xticks(ax2,[]);yticks(ax2,[]);
                        xticks(ax3,[]);yticks(ax3,[]);
                        title(ax1,Ltot);
                        title(ax2,L1,'color',handles.color1);
                        title(ax3,L2,'color',handles.color2);
                            set_minutes_time_axis_for_imagesc(ax1,taxis,b1,step,step_val);
                            set_minutes_time_axis_for_imagesc(ax2,taxis,b1,step,step_val);
                            set_minutes_time_axis_for_imagesc(ax3,taxis,b1,step,step_val);
                            caxis(ax1,[handles.minmax_raw_FRET_vals(1),handles.minmax_raw_FRET_vals(2)]);
                            caxis(ax2,[handles.minmax_raw_FRET_vals(1),handles.minmax_raw_FRET_vals(2)]);
                            caxis(ax3,[handles.minmax_raw_FRET_vals(1),handles.minmax_raw_FRET_vals(2)]);
                    case 2 % tot 
                        ax1 = subplot(1,1,1,'Parent',handles.vis_panel);
                        imagesc(ax1,profiles);        
                        yticks(ax1,[]);
                        title(ax1,Ltot);                        
                        %
                        set_minutes_time_axis_for_imagesc(ax1,taxis,b1,step,step_val);
                        caxis(ax1,[handles.minmax_raw_FRET_vals(1),handles.minmax_raw_FRET_vals(2)]);                        
                        colorbar(ax1,'EastOutside');
                    case 3 % types       
                        ax2 = subplot(1,2,1,'Parent',handles.vis_panel);
                        imagesc(ax2,profiles1);
                        ax3 = subplot(1,2,2,'Parent',handles.vis_panel);
                        imagesc(ax3,profiles2);
                        xticks(ax2,[]);yticks(ax2,[]);
                        xticks(ax3,[]);yticks(ax3,[]);
                        title(ax2,L1,'color',handles.color1);
                        title(ax3,L2,'color',handles.color2);
                        caxis(ax2,[handles.minmax_raw_FRET_vals(1),handles.minmax_raw_FRET_vals(2)]);
                        caxis(ax3,[handles.minmax_raw_FRET_vals(1),handles.minmax_raw_FRET_vals(2)]);
                        set_minutes_time_axis_for_imagesc(ax2,taxis,b1,step,step_val);
                        set_minutes_time_axis_for_imagesc(ax3,taxis,b1,step,step_val);

                end
                                                                            
    end
    
function set_minutes_time_axis_for_imagesc(ax1,taxis,b1,step,step_val)
                            xticks(ax1,0:step:(length(taxis)));
                            xticklabels(ax1,b1+step_val*(0:4));                            
                            xlabel(ax1,'time [min]');
                            grid(ax1,'on');         

                               
% --- Executes on selection change in show_curves.
function show_curves_Callback(hObject, eventdata, handles)
%
update_visuals_axis(hObject,handles);

% --- Executes during object creation, after setting all properties.
function show_curves_CreateFcn(hObject, eventdata, handles)
% hObject    handle to show_curves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clustering_details.
function clustering_details_Callback(hObject, ~, handles)
    handles.IDX = do_clustering(handles,true);
    guidata(hObject, handles);
    visualize_platemap(hObject,handles);
    update_visuals_axis(hObject,handles);

    
function [selected_wells,D] = create_discrimination_chart_data(handles)
 
    [selected_wells,~] = select_wells(handles);
    N = numel(selected_wells);
    
    sample = cell(N,1); % statistical samples
    
    tracks = handles.TrackPlotter_handles.MI_tracks; 
        
    dt = handles.TrackPlotter_handles.dt;
    
    for k = 1:numel(tracks)
        track_k = tracks{k};  
        index = find(ismember(selected_wells,handles.raw_data_tokens{k}));
        if isempty(index), continue, end
        v = nan;
        switch get(handles.plate_map_display,'Value')
            case 2
                v = (track_k(1,1)+handles.DF-1)*dt; % mitosis_time 
            case 4
                FRET_values = track_k(handles.DF:size(track_k,1),4);
                v = median(FRET_values);
            case 5                
                FRET_values = track_k(handles.DF+2:size(track_k,1),4); % to ensure it is well down
                v = max(FRET_values)-min(FRET_values);                
        end
        if ~isnan(v)
            sample{index} = [sample{index} v];
        end
    end

    D = nan(N);
    
    hw = waitbar(0,'calculating statistics..','WindowStyle','modal');
    for k=1:numel(sample)
        for m=1:numel(sample)
                        x1 = sample{k};
                        x2 = sample{m};
                        
            mode = get(handles.statistic,'Value');
                        
            if ~isempty(x1) && ~isempty(x2) && k<m
                switch mode
                    case 1
                        [~,P] = kstest2(x1,x2);
                    case 2
                        [~,P] = ttest2(x1,x2,'vartype','unequal'); % not sure                        
                    case 3  
                        [P,~] = ranksum(x1,x2);
                    case 4 % Cohen's d
                        N1 = numel(x1);
                        N2 = numel(x2);
                        s = sqrt( 1/(N1+N2)*( (N1-1)*var(x1) + (N2-1)*var(x2) ) );
                        P = abs( mean(x1) - mean(x2) )/s;
                    case 5
                        P = abs(median(x1)-median(x2));
                end
                D(k,m) = P; 
            end            
        end
        if ~isempty(hw), waitbar(k/numel(sample),hw); drawnow, end   
    end
    if ~isempty(hw), delete(hw), drawnow; end
    
    
    
    
    
    
    
    
    
    


% --- Executes on selection change in statistic.
function statistic_Callback(hObject, eventdata, handles)
visualize_platemap(hObject,handles);

% --- Executes during object creation, after setting all properties.
function statistic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to statistic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in platemap_visualization_mode.
function platemap_visualization_mode_Callback(hObject, eventdata, handles)
visualize_platemap(hObject,handles);

% --- Executes during object creation, after setting all properties.
function platemap_visualization_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to platemap_visualization_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
