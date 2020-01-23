function varargout = t_dependent_Nuclei_ratio_FRET_TrackPlotter(varargin)
% T_DEPENDENT_NUCLEI_RATIO_FRET_TRACKPLOTTER MATLAB code for t_dependent_Nuclei_ratio_FRET_TrackPlotter.fig
%      T_DEPENDENT_NUCLEI_RATIO_FRET_TRACKPLOTTER, by itself, creates a new T_DEPENDENT_NUCLEI_RATIO_FRET_TRACKPLOTTER or raises the existing
%      singleton*.
%
%      H = T_DEPENDENT_NUCLEI_RATIO_FRET_TRACKPLOTTER returns the handle to a new T_DEPENDENT_NUCLEI_RATIO_FRET_TRACKPLOTTER or the handle to
%      the existing singleton*.
%
%      T_DEPENDENT_NUCLEI_RATIO_FRET_TRACKPLOTTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in T_DEPENDENT_NUCLEI_RATIO_FRET_TRACKPLOTTER.M with the given input arguments.
%
%      T_DEPENDENT_NUCLEI_RATIO_FRET_TRACKPLOTTER('Property','Value',...) creates a new T_DEPENDENT_NUCLEI_RATIO_FRET_TRACKPLOTTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before t_dependent_Nuclei_ratio_FRET_TrackPlotter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to t_dependent_Nuclei_ratio_FRET_TrackPlotter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help t_dependent_Nuclei_ratio_FRET_TrackPlotter

% Last Modified by GUIDE v2.5 31-Oct-2018 21:19:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @t_dependent_Nuclei_ratio_FRET_TrackPlotter_OpeningFcn, ...
                   'gui_OutputFcn',  @t_dependent_Nuclei_ratio_FRET_TrackPlotter_OutputFcn, ...
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


% --- Executes just before t_dependent_Nuclei_ratio_FRET_TrackPlotter is made visible.
function t_dependent_Nuclei_ratio_FRET_TrackPlotter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to t_dependent_Nuclei_ratio_FRET_TrackPlotter (see VARARGIN)

% Choose default command line output for t_dependent_Nuclei_ratio_FRET_TrackPlotter
handles.output = hObject;

handles.figureName = get(handles.figure1,'Name');

handles.features = {'duration [h]', ...
                    'XY speed [um/min]', ...
                    'XY directionality', ...
                    '#neighbours', ...
                    'cell density', ...
                    'FRET ratio', ...
                    'FRET ratio variability', ...
                    'Donor intensity',...
                    'Acc. intensity',...
                    'nucleus size [um^2]',...
                    'D/A Pearson corr.',... 
                    't(start) [h]',...
                    'FRET ratio c-time [min]',...
                    'FRET molar fraction', ...
                    'Cell area (ref) [um^2]',...
                    'Cell intensity (ref)',...
                    'Nucleus intensity (ref)',...
                    'Nuc/Cyt intensity (ref)',...
                    'Nuc/Cell area (ref)',...                    
                    };

handles.mask = [];    
handles.track_data = [];

handles.velocity_t = [];
handles.nuc_cell_area_ratio_t = [];

% ranges for visualization
handles.rng_duration = [0 0]; % defined in calculate_track_data
handles.rng_start_time = [0 0]; % ditto
handles.rng_speed = [0 10];
handles.rng_directionality = [-1 1];
handles.rng_nghbrs = [1 12];
handles.rng_cell_density = [1e-7 0.008];
handles.rng_FRET_ratio = [0.3 7];
handles.rng_FRET_ratio_variability = [0 0.5];
handles.rng_intensity = [0 1000];
handles.rng_nucleus_size = [20 750];
handles.rng_Pearson_corr = [-1 1];
handles.rng_autocorr_time = [0 120];
handles.rng_FRET_molar_fraction = [0 1];
%
handles.rng_cell_area = [20 750*10];
handles.rng_cell_intensity_ref = [0 1000];
handles.rng_nuc_intensity_ref = [0 1000];
handles.rng_intensity_ref_nuc_cyt_ratio = [0 7];
handles.rng_nuc_cell_area_ratio = [0 1];
% ranges for visualization

%handles.mitotic_intervals_option = 'EXCLUDE';
handles.mitotic_intervals_option = 'LEAVE ONLY';

if 1 == nargin-3
    handles.raw_data = [];
    handles.dt = 1/12; % 5 minutes
    handles.pixelsize = 2; % microns
    handles.track_breaking_flag = varargin{1};
elseif 5 == nargin-3
    tracks = varargin{1};
    handles.dt = varargin{2};
    handles.pixelsize = varargin{3};
    handles.track_breaking_flag = varargin{5};
        % convention - the data saved by ALYtools are not refined, so one needs
        % to refine it now
        handles.raw_data = refine_tracks_by_sorting_mitotic_intervals(handles,tracks);
        %handles.raw_data = tracks;
    [handles.track_data,handles.velocity_t,handles.nuc_cell_area_ratio_t] = calculate_track_data(hObject,handles);
    %
    % set filename in window title
    set(handles.figure1, 'Name', [handles.figureName ' : ' varargin{4}]);
    %
    % single mat file full name
    % CALL SYNTAX - t_dependent_Nuclei_ratio_FRET_TrackPlotter({fname});
elseif 2 == nargin-3 
    if 2==exist(char(varargin{1}))
    handles.track_breaking_flag = varargin{2};
    [filepath,name,ext] = fileparts(char(varargin{1}));
    filename = [name ext];
    %
    load(char(varargin{1}));    
    handles.fullfilename = char(varargin{1});
    %
    if ~exist('microns_per_pixel','var'), return, end
    
    if exist('NUC_STATS','var')
        handles.NUC_STATS = NUC_STATS;
        handles.cell_nums = cell_nums;
    end
                
        handles.dt = dt;
        handles.pixelsize = microns_per_pixel;
            set(handles.pixel_size_edit,'String',handles.pixelsize);
            set(handles.delta_t_edit,'String',handles.dt);    
        %
        % convention - the data saved by ALYtools are not refined, so one needs
        % to refine it now optionally
        if handles.track_breaking_flag    
            handles.raw_data = refine_tracks_by_sorting_mitotic_intervals(handles,tracks);
            if isempty(handles.raw_data)
                disp('refine_tracks_by_sorting_mitotic_intervals - failed');
                handles.raw_data = tracks;                
            end
        else
            handles.raw_data = tracks;
        end
        % this object is for visualizing       
        [handles.track_data,handles.velocity_t,handles.nuc_cell_area_ratio_t] = calculate_track_data(hObject,handles);

        set(handles.figure1, 'Name', [handles.figureName ' : ' filename]);

        str = [{'time'} handles.features];
        minmaxlimits = zeros(numel(str),2);
        minmaxlimits(:,1)=-Inf;
        minmaxlimits(:,2)=Inf;
        if ~isempty(handles.track_data)
        minmaxlimits(1,1)=min(squeeze(handles.track_data(:,1)));
        minmaxlimits(1,2)=max(squeeze(handles.track_data(:,2)));
            for k=1:numel(handles.features)
                minmaxlimits(k+1,1)=min(squeeze(handles.track_data(:,k+2)));
                minmaxlimits(k+1,2)=max(squeeze(handles.track_data(:,k+2)));
            end
        end
        set(handles.filter_table, 'Data', minmaxlimits);

        handles.mask = calculate_mask(hObject,handles);                       
    end
else
    disp('wrong arguments, cant start');
    return;
end

guidata(hObject, handles);

set(handles.time_plot_Y_feature,'String',handles.features);
set(handles.time_plot_colour_feature,'String',handles.features);
set(handles.histo2_Y_feature,'String',handles.features);
set(handles.histo2_X_feature,'String',handles.features);

set(handles.pixel_size_edit,'String',handles.pixelsize);
set(handles.delta_t_edit,'String',handles.dt);
set(handles.histo2_mode,'String',{'scatter','histo2'});

str = [{'time'} handles.features];
minmaxlimits = zeros(numel(str),2);
minmaxlimits(:,1)=-Inf;
minmaxlimits(:,2)=Inf;

if ~isempty(handles.track_data)
minmaxlimits(1,1)=min(squeeze(handles.track_data(1,:)));
minmaxlimits(1,2)=max(squeeze(handles.track_data(2,:)));
       for k=1:numel(handles.features)
            minmaxlimits(k+1,1)=min(squeeze(handles.track_data(:,k+2)));
            minmaxlimits(k+1,2)=max(squeeze(handles.track_data(:,k+2)));
       end
set(handles.filter_table, 'Data', minmaxlimits);       
handles.mask = calculate_mask(hObject,handles);        
end

set(handles.filter_table, 'RowName', str);
set(handles.filter_table, 'ColumnName', {'min','max'});
set(handles.filter_table,'CellEditCallback',@filter_check_callback);

set(handles.time_plot_axes, 'xticklabel', [], 'yticklabel', []);
set(handles.histo2_axes, 'xticklabel', [], 'yticklabel', []);

set(handles.histo2_X_feature,'Value',6);
set(handles.histo2_Y_feature,'Value',7);
set(handles.histo2_mode,'Value',2);

set(handles.time_plot_Y_feature,'Value',6);
set(handles.time_plot_colour_feature,'Value',7);

set(handles.actual_dependence_checkbox,'Value',true);

visualize_histo2(hObject,handles);
visualize_time_dependence(hObject,handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes t_dependent_Nuclei_ratio_FRET_TrackPlotter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line
function varargout = t_dependent_Nuclei_ratio_FRET_TrackPlotter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning out   put args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function pixel_size_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pixel_size_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixel_size_edit as text
%        str2double(get(hObject,'String')) returns contents of pixel_size_edit as a double
if isempty(handles.raw_data), return, end;

value = str2double(get(hObject,'String'));
if ~isnan(value) && value > 0
    handles.pixelsize = value;
    guidata(hObject,handles);
else 
    set(hObject,'String','1');
    guidata(hObject,handles);
    return;
end
[handles.track_data,handles.velocity_t,handles.nuc_cell_area_ratio_t] = calculate_track_data(hObject,handles);
minmaxlimits(1,1)=min(squeeze(handles.track_data(1,:)));
minmaxlimits(1,2)=max(squeeze(handles.track_data(2,:)));
       for k=1:numel(handles.features)
            minmaxlimits(k+1,1)=min(squeeze(handles.track_data(:,k+2)));
            minmaxlimits(k+1,2)=max(squeeze(handles.track_data(:,k+2)));
       end
set(handles.filter_table, 'Data', minmaxlimits);
visualize_histo2(hObject,handles);
visualize_time_dependence(hObject,handles);
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function pixel_size_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixel_size_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function delta_t_edit_Callback(hObject, eventdata, handles)
% hObject    handle to delta_t_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delta_t_edit as text
%        str2double(get(hObject,'String')) returns contents of delta_t_edit as a double
if isempty(handles.raw_data), return, end;

value = str2double(get(hObject,'String'));
if ~isnan(value) && value > 0 
    handles.dt = value;
    guidata(hObject,handles);
else 
    set(hObject,'String','0.1');
    guidata(hObject,handles);
    return;
end
[handles.track_data,handles.velocity_t,handles.nuc_cell_area_ratio_t] = calculate_track_data(hObject,handles);
minmaxlimits(1,1)=min(squeeze(handles.track_data(1,:)));
minmaxlimits(1,2)=max(squeeze(handles.track_data(2,:)));
       for k=1:numel(handles.features)
            minmaxlimits(k+1,1)=min(squeeze(handles.track_data(:,k+2)));
            minmaxlimits(k+1,2)=max(squeeze(handles.track_data(:,k+2)));
       end
set(handles.filter_table, 'Data', minmaxlimits);
visualize_histo2(hObject,handles);
visualize_time_dependence(hObject,handles);
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function delta_t_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delta_t_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in histo2_mode.
function histo2_mode_Callback(hObject, eventdata, handles)
% hObject    handle to histo2_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns histo2_mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from histo2_mode
visualize_histo2(hObject,handles);

% --- Executes during object creation, after setting all properties.
function histo2_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to histo2_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in histo2_Y_feature.
function histo2_Y_feature_Callback(hObject, eventdata, handles)
% hObject    handle to histo2_Y_feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns histo2_Y_feature contents as cell array
%        contents{get(hObject,'Value')} returns selected item from histo2_Y_feature
visualize_histo2(hObject,handles);

% --- Executes during object creation, after setting all properties.
function histo2_Y_feature_CreateFcn(hObject, eventdata, handles)
% hObject    handle to histo2_Y_feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in histo2_X_feature.
function histo2_X_feature_Callback(hObject, eventdata, handles)
% hObject    handle to histo2_X_feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns histo2_X_feature contents as cell array
%        contents{get(hObject,'Value')} returns selected item from histo2_X_feature
visualize_histo2(hObject,handles);

% --- Executes during object creation, after setting all properties.
function histo2_X_feature_CreateFcn(hObject, eventdata, handles)
% hObject    handle to histo2_X_feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in time_plot_Y_feature.
function time_plot_Y_feature_Callback(hObject, eventdata, handles)
% hObject    handle to time_plot_Y_feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns time_plot_Y_feature contents as cell array
%        contents{get(hObject,'Value')} returns selected item from time_plot_Y_feature
visualize_time_dependence(hObject,handles);

% --- Executes during object creation, after setting all properties.
function time_plot_Y_feature_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_plot_Y_feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in time_plot_colour_feature.
function time_plot_colour_feature_Callback(hObject, eventdata, handles)
% hObject    handle to time_plot_colour_feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns time_plot_colour_feature contents as cell array
%        contents{get(hObject,'Value')} returns selected item from time_plot_colour_feature
visualize_time_dependence(hObject,handles);

% --- Executes during object creation, after setting all properties.
function time_plot_colour_feature_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_plot_colour_feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%-------------------------------------------------------------------------%    
function filter_check_callback(hObject,callbackdata)
    handles = guidata(hObject);
    handles.mask = calculate_mask(hObject,handles);
        % Update handles structure
        guidata(hObject, handles);    
    visualize_histo2(hObject,handles);
    visualize_time_dependence(hObject,handles);
    update_possible_visualizer(hObject,handles);

% --------------------------------------------------------------------
function load_trackmate_plus_data_Callback(hObject, eventdata, handles)
% hObject    handle to load_trackmate_plus_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [filename,pathname] = uigetfile({'*.mat','TrackMate+ Files'}, ...
                    'Select data file',pwd);
    if filename == 0, return, end;
    load_trackmate_plus_data(pathname,filename,hObject,handles);
            
% --------------------------------------------------------------------    
function load_trackmate_plus_data(pathname,filename,hObject,handles)

    load([pathname filesep filename]);
    handles.fullfilename = [pathname filesep filename];

    if ~exist('microns_per_pixel','var'), return, end
    
    if exist('NUC_STATS','var')
        handles.NUC_STATS = NUC_STATS;
        handles.cell_nums = cell_nums;
    end
                
    handles.dt = dt;
    handles.pixelsize = microns_per_pixel;
        set(handles.pixel_size_edit,'String',handles.pixelsize);
        set(handles.delta_t_edit,'String',handles.dt);    
    %
    % convention - the data saved by ALYtools are not refined, so one needs
    % to refine it now - optionally
    if handles.track_breaking_flag
        handles.raw_data = refine_tracks_by_sorting_mitotic_intervals(handles,tracks);
        if isempty(handles.raw_data)
            disp('refine_tracks_by_sorting_mitotic_intervals - failed');
            handles.raw_data = tracks;
        end            
    else
        handles.raw_data = tracks;
    end
    % this object is for visualizing       
    [handles.track_data,handles.velocity_t,handles.nuc_cell_area_ratio_t] = calculate_track_data(hObject,handles);
       
    set(handles.figure1, 'Name', [handles.figureName ' : ' filename]);

    str = [{'time'} handles.features];
    minmaxlimits = zeros(numel(str),2);
    minmaxlimits(:,1)=-Inf;
    minmaxlimits(:,2)=Inf;
    if ~isempty(handles.track_data)
    minmaxlimits(1,1)=min(squeeze(handles.track_data(:,1)));
    minmaxlimits(1,2)=max(squeeze(handles.track_data(:,2)));
        for k=1:numel(handles.features)
            minmaxlimits(k+1,1)=min(squeeze(handles.track_data(:,k+2)));
            minmaxlimits(k+1,2)=max(squeeze(handles.track_data(:,k+2)));
        end
    end
    set(handles.filter_table, 'Data', minmaxlimits);
    
    handles.mask = calculate_mask(hObject,handles);

    visualize_histo2(hObject,handles);
    visualize_time_dependence(hObject,handles);
    
    % Update handles structure
    guidata(hObject, handles);

% --------------------------------------------------------------------
function [track_data,velocity_t,nuc_cell_area_ratio_t] = calculate_track_data(hObject,handles)
D = handles.raw_data;
if isempty(D), return, end;

velocity_t = cell(numel(D),1);
nuc_cell_area_ratio_t = cell(numel(D),1);

track_data = zeros(numel(D),2+numel(handles.features)); % first 2 are reserved for begin and end time
for k=1:numel(D)
    track = D{k};
    track_data(k,1) = (track(1,1)-1)*handles.dt;
    track_data(k,2) = track(size(track,1)-1,1)*handles.dt;
    track_data(k,1+2) = size(track,1)*handles.dt;
    FRET_ratio = squeeze(track(:,4));   
    donor_intensity = squeeze(track(:,5));  
    acceptor_intensity = squeeze(track(:,6));
    nucleus_size = squeeze(track(:,7));
    Pearson_correlation = squeeze(track(:,8));
    %
    nnghb = zeros(size(FRET_ratio));
    cell_density = zeros(size(FRET_ratio));
    if size(track,2) >= 10 
        nnghb = squeeze(track(:,9));
        cell_density = squeeze(track(:,10));            
    end
    beta_FRET = 0;            
    if size(track,2)>=11
        nnghb = squeeze(track(:,9));
        cell_density = squeeze(track(:,10));
        beta_FRET = squeeze(track(:,11));        
    end    
    %
    mean_FRET_ratio = mean(FRET_ratio);
    mean_donor_intensity = mean(donor_intensity);
    mean_acceptor_intensity = mean(acceptor_intensity);
    mean_nucleus_size = mean(nucleus_size);
    mean_Pearson_correlation = mean(Pearson_correlation);
    mean_nnghb = mean(nnghb);
    mean_cell_density = mean(cell_density);    
    mean_beta_FRET = mean(beta_FRET);    
    %
    track_data(k,6+2) = mean_FRET_ratio;
    %
    %smoothed_FRET_ratio = medfilt2(FRET_ratio,[5 1]);
    %track_data(k,7+2) = norm(smoothed_FRET_ratio - mean_FRET_ratio);
    %
%     figure(22);
%     if length(FRET_ratio)>100 && length(FRET_ratio)<1000
%         s = FRET_ratio;
%         W = 40;
%         sm_s = medfilt2(FRET_ratio,[5 1]);
%         [~,lf] = TD_high_pass_filter(sm_s,W);        
%         subplot(3,1,1);
%         plot(1:length(s),s,'b.-',1:length(s),lf,'r-',1:length(s),sm_s,'k-');
%         axis([0 length(s)-1 0.5 1.8])
%         legend({'FRET ratio - original','LF','smoothed'});                
%         grid on;
%         
%         subplot(3,1,2);
%         plot(1:length(s),sm_s-lf,'b.-');
%         axis([1 1000 -0.4 0.4]);
%         legend('FRET ratio - filtered');        
%         grid on;
%         
%         subplot(3,1,3);
%         s1 = medfilt2(nucleus_size,[5 1]);
%         s2 = sm_s-lf;
%         s1 = (s1-mean(s1(:)))/(max(s1(:)) - min(s1(:)));
%         s2 = (s2-mean(s2(:)))/(max(s2(:)) - min(s2(:)));
%         plot(1:length(s1),s1,'r.-',1:length(s2),s2,'b.-');
%         axis([1 1000 -1 1]);
%         grid on;
%         legend({'nucleus size','FRET ratio'});
%         disp(k);
%      end

        dt = handles.dt*60; % interval between frames in minutes 
        big_smoothing_window  = round((40*7)/dt);
        small_smoothing_window  = round((5*7)/dt);
        sm_s = medfilt2(FRET_ratio,[small_smoothing_window 1]);
        [~,lf] = TD_high_pass_filter(sm_s,big_smoothing_window);
        s = sm_s-lf;
        track_data(k,7+2) = sqrt(mean(s.*s)); % FRET ratio variability   
        
    %quantify tracks XY trajectory
    x = squeeze(track(:,2));
    y = squeeze(track(:,3));
    z = zeros(size(x));
    [directionality,velocity,velocity_sd] = quantify_track(x',y',z','noZ');
    track_data(k,3+2) = directionality;
    %
    DD = [x';y'];
    velocity_k = zeros(size(x));
    for m=1:numel(x)-1
        velocity_k(m) = norm(DD(1:2,m+1)-DD(1:2,m));
    end
    velocity_k(end) = velocity_k(end-1);
    velocity_t{k} = velocity_k*handles.pixelsize/dt;    
    track_data(k,2+2) = mean(velocity_t{k}); % :0        
    track_data(k,8+2) = mean_donor_intensity;
    track_data(k,9+2) = mean_acceptor_intensity;
    try
    track_data(k,10+2) = mean_nucleus_size*(handles.pixelsize)^2;
    track_data(k,11+2) = mean_Pearson_correlation;
    track_data(k,12+2) = track_data(k,1); % start time
    catch
    end
    %
    track_data(k,4+2) = mean_nnghb;
    track_data(k,5+2) = mean_cell_density/(handles.pixelsize)^2;        
    %
    % autocorr. time
    [ac,lags,bounds] = autocorr(s);
    t=max(bounds);  
        for mm=1:length(lags)
            if ac(mm)<t, break, end
        end
    index = mm-1;
    deix=(ac(index)-t)/(ac(index)-ac(index+1));
    critlag = lags(index)+deix;
    critlag=critlag*dt;
    if critlag < 1, critlag = 1; end
    if critlag > 120, critlag = 120; end
    %
    track_data(k,13+2) = critlag; % autocorr. time
    track_data(k,14+2) = mean_beta_FRET;
    %
%         figure(22);
%         plot(lags,ac/ac(1),'k.-',[0 length(lags)-1],[t/ac(1) t/ac(1)],'r:','linewidth',2)
%         axis([0 length(lags) -1 1]);
%         grid on;
%         legend({'autocorrelation','up.confidence bound'},'fontsize',16);
%         disp(k);
    %
    if 19 == numel(handles.features)
        try
        track_data(k,15+2) = mean(squeeze(track(:,12)))*(handles.pixelsize)^2;
        track_data(k,16+2) = mean(squeeze(track(:,13)));
        track_data(k,17+2) = mean(squeeze(track(:,14)));
        track_data(k,18+2) = mean(squeeze(track(:,15)));
        nuc_cell_area_ratio_t{k} = nucleus_size./squeeze(track(:,12));
        track_data(k,19+2) = mean(nuc_cell_area_ratio_t{k}); % nucleus to cell area ratio
        catch
        end
    end           
end

dur = squeeze(track_data(:,1+2));
stt = squeeze(track_data(:,12+2));
handles.rng_duration = [min(dur) max(dur)];
handles.rng_start_time = [min(stt) max(stt)];

%
% proper place to handle trend curves
%     amin=[];
%     amax=[];
%     for k=1:numel(D)
%         track = D{k};
%         min_frame = track(1,1);
%         max_frame = track(size(track,1),1);        
%         amin = [amin min_frame];
%         amax = [amax max_frame];        
%     end
%     f_min = min(amin);
%     f_max = max(amax);
%     frames = f_min:f_max;
%     %
%                     %1   'duration [h]', ...
%                     %2   'XY speed [um/min]', ...
%                     %3   'XY directionality', ...
%                     %4*   '#neighbours', ...
%                     %5*   'cell density', ...
%                     %6*   'FRET ratio', ...
%                     %7   'FRET ratio variability', ...
%                     %8*   'Donor intensity',...
%                     %9*   'Acc. intensity',...
%                     %10*   'nucleus size [um^2]',...
%                     %11*   'D/A Pearson corr.',... 
%                     %12   't(start) [h]'
%     s_nnghb = cell(1,numel(frames));
%     s_cell_density = cell(1,numel(frames));
%     s_FRET_ratio = cell(1,numel(frames));
%     s_donor_intensity = cell(1,numel(frames));
%     s_acc_intensity = cell(1,numel(frames));
%     s_nucleus_size = cell(1,numel(frames));
%     s_Pearson = cell(1,numel(frames));
%     for k=1:numel(D)
%         track = D{k};
%         for m=1:size(track,1)
%             f = track(m,1);
%                 nnghb = track(m,9);
%                 cell_density = track(m,10);
%                 FRET_ratio = track(m,4);
%                 donor_intensity = track(m,5);
%                 acc_intensity = track(m,6);
%                 nucleus_size = track(m,7);
%                 Pearson = track(m,8);
%                 %
%                 s_nnghb{f} = [s_nnghb{f} nnghb];
%                 s_cell_density{f} = [s_cell_density{f} cell_density];
%                 s_FRET_ratio{f} = [s_FRET_ratio{f} FRET_ratio];
%                 s_donor_intensity{f} = [s_donor_intensity{f} donor_intensity];
%                 s_acc_intensity{f} = [s_acc_intensity{f} acc_intensity];
%                 s_nucleus_size{f} = [s_nucleus_size{f} nucleus_size];
%                 s_Pearson{f} = [s_Pearson{f} Pearson];                
%         end
%     end
% proper place to handle trend curves

% --------------------------------------------------------------------
function visualize_histo2(hObject,handles)

D = handles.track_data;
if isempty(D), return, end;

% need to use filtered data here
x_ind = get(handles.histo2_X_feature,'Value');
y_ind = get(handles.histo2_Y_feature,'Value');
x_data = squeeze(D(:,x_ind+2));
y_data = squeeze(D(:,y_ind+2));

mask = handles.mask; % shouldn't be empty
x_data=x_data(mask==1);
y_data=y_data(mask==1);

str = get(handles.histo2_mode,'String');
mode = str{get(handles.histo2_mode,'Value')};

[x_min_val, x_max_val] = visualization_range(handles,x_ind);
[y_min_val, y_max_val] = visualization_range(handles,y_ind);
                                                         
if strcmp(mode,'scatter')
    plot(handles.histo2_axes,x_data,y_data,'r.');
    if ismember(x_ind,[2:11 13 14]) && ismember(y_ind,[2:11 13 14])
        axis(handles.histo2_axes,[x_min_val x_max_val y_min_val y_max_val]);
    end
    grid(handles.histo2_axes,'on');
elseif strcmp(mode,'histo2')
    corr_map_W = 100;    
    if ismember(x_ind,[2:11 13 14]) && ismember(y_ind,[2:11 13 14])
        x_data = [x_data; x_min_val; x_max_val];
        x_data(x_data<x_min_val)=x_min_val;
        x_data(x_data>x_max_val)=x_max_val;
        %
        y_data = [y_data; y_min_val; y_max_val];
        y_data(y_data<y_min_val)=y_min_val;
        y_data(y_data>y_max_val)=y_max_val;        
    end        
    corrmap = correlation_map(x_data,y_data,corr_map_W);       
    AXES = handles.histo2_axes;
    imagesc(corrmap,'Parent',AXES);
    
    cmap = jet(256);       
    cmap(1,:)=[0,0,0];
    colormap(AXES,cmap);         
        
    daspect(AXES,[1 1 1]);
    set(AXES, 'xticklabel', [], 'yticklabel', []);
end

% --------------------------------------------------------------------
function visualize_time_dependence(hObject,handles)
D = handles.track_data;
if isempty(D), return, end;

DR = handles.raw_data;

% need to use filtered data here
c_ind = get(handles.time_plot_colour_feature,'Value');
y_ind = get(handles.time_plot_Y_feature,'Value');
c_data = squeeze(D(:,c_ind+2));
y_data = squeeze(D(:,y_ind+2));
tb_data = squeeze(D(:,1));
te_data = squeeze(D(:,2));

mask = handles.mask;

n1=numel(y_data);
n2=sum(mask);
set(handles.track_numbers,'String',[num2str(n1) ' : ' num2str(n2)]);

% define colors
Ngrades = 256;
Colors = jet(Ngrades);
min_val = min(c_data);
max_val = max(c_data);
% ??? [min_val, max_val] = visualization_range(handles,c_ind);

show_actual_dependence = false;
%if ismember(y_ind,[4 5 6 8 9 10 11 14 15 16 17 18])
if ismember(y_ind,[4 5 6 8 9 10 11 14 15 16 17 18 2 19])
    show_actual_dependence = true;
end

cla(handles.time_plot_axes,'reset');
axes(handles.time_plot_axes);

for k = 1:numel(y_data)
    if 1==mask(k)
        index = max(1,round((c_data(k) - min_val)/(max_val - min_val)*Ngrades));
        X = [tb_data(k) te_data(k)];
        Y = [y_data(k) y_data(k)];
        if show_actual_dependence && get(handles.actual_dependence_checkbox,'Value') 
            %
            % infer the X,Y data for plotting depending on what it is
            track = DR{k};            
            X = track(:,1)*handles.dt;            
            Y = [];
                    %1   'duration [h]', ...
                    %2   'XY speed [um/min]', ...
                    %3   'XY directionality', ...
                    %4*   '#neighbours', ...
                    %5*   'cell density', ...
                    %6*   'FRET ratio', ...
                    %7   'FRET ratio variability', ...
                    %8*   'Donor intensity',...
                    %9*   'Acc. intensity',...
                    %10*   'nucleus size [um^2]',...
                    %11*   'D/A Pearson corr.',... 
                    %12   't(start) [h]'            
                switch y_ind                    
                    case 4 % #nghbrs                    
                        if ismember(size(track,2),[10 11 15])
                            Y = squeeze(track(:,9));
                        end
                    case 5 % cell density
                        if ismember(size(track,2),[10 11 15])
                            Y = squeeze(track(:,10))/(handles.pixelsize)^2;
                        end
                    case 6 % FRET ratio
                        Y = squeeze(track(:,4));
                    case 8 % donor intensity                     
                        Y = squeeze(track(:,5));
                    case 9 % acceptor intensity
                        Y = squeeze(track(:,6));
                    case 10 % nucleus size
                        Y = squeeze(track(:,7))*(handles.pixelsize)^2;
                    case 11 % Pearson                      
                        Y = squeeze(track(:,8));
                    case 14 % FRET molar fraction
                        try
                        Y = squeeze(track(:,11));
                        catch
                        end
                    case 15 % cell size
                        try
                        Y = squeeze(track(:,12))*(handles.pixelsize)^2;
                        catch
                        end
                    case 16 % 
                        try
                        Y = squeeze(track(:,13));
                        catch
                        end
                    case 17 % 
                        try
                        Y = squeeze(track(:,14));
                        catch
                        end
                    case 18 % 
                        try
                        Y = squeeze(track(:,15));
                        catch
                        end                                                
                    case 2 % speed                        
                        Y = squeeze(handles.velocity_t{k});
                    case 19 % 
                        Y = squeeze(handles.nuc_cell_area_ratio_t{k});                        
                end                                        
            %
        end
        if isempty(Y) % back
            X = [tb_data(k) te_data(k)];            
            Y = [y_data(k) y_data(k)];
        end
        plot(handles.time_plot_axes,X,Y,'Color',Colors(index,:));
        hold(handles.time_plot_axes,'on');
        plot(handles.time_plot_axes,X(1),Y(1),'k.');
        hold(handles.time_plot_axes,'on');        
    end
end

if ~(show_actual_dependence && get(handles.actual_dependence_checkbox,'Value')) 
    plot(handles.time_plot_axes,tb_data(mask==1),y_data(mask==1),'k.');
end

% try to plot stats curve if possible
if isfield(handles,'NUC_STATS')    
              index = 0;
              switch y_ind                    
                    case 4 % #nghbrs
                        index = 6;
                    case 5 % cell density
                        index = 7;
                    case 6 % FRET ratio
                        index = 5;
                    case 8 % donor intensity                     
                        index = 3;
                    case 9 % acceptor intensity
                        index = 4;
                    case 10 % nucleus size
                        index = 1;
                    case 11 % Pearson
                        index = 2;
                    case 14 % FRET molar fraction
                        index = 8;                                                                                                
                    case 15 % 
                        index = 9;  %12; % cell size                                                                       
                    case 16 % 
                        index = 10; %13; % cell intensity ref                                                                       
                    case 17 % 
                        index = 11; %14; % nuke intensity ref                                                                       
                    case 18 % 
                        index = 12; %15; % nuc/cell ref intensity ratio                                                                                                                                               
                    case 19 % 
                        index = 13; %16; % nuc/cell ref area ratio                                                                                                                       
              end
              
    if 0~=index  
        mean_std = get(handles.show_per_frame_mean_std,'Value');
        if mean_std && index <= size(handles.NUC_STATS,2)
            meanvals = squeeze(handles.NUC_STATS(:,index,1)); % mean
            stdvals = squeeze(handles.NUC_STATS(:,index,2)); % std
            taxis = handles.dt*(1:numel(meanvals));
            if ~isempty(intersect(index,[1 9]))
                meanvals = meanvals*(handles.pixelsize)^2;
                stdvals = stdvals*(handles.pixelsize)^2;        
            end
            if 7==index
                meanvals = meanvals/(handles.pixelsize)^2;
                stdvals = stdvals/(handles.pixelsize)^2;      
            end                
            errorbar(taxis,meanvals,stdvals,'Color','black','Marker','o','MarkerFaceColor','magenta','linewidth',1); %    
        else
    %         medvals = squeeze(handles.NUC_STATS(:,index,3)); % median
    %         negvals = squeeze(handles.NUC_STATS(:,index,4)); % 025Q
    %         posvals = squeeze(handles.NUC_STATS(:,index,5)); % 075Q        
    %         taxis = handles.dt*(1:numel(medvals));
    %         if 1==index
    %             medvals = medvals*(handles.pixelsize)^2;
    %             negvals = negvals*(handles.pixelsize)^2;
    %             posvals = posvals*(handles.pixelsize)^2;
    %         end
    %         if 7==index
    %             medvals = medvals/(handles.pixelsize)^2;
    %             negvals = negvals/(handles.pixelsize)^2;
    %             posvals = posvals/(handles.pixelsize)^2;
    %         end                
    %         errorbar(taxis,medvals,negvals,posvals,'Color','black','Marker','o','MarkerFaceColor','magenta','linewidth',1); %            
        end
    end
    hold(handles.time_plot_axes,'on');    
end

hold(handles.time_plot_axes,'off');

[min_val, max_val] = visualization_range(handles,c_ind);
                            
try % calm down if there is no data,     
    c = colorbar(handles.time_plot_axes,'TickLabels',{linspace(min_val,max_val,11)});        
    caxis(handles.time_plot_axes,[min_val,max_val]);
    cmap = jet(256);       
    colormap(handles.time_plot_axes,cmap);                
    c.Label.String = handles.features(c_ind);            
catch
end

try
axis(handles.time_plot_axes,[min(tb_data) max(te_data) min(y_data) max(y_data)]);
catch 
end

[y_min_val, y_max_val] = visualization_range(handles,y_ind);
axis(handles.time_plot_axes,[min(tb_data) max(te_data) y_min_val y_max_val]);                                      

xlabel(handles.time_plot_axes,'time [h]');
%str = get(handles.time_plot_Y_feature,'String') 
%ylabel(handles.time_plot_axes,str{y_ind});
grid(handles.time_plot_axes,'on');


% --------------------------------------------------------------------
function mask = calculate_mask(hObject,handles)
    D = handles.track_data;
    if isempty(D), return, end

    minmaxlimits = get(handles.filter_table, 'Data');
    minmaxlimits = minmaxlimits(2:size(minmaxlimits,1),:);

    D = D(:,3:size(D,2));

    mask = ones(size(D,1),1);
    
    for k=1:numel(handles.features)
        cur_min_val = minmaxlimits(k,1);
        cur_max_val = minmaxlimits(k,2);
        cur_vals = squeeze(D(:,k));
        cur_vals(isnan(cur_vals)) = cur_min_val; % safety
        cur_mask = cur_vals>=cur_min_val & cur_vals<=cur_max_val;
        mask = mask & cur_mask;
    end
    
function update_possible_visualizer(hObject,handles)  
if ~get(handles.update_vidi_immediately,'Value'), return, end;

FIGS = findobj(0, 'type', 'figure');
for k=1:size(FIGS,1)
    h = FIGS(k,1);
    cur_handles = guidata(h);
    if ~strcmp(cur_handles.figure1.Name,handles.figure1.Name)
        try
        if strcmp(cur_handles.TrackPlotter_handles.fullfilename,handles.fullfilename) % mine
            %@(hObject,eventdata)t_dependent_Nuclei_ratio_FRET_visualizer('update_button_Callback',hObject,eventdata,guidata(hObject));  
            t_dependent_Nuclei_ratio_FRET_visualizer('update_button_Callback',cur_handles.update_button,[],cur_handles);
        end
        catch
        end
    end
end
     
% --------------------------------------------------------------------
function save_settings_as_matfile_Callback(hObject, eventdata, handles)
% hObject    handle to save_settings_as_matfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [fname, fpath] = uiputfile('*.mat','Save Settings as..',[pwd filesep 'TrackPlotter_settings']);
                if fpath == 0; return; end
                filespec = fullfile(fpath,fname);

    delta_t = handles.dt;
    pixelsize = handles.pixelsize; %microns
    filter_table_data = get(handles.filter_table, 'Data');

    try
        save(filespec,'delta_t','pixelsize','filter_table_data');
    catch
        errordlg('Error while trying to save settings - not saved');
    end

% --------------------------------------------------------------------
function load_settings_from_matfile_Callback(hObject, eventdata, handles)
% hObject    handle to load_settings_from_matfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [filename,pathname] = uigetfile({'*.mat','TrackPlotter settings'}, ...
                    'Select the settings file',pwd);
    if filename == 0, return, end;
    load([pathname filesep filename]);

    try
        handles.dt = delta_t;
        handles.pixelsize = pixelsize; %microns
        set(handles.filter_table, 'Data', filter_table_data);
    catch
        errordlg('Error while trying to load settings - not set up');
    end


% --- Executes on button press in actual_dependence_checkbox.
function actual_dependence_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to actual_dependence_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of actual_dependence_checkbox
visualize_time_dependence(hObject,handles);

% this is highly specific function that exploits correlation between
% FRET_ratio(t) and nuclear_size(t) during mitosis, relying on many hardcoded params
function ref_tracks = refine_tracks_by_sorting_mitotic_intervals(handles,tracks)
ref_tracks = {};
dt = handles.dt*60; % interval between frames in minutes 
for k=1:numel(tracks)
    track = tracks{k};
    FRET_ratio = squeeze(track(:,4));
    nucleus_size = squeeze(track(:,7));
    %
    big_smoothing_window  = round((40*7)/dt);
    small_smoothing_window  = round((5*7)/dt);
    correlation_window  = round((15*7)/dt);
    t = 0.08; % threshold
    fill_little_gaps_size = round((3*7)/dt);
    pre_mit = round((3*7)/dt);
    post_mit = round((5*7)/dt);
    min_track_length = round((10*7)/dt);
    
    % interested in mitotic intervals
    if ~strcmpi('exclude',handles.mitotic_intervals_option)
        pre_mit = pre_mit*2;
        post_mit = post_mit*2;        
    end
    %
    sm_s = medfilt2(FRET_ratio,[small_smoothing_window 1]); % small smoothing window - HARDCODED
    [~,lf] = TD_high_pass_filter(sm_s,big_smoothing_window);
    s1 = sm_s-lf;
    % repeat for nuclei size
    sm_s = medfilt2(nucleus_size,[small_smoothing_window 1]);
    [~,lf] = TD_high_pass_filter(sm_s,big_smoothing_window);
    s2 = sm_s-lf;
    %
    % some black magic ...
    %
    % normalize signals
    s1 = (s1-mean(s1(:)))/(max(s1(:)) - min(s1(:)));
    s2 = (s2-mean(s2(:)))/(max(s2(:)) - min(s2(:)));    
    r = movcorr(s1,s2,correlation_window);    
    r = r.*s1.*s2.*(s1>0 & s2>0); % where both local signals are positive
    r = r > t; % thresholding
    %
    r = imclose(r,ones(fill_little_gaps_size,1)); % to fill little gaps - HARDCODED
    
%      if length(s1)>200
%         figure(22);
%         f=1:length(s1);
%         h1=gca;
%         plot(f,s1,'k.-',f,s2,'b.-',f,r,'g-');   
%         grid(h1,'on');
%         pause(1e-3); % put your breakpoint here
%     end   
    
    
%      if length(s1)>200
%         figure(22);
%         f=1:length(s1);
%         h1 = subplot(2,1,1);
%         plot(f,s1,'k.-');   
%         grid(h1,'on');
%         legend(h1,{'trend-subtracted FRET ratio'},'fontsize',14);
%         axis(h1,[0 199 -0.8 0.8]);
%         h3=subplot(2,1,2);
%         [ac,lags,bounds] = autocorr(s1);
%         plot(lags,ac,'k.-',lags,max(bounds)*ones(size(lags)),'r:','linewidth',2);
%         t=max(bounds);  
%          for mm=1:length(lags)
%             if ac(mm)<t, break, end
%         end
%         index = mm-1;
%         deix=(ac(index)-t)/(ac(index)-ac(index+1));
%         critlag = lags(index)+deix;
%         title(h3,[ num2str(critlag) ' [frame]']);
%         grid(h3,'on');
%         legend(h3,{'autocorrelation','up. confidence bound'},'fontsize',14);
%         pause(1e-3); % put your breakpoint here
%     end

    excl = zeros(size(r)); % exclusion mask   
    L = bwlabel(r);
    s = regionprops(L,'Centroid');
    for m=1:numel(s)
        c = s(m).Centroid;
        x = round(c(2));        
        min_ind = max(1,x-pre_mit);
        max_ind = min(length(r),x+post_mit);
        excl(min_ind:max_ind) = 1;
    end    
    %
%      if length(s1)>200
%         figure(22);
%         f=1:length(s1);
%         h1=gca;
%         plot(f,s1,'k.-',f,s2,'b.-',f,excl,'g-');   
%         grid(h1,'on');
%         pause(1e-3); % put your breakpoint here
%     end       
    
    if strcmpi('exclude',handles.mitotic_intervals_option)
        continue;
    else
        excl = ~excl;
    end
    %
    indices = 1:length(excl);
    if 0~=sum(excl(:))
        L = bwlabel(~excl);
        for z=1:max(L)
            s = indices'.*(L==z);
            s=s(s~=0);
            part_track = track(min(s):max(s),:);
            if size(part_track,1) > min_track_length % minimal length of a track - HARDCODED
                ref_tracks = [ref_tracks; part_track];
            end
        end
    end          
end


% --------------------------------------------------------------------
function Tools_Callback(hObject, eventdata, handles)
% hObject    handle to Tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function show_cell_numbers_Callback(hObject, eventdata, handles)
% hObject    handle to show_cell_numbers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'cell_nums'), return, end
t = handles.dt*(0:numel(handles.cell_nums)-1);
h=figure;
plot(t,handles.cell_nums,'k.-');
axis([t(1) t(numel(handles.cell_nums)) min(handles.cell_nums) max(handles.cell_nums)]);
    xlabel('time [h]');
    ylabel('#cells');
grid on;    
%
figurename = get(handles.figure1,'Name');
str = strsplit(figurename,(' : '));
figurename = char(str(2));
set(h,'Name',figurename);

% --------------------------------------------------------------------
function average_pixel_brightness_Callback(hObject, eventdata, handles)
% hObject    handle to average_pixel_brightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'NUC_STATS'), return, end
don  = squeeze(handles.NUC_STATS(:,3,1));
acc  = squeeze(handles.NUC_STATS(:,4,1));
t= handles.dt*(0:numel(acc)-1);
h = figure;
if 19 ~= numel(handles.features)
        plot(t,don,'b.-',t,acc,'r.-');
        xlabel('time [h]');
        ylabel('intensity');
        legend({'donor','acceptor'});
else
        ref = squeeze(handles.NUC_STATS(:,10,1));
        plot(t,don,'b.-',t,acc,'r.-',t,ref,'k.-');
        xlabel('time [h]');
        ylabel('intensity');
        legend({'donor','acceptor','ref'});        
end               
grid on;
ax_new=gca;
set(ax_new,'Position','default');
legend(ax_new,'-DynamicLegend');
%
figurename = get(handles.figure1,'Name');
str = strsplit(figurename,(' : '));
figurename = char(str(2));
set(h,'Name',figurename);

% --------------------------------------------------------------------
% --- Executes on button press in show_per_frame_mean_std.
function show_per_frame_mean_std_Callback(hObject, eventdata, handles)
% hObject    handle to show_per_frame_mean_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_per_frame_mean_std
visualize_time_dependence(hObject,handles);


% --------------------------------------------------------------------
function visualize_selection_Callback(hObject, eventdata, handles)
% hObject    handle to visualize_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
t_dependent_Nuclei_ratio_FRET_visualizer(handles);


% --- Executes on button press in update_vidi_immediately.
function update_vidi_immediately_Callback(hObject, eventdata, handles)
% hObject    handle to update_vidi_immediately (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of update_vidi_immediately


function [minval, maxval] = visualization_range(handles,index)
             switch index
                     case 2 % #speed
                        minval = min(handles.rng_speed);
                        maxval = max(handles.rng_speed);
                    case 3 % #directionality
                        minval = min(handles.rng_directionality);
                        maxval = max(handles.rng_directionality);                  
                    case 4 % #nghbrs
                        minval = min(handles.rng_nghbrs);
                        maxval = max(handles.rng_nghbrs);  
                    case 5 % cell density
                        minval = min(handles.rng_cell_density);
                        maxval = max(handles.rng_cell_density);  
                    case 6 % FRET ratio
                        minval = min(handles.rng_FRET_ratio);
                        maxval = max(handles.rng_FRET_ratio);  
                    case 7 % FRET ratio variability
                        minval = min(handles.rng_FRET_ratio_variability);
                        maxval = max(handles.rng_FRET_ratio_variability);  
                    case {8,9} % donor intensity                     
                        minval = min(handles.rng_intensity);
                        maxval = max(handles.rng_intensity);  
                    case 10 % nucleus size
                        minval = min(handles.rng_nucleus_size);
                        maxval = max(handles.rng_nucleus_size);  
                    case 11 % Pearson
                        minval = min(handles.rng_Pearson_corr);
                        maxval = max(handles.rng_Pearson_corr);  
                    case 13 % autocorr. time
                        minval = min(handles.rng_autocorr_time);
                        maxval = max(handles.rng_autocorr_time);  
                    case 14 % FRET molar fraction
                        minval = min(handles.rng_FRET_molar_fraction);
                        maxval = max(handles.rng_FRET_molar_fraction);
                     case 1 % duration
                        minval = min(handles.rng_duration);
                        maxval = max(handles.rng_duration);                    
                    case 12 % start time
                        minval = min(handles.rng_start_time);
                        maxval = max(handles.rng_start_time); 
                    case 15 % 
                        minval = min(handles.rng_cell_area);
                        maxval = max(handles.rng_cell_area); 
                    case 16 % 
                        minval = min(handles.rng_cell_intensity_ref);
                        maxval = max(handles.rng_cell_intensity_ref); 
                    case 17 % 
                        minval = min(handles.rng_nuc_intensity_ref);
                        maxval = max(handles.rng_nuc_intensity_ref); 
                    case 18 % 
                        minval = min(handles.rng_intensity_ref_nuc_cyt_ratio);
                        maxval = max(handles.rng_intensity_ref_nuc_cyt_ratio); 
                    case 19 % 
                        minval = min(handles.rng_nuc_cell_area_ratio);
                        maxval = max(handles.rng_nuc_cell_area_ratio);
             end


% --------------------------------------------------------------------
function mitotic_interval_FRET_ratio_heatmap_Callback(hObject, eventdata, handles)
% hObject    handle to mitotic_interval_FRET_ratio_heatmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.mitotic_intervals_option,'LEAVE ONLY')
    t_dependent_Nuclei_ratio_FRET_ratio_heatmapper(handles);
end

