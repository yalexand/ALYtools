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

% Last Modified by GUIDE v2.5 02-Mar-2021 14:11:21

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
                    'Nucleus solidity',...
                    'Nucleus shape factor',...
                    'Nucleus eccentricity',...
                    'Nucleus orientation',...
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
handles.rng_intensity_ref_nuc_cyt_ratio = [0 20];
handles.rng_nuc_cell_area_ratio = [0 1];
%
handles.rng_nuc_solidity = [0 1];
handles.rng_nuc_shape_factor = [0.95 2];
handles.rng_nuc_eccentricity = [0 1];
handles.rng_nuc_orientation = [-90 90];

% ranges for visualization

handles.timestamp = datestr(now,'yyyy-mm-dd HH-MM-SS'); %to ID from outside 

set(handles.histo2_mode,'String',{'scatter','histo2'});
set(handles.tracks_to_show,'String',{'tracks','mitotic intervals'});

handles.MI_LLeft = 120;
handles.MI_LRight = 140;
     
handles.MED_t_nucsize_drv_uncond     = 0.4;     % threshold for unconditional acceptance (derivative strength)
handles.MED_t_nucsize_drv_cond       = 0.25;    % threshold for conditional acceptance (derivative strength)
handles.MED_dt_twin_max              = 76;      % [min] beyond this value, can consider as next mitotic      
handles.MED_T_min                    = 5;       % don't consider signals shorter than T_min     
handles.MED_big_smoothing_window     = 280;     % [min]
handles.MED_small_smoothing_window   = 35;      % [min]
handles.MED_t_nucsize_ratio          = 1.75;    % dimensionless
handles.MED_DT                       = 10;      % [min] time interval either to the past or future to calculate average nuc_size there
handles.ME                           = [];      % mitotic events storage

guidata(hObject, handles);

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
%         % convention - the data saved by ALYtools are not refined, so one needs
%         % to refine it now
%         handles.raw_data = refine_tracks_by_sorting_mitotic_intervals(handles,tracks);
%         %handles.raw_data = tracks;
    [handles.track_data,handles.velocity_t,handles.nuc_cell_area_ratio_t] = calculate_track_data(handles);
    %
    % set filename in window title
    set(handles.figure1, 'Name', [handles.figureName ' : ' varargin{4}]);
    %
    % single mat file full name
    % CALL SYNTAX - t_dependent_Nuclei_ratio_FRET_TrackPlotter({fname});
elseif 2 == nargin-3 
    if 2==exist(char(varargin{1}))
        handles.track_breaking_flag = varargin{2};
        [pathname,name,ext] = fileparts(char(varargin{1}));
        filename = [name ext];    

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
            set(handles.pixel_size_edit,'String',num2str(handles.pixelsize));
            set(handles.delta_t_edit,'String',num2str(handles.dt));    
            
        handles.raw_data = tracks;

        % this object is for visualizing       
        [handles.track_data,handles.velocity_t,handles.nuc_cell_area_ratio_t] = calculate_track_data(handles);

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

set(handles.time_plot_Y_feature,'String',handles.features);
set(handles.time_plot_colour_feature,'String',handles.features);
set(handles.histo2_Y_feature,'String',handles.features);
set(handles.histo2_X_feature,'String',handles.features);

set(handles.pixel_size_edit,'String',num2str(handles.pixelsize));
set(handles.delta_t_edit,'String',num2str(handles.dt));

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

%visualize_histo2(hObject,handles);
%visualize_time_dependence(hObject,handles);

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
[handles.track_data,handles.velocity_t,handles.nuc_cell_area_ratio_t] = calculate_track_data(handles);
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
[handles.track_data,handles.velocity_t,handles.nuc_cell_area_ratio_t] = calculate_track_data(handles);
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
    if filename == 0, return, end
    load_trackmate_plus_data(pathname,filename,hObject,handles);
            
% --------------------------------------------------------------------    
function load_trackmate_plus_data(pathname,filename,hObject,handles)

    single_FOVs_options(handles,'on');

    if isfield(handles,'filenames') && ~isempty(handles.filenames)
        handles.filenames = [];
    end
    
    load([pathname filesep filename]);
    handles.fullfilename = [pathname filesep filename];

    if ~exist('microns_per_pixel','var'), return, end
    
    if exist('NUC_STATS','var')
        handles.NUC_STATS = NUC_STATS;
        handles.cell_nums = cell_nums;
    end
                
    handles.dt = dt;
    handles.pixelsize = microns_per_pixel;
        set(handles.pixel_size_edit,'String',num2str(handles.pixelsize));
        set(handles.delta_t_edit,'String',num2str(handles.dt));
    %   
    handles.ST_raw_data = tracks;
    [handles.features_lut,handles.features_void,handles.features_coeff] = ... 
                                    set_features_to_data_correspondence(handles);
    guidata(hObject,handles);
  
    handles.MI_tracks = cell(0);
    handles.MI_fnames = cell(0);
    handles.MI_track_indices = [];
    handles.MI_norm_FRET_ratio = [];
    handles.MI_norm_nuc_size = [];
    handles.MI_peak_shift = [];
    
    handles.ST_raw_data_filenames = repmat({filename},[size(tracks,1) 1]); % corresponding filenames
        
    handles.ME = detect_mitotic_events(handles,'on'); % verbose  
    %
    [MI_tracks, ... 
    MI_track_indices, ...
    MI_norm_FRET_ratio, ...
    MI_norm_nuc_size, ...
    MI_peak_shift] = get_mitotic_intervals_by_using_mitotic_events(handles,'on');
    %
    if ~isempty(MI_tracks)
        handles.MI_tracks = MI_tracks;
        handles.MI_norm_FRET_ratio = MI_norm_FRET_ratio;
        handles.MI_norm_nuc_size = MI_norm_nuc_size; 
        handles.MI_fnames = repmat(filename,[size(MI_tracks,1) 1]);
        handles.MI_track_indices = MI_track_indices;
        handles.MI_peak_shift = MI_peak_shift;
    end
                                        
    tracks_to_show_Callback(hObject, [], handles);

% --------------------------------------------------------------------
function [track_data,velocity_t,nuc_cell_area_ratio_t] = calculate_track_data(handles)
D = handles.raw_data;
if isempty(D), return, end

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
    critlag = 0;
    try
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
    catch
        disp('autocorr time - failed to calculate');
    end
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
       
% size of params may be or 15 or 16 or 20
    switch size(track,2)
        case 15 % 3-rd channel only
            track_data(k,15+2) = mean(squeeze(track(:,12)))*(handles.pixelsize)^2;
            track_data(k,16+2) = mean(squeeze(track(:,13)));
            track_data(k,17+2) = mean(squeeze(track(:,14)));
            track_data(k,18+2) = mean(squeeze(track(:,15)));
            nuc_cell_area_ratio_t{k} = nucleus_size./squeeze(track(:,12));
            track_data(k,19+2) = mean(nuc_cell_area_ratio_t{k}); % nucleus to cell area ratio
        case 16 % % no 3rd channel, with nuc shapes
            track_data(k,20+2) = mean(squeeze(track(:,12))); % solidity
            track_data(k,21+2) = mean(squeeze(track(:,13))); % shape factor
            track_data(k,22+2) = mean(squeeze(track(:,14))); % eccentricity
            track_data(k,23+2) = mean(squeeze(track(:,15))); % orientation        
        case 20 % with 3rd channel, with nuc shapes
            track_data(k,15+2) = mean(squeeze(track(:,12)))*(handles.pixelsize)^2;
            track_data(k,16+2) = mean(squeeze(track(:,13)));
            track_data(k,17+2) = mean(squeeze(track(:,14)));
            track_data(k,18+2) = mean(squeeze(track(:,15)));
            nuc_cell_area_ratio_t{k} = nucleus_size./squeeze(track(:,12));
            track_data(k,19+2) = mean(nuc_cell_area_ratio_t{k}); % nucleus to cell area ratio
            track_data(k,20+2) = mean(squeeze(track(:,16))); % solidity
            track_data(k,21+2) = mean(squeeze(track(:,17))); % shape factor
            track_data(k,22+2) = mean(squeeze(track(:,18))); % eccentricity
            track_data(k,23+2) = mean(squeeze(track(:,19))); % orientation                
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

cla(handles.histo2_axes,'reset');
set(handles.histo2_axes, 'xticklabel', [], 'yticklabel', []);
axes(handles.histo2_axes);

D = handles.track_data;
if isempty(D), return, end

% need to use filtered data here
x_ind = get(handles.histo2_X_feature,'Value');
y_ind = get(handles.histo2_Y_feature,'Value');
x_data = squeeze(D(:,x_ind+2));
y_data = squeeze(D(:,y_ind+2));

if isempty(x_data) || 0==sum(x_data(:)) || isempty(y_data) || 0==sum(y_data(:))
    return,
end

mask = handles.mask; % shouldn't be empty
x_data=x_data(mask==1);
y_data=y_data(mask==1);

str = get(handles.histo2_mode,'String');
mode = str{get(handles.histo2_mode,'Value')};

[x_min_val, x_max_val] = visualization_range(handles,x_ind);
[y_min_val, y_max_val] = visualization_range(handles,y_ind);
                                                         
if strcmp(mode,'scatter')
    plot(handles.histo2_axes,x_data,y_data,'r.');
    if ismember(x_ind,[2:11 13:23]) && ismember(y_ind,[2:11 13:23])
        axis(handles.histo2_axes,[x_min_val x_max_val y_min_val y_max_val]);
    end
    grid(handles.histo2_axes,'on');
elseif strcmp(mode,'histo2')
    corr_map_W = 100;    
    if ismember(x_ind,[2:11 13:23]) && ismember(y_ind,[2:11 13:23])
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

% handles.features = {'duration [h]', ...                   1
%                     'XY speed [um/min]', ...              2
%                     'XY directionality', ...              3
%                     '#neighbours', ...                    4
%                     'cell density', ...                   5
%                     'FRET ratio', ...                     6
%                     'FRET ratio variability', ...         7
%                     'Donor intensity',...                 8
%                     'Acc. intensity',...                  9
%                     'nucleus size [um^2]',...             10
%                     'D/A Pearson corr.',...               11
%                     't(start) [h]',...                    12
%                     'FRET ratio c-time [min]',...         13
%                     'FRET molar fraction', ...            14
%                     'Cell area (ref) [um^2]',...          15
%                     'Cell intensity (ref)',...            16
%                     'Nucleus intensity (ref)',...         17
%                     'Nuc/Cyt intensity (ref)',...         18
%                     'Nuc/Cell area (ref)',...             19
%                     'Nucleus solidity',...                20
%                     'Nucleus shape factor',...            21
%                     'Nucleus eccentricity',...            22
%                     'Nucleus orientation',...             23
%                     };
% --------------------------------------------------------------------
function visualize_time_dependence(hObject,handles)

cla(handles.time_plot_axes,'reset');
axes(handles.time_plot_axes);

c_ind = get(handles.time_plot_colour_feature,'Value');
y_ind = get(handles.time_plot_Y_feature,'Value');

% some may be not available
if ismember(y_ind,handles.features_void) || ismember(c_ind,handles.features_void), return, end

D = handles.track_data;
if isempty(D), return, end

c_data = squeeze(D(:,c_ind+2));
y_data = squeeze(D(:,y_ind+2));
tb_data = squeeze(D(:,1));
te_data = squeeze(D(:,2));

% need to use filtered data here
mask = handles.mask;

n1=numel(y_data);
n2=sum(mask);
set(handles.track_numbers,'String',[num2str(n1) ' : ' num2str(n2)]);

% define colors
Ngrades = 256;
Colors = jet(Ngrades);
c_min_val = min(c_data);
c_max_val = max(c_data);
%[c_min_val, c_max_val] = visualization_range(handles,c_ind);

can_show_actual_dependence = false;
%
if ismember(y_ind,[2 4:6 8:11 14:23])
    can_show_actual_dependence = true;
end

for k = 1:numel(y_data)
    if 1==mask(k)
        index = max(1,round((c_data(k) - c_min_val)/(c_max_val - c_min_val)*Ngrades));
        if can_show_actual_dependence && get(handles.actual_dependence_checkbox,'Value') 
            %
            track = handles.raw_data{k};
            X = track(:,1)*handles.dt;
            %
            if 2==y_ind
               Y = squeeze(handles.velocity_t{k}); 
            elseif 19==y_ind
               Y = squeeze(handles.nuc_cell_area_ratio_t{k});
            else
               % infer the X,Y data for plotting depending on what it is
               Y = squeeze(track(:,handles.features_lut(y_ind)))*handles.features_coeff(y_ind);
            end
        else
            X = [tb_data(k) te_data(k)];
            Y = [y_data(k) y_data(k)];            
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

if ~(can_show_actual_dependence && get(handles.actual_dependence_checkbox,'Value')) 
    plot(handles.time_plot_axes,tb_data(mask==1),y_data(mask==1),'k.');
end

% try to plot stats curve if possible
if isfield(handles,'NUC_STATS') && ~ismember(y_ind,[7 12 13]) && strcmp('on',get(handles.show_per_frame_mean_std,'Enable'))   
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
              end              

            if      15 == size(track,2) % +3-rd channel, no nuc.shapes, 12 size of NUC_STATS
                switch y_ind
                    case 15 %
                        index = 9;  %12; % cell size                                                                       
                    case 16 % 
                        index = 10; %13; % cell intensity ref                                                                       
                    case 17 % 
                        index = 11; %14; % nuke intensity ref                                                                       
                    case 18 % 
                        index = 12; %15; % nuc/cell ref intensity ratio
                end
            elseif  16 == size(track,2) % no 3-rd channel +nuc.shapes, 13 size of NUC_STATS
                switch y_ind
                    case 20 %
                        index = 9;  %12 solidity                               
                    case 21 % 
                        index = 10; %13 shape factor
                    case 22 % 
                        index = 11; %14 eccentricity
                    case 23 % 
                        index = 12; %15 orientation
                end                
            elseif  20 == size(track,2) % +3-rd channel, +nuc.shapes 17 size of NUC_STATS
                switch y_ind
                    case 15 %
                        index = 9;  %12; % cell size                                                                       
                    case 16 % 
                        index = 10; %13; % cell intensity ref                                                                       
                    case 17 % 
                        index = 11; %14; % nuke intensity ref                                                                       
                    case 18 % 
                        index = 12; %15; % nuc/cell ref intensity ratio
                    case 20 %
                        index = 13;  %12 solidity                               
                    case 21 % 
                        index = 14; %13 shape factor
                    case 22 % 
                        index = 15; %14 eccentricity
                    case 23 % 
                        index = 16; %15 orientation
                end                
            end
                                          
    if 0~=index  
        mean_std = get(handles.show_per_frame_mean_std,'Value');
        if mean_std && index <= size(handles.NUC_STATS,2)
            meanvals = squeeze(handles.NUC_STATS(:,index,1)); % mean
            stdvals = squeeze(handles.NUC_STATS(:,index,2)); % std
            taxis = handles.dt*(1:numel(meanvals));
            if ~isempty(intersect(y_ind,[10 15])) % nucleus or cell size
                meanvals = meanvals*(handles.pixelsize)^2;
                stdvals = stdvals*(handles.pixelsize)^2;        
            end
            if 5==y_ind
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

% may be, may be not ?
[c_min_val, c_max_val] = visualization_range(handles,c_ind);
%                            
try % calm down if there is no data,     
    c = colorbar(handles.time_plot_axes,'TickLabels',{linspace(c_min_val,c_max_val,11)});        
    caxis(handles.time_plot_axes,[c_min_val,c_max_val]);
    cmap = jet(256);       
    colormap(handles.time_plot_axes,cmap);                
    c.Label.String = handles.features(c_ind);            
catch
end
%
% try ??
% axis(handles.time_plot_axes,[min(tb_data) max(te_data) min(y_data) max(y_data)]);
% catch 
% end
%
[y_min_val, y_max_val] = visualization_range(handles,y_ind);
axis(handles.time_plot_axes,[min(tb_data) max(te_data) y_min_val y_max_val]);                                      
%
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
legend(['{\Delta}N = ' num2str(max(handles.cell_nums)-min(handles.cell_nums))]);
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
if 8==size(handles.NUC_STATS,2)
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
                    case 20 % 
                        minval = min(handles.rng_nuc_solidity);
                        maxval = max(handles.rng_nuc_solidity);
                    case 21 % 
                        minval = min(handles.rng_nuc_shape_factor);
                        maxval = max(handles.rng_nuc_shape_factor);
                    case 22 % 
                        minval = min(handles.rng_nuc_eccentricity);
                        maxval = max(handles.rng_nuc_eccentricity);
                    case 23 % 
                        minval = min(handles.rng_nuc_orientation);
                        maxval = max(handles.rng_nuc_orientation);                                                
             end
                          
% --------------------------------------------------------------------
function mitotic_interval_FRET_ratio_heatmap_Callback(hObject, eventdata, handles)
% hObject    handle to mitotic_interval_FRET_ratio_heatmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if ~isfield(handles,'MI_tracks') || isempty(handles.MI_tracks), return, end    
    mode = get(handles.tracks_to_show,'Value');
    if 2~=mode
        set(handles.tracks_to_show,'Value',2); % mitotics
        tracks_to_show_Callback(hObject, eventdata, handles);
    end
    t_dependent_Nuclei_ratio_FRET_ratio_heatmapper(handles);

function single_FOVs_options(handles,flag)  % 'off' or 'on'
    set(handles.show_per_frame_mean_std,'Enable',flag);
    set(handles.update_vidi_immediately,'Enable',flag);
    
    set(handles.show_cell_numbers,'Enable',flag);
    set(handles.average_pixel_brightness,'Enable',flag);
    set(handles.visualize_selection,'Enable',flag);
    

% --------------------------------------------------------------------
function load_trackmate_plus_data_multiple_Callback(hObject, eventdata, handles)
%
single_FOVs_options(handles,'off');

[filenames,pathname] = uigetfile('*.mat','Select track data files',pwd,'MultiSelect','on');                
if isempty(filenames), return, end       
if isnumeric(filenames) && 0==filenames, return, end
if ~iscell(filenames), filenames = cellstr(filenames); end % if single FOV

handles.MI_tracks = cell(0);
handles.MI_fnames = cell(0);
handles.MI_track_indices = [];
handles.MI_norm_FRET_ratio = [];
handles.MI_norm_nuc_size = [];
handles.MI_peak_shift = [];
handles.filenames = filenames;

handles.ST_raw_data = cell(0); % "storage"
handles.ST_raw_data_filenames = cell(0); % corresponding filenames

%hw = waitbar(0,'finding mitotic intervals.. please wait');
%warning('off');
%for k=1:numel(tracks)
%if ~isempty(hw), waitbar(k/numel(tracks),hw); drawnow, end    

hw = waitbar(0,'loading multiple FOVs data.. please wait');
for k=1:numel(filenames)
    if ~isempty(hw), waitbar(k/numel(filenames),hw,['loading multiple FOVs data.. please wait : ' num2str(k) ' , ' num2str(numel(filenames))]); drawnow, end        
    load([pathname filesep filenames{k}]);
    handles.ST_raw_data = [handles.ST_raw_data; tracks];    
    handles.ST_raw_data_filenames = [handles.ST_raw_data_filenames; repmat(filenames(k),numel(tracks),1) ];

    handles.dt = dt;
    handles.pixelsize = microns_per_pixel;     
end
if ~isempty(hw), delete(hw), drawnow; end

    [handles.features_lut,handles.features_void,handles.features_coeff] = ... 
                                    set_features_to_data_correspondence(handles);
    guidata(hObject,handles);

    handles.ME = detect_mitotic_events(handles,'on'); % verbose
    %
    [MI_tracks, ... 
    MI_track_indices, ...
    MI_norm_FRET_ratio, ...
    MI_norm_nuc_size, ...
    MI_peak_shift, ...
    MI_fnames] = get_mitotic_intervals_by_using_mitotic_events(handles,'on');
    %
    if ~isempty(MI_tracks)
        handles.MI_tracks = MI_tracks;
        handles.MI_norm_FRET_ratio = MI_norm_FRET_ratio;
        handles.MI_norm_nuc_size = MI_norm_nuc_size; 
        handles.MI_fnames = MI_fnames;
        handles.MI_track_indices = MI_track_indices;
        handles.MI_peak_shift = MI_peak_shift;
    end
    
set(handles.pixel_size_edit,'String',num2str(handles.pixelsize));
set(handles.delta_t_edit,'String',num2str(handles.dt)); 

guidata(hObject, handles);

tracks_to_show_Callback(hObject, eventdata, handles);
% sic!

% --- Executes on selection change in tracks_to_show.
function tracks_to_show_Callback(hObject, eventdata, handles)

    if ~isfield(handles,'MI_tracks'), return, end    

    str = get(handles.tracks_to_show,'String');
    mode = str{get(handles.tracks_to_show,'Value')};

    if strcmp(mode,'tracks')
        handles.raw_data = handles.ST_raw_data;
    elseif strcmp(mode,'mitotic intervals')
        if ~isempty(handles.MI_tracks)
            handles.raw_data = handles.MI_tracks;
        else
            handles.raw_data = handles.ST_raw_data;
            set(handles.tracks_to_show,'Value',1); % tracks only
        end
    end

    % this object is for visualizing       
    [handles.track_data,handles.velocity_t,handles.nuc_cell_area_ratio_t] = calculate_track_data(handles);

    if isfield(handles,'filenames') && ~isempty(handles.filenames)
        handles.fullfilename = [handles.filenames{1} ' .. ' handles.filenames{numel(handles.filenames)}];
    end
    set(handles.figure1, 'Name', [handles.figureName ' : ' handles.fullfilename]);    
    
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
            
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tracks_to_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tracks_to_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function mitotic_event_detection_setup_Callback(hObject, eventdata, handles)
t_dependent_Nuclei_ratio_FRET_ME_MI_settings(handles);

% --------------------------------------------------------------------
function mitotic_interval_FRET_save_intervals_Callback(hObject, eventdata, handles)
    if isfield(handles,'MI_tracks') && ~isempty(handles.MI_tracks)    
        [fname, fpath] = uiputfile('*.mat','Save Intervals as..','_mitotic_intervals.mat');
        if fpath == 0; return; end
        filespec = fullfile(fpath,fname);
        %                
        MI_tracks = handles.MI_tracks;
        MI_norm_FRET_ratio = handles.MI_norm_FRET_ratio;
        MI_norm_nuc_size = handles.MI_norm_nuc_size; 
        MI_fnames = handles.MI_fnames;
        MI_track_indices = handles.MI_track_indices;
        MI_peak_shift = handles.MI_peak_shift;  
        dt = handles.dt;
        pixelsize = handles.pixelsize;
        %    
        save(filespec,'MI_tracks','MI_norm_FRET_ratio','MI_norm_nuc_size','MI_fnames','MI_track_indices','MI_peak_shift','dt','pixelsize');
    end 
                        
% --------------------------------------------------------------------
function HTS_heatmapper_Callback(hObject, eventdata, handles)
    %
    t_dependent_Nuclei_ratio_FRET_HTS_heatmapper(handles);
    % detect_mitotic_events(handles,'on')

% --------------------------------------------------------------------
function MI_HTS_mapper_Callback(hObject, eventdata, handles)
    %
    t_dependent_Nuclei_ratio_FRET_HTS_MI_heatmapper(handles);
    
    
function mevents = detect_mitotic_events(handles,verbose)
            %handles.ST_raw_data = [handles.ST_raw_data; tracks];    
            %handles.ST_raw_data_filenames = [handles.ST_raw_data_filenames; repmat(filenames(k),numel(tracks),1) ];

t_nucsize_drv_uncond    = handles.MED_t_nucsize_drv_uncond; % 0.4; % threshold for unconditional acceptance (derivative strength)
t_nucsize_drv_cond      = handles.MED_t_nucsize_drv_cond; % 0.25; % threshold for conditional acceptance (derivative strength)
dt_twin_max             = handles.MED_dt_twin_max; % 76; % beyond 60 minutes, can consider as next mitotic      
T_min                   = handles.MED_T_min; % 5; % don't consider signals shorter than T_min     
big_smoothing_window    = handles.MED_big_smoothing_window; % 280; % min
small_smoothing_window  = handles.MED_small_smoothing_window; % 35; %min;
t_nucsize_ratio         = handles.MED_t_nucsize_ratio; % 1.75;
DT                      = handles.MED_DT; % 10; % time interval (min) either to the past or future to calculate average nuc_size there
                        
L_min = fix(T_min/handles.dt);

dt_min = handles.dt*60; % interval between frames in minutes
% convert to frames to use in smoothing
big_smoothing_window  = round(big_smoothing_window/dt_min);
small_smoothing_window  = round(small_smoothing_window/dt_min);

%verbose = get(handles.show_per_frame_mean_std,'enable');

if strcmp('on',verbose)
    hw = waitbar(0,'finding mitotic events.. please wait');
end

warning('off');

mevents = [];

tracks = handles.ST_raw_data;
for track_ind=1:numel(tracks)
if strcmp('on',verbose) && ~isempty(hw), waitbar(track_ind/numel(tracks),hw); drawnow, end    
    track = tracks{track_ind};    
    %
    FRET_ratio = squeeze(track(:,4));
    nucleus_size = squeeze(track(:,7));
    %
     sm_s = medfilt2(FRET_ratio,[small_smoothing_window 1]);
     [~,lf_FRET] = TD_high_pass_filter(sm_s,big_smoothing_window);
     s1 = sm_s-lf_FRET;
    % repeat for nuclei size
    sm_s = medfilt2(nucleus_size,[small_smoothing_window 1]);
    [~,lf] = TD_high_pass_filter(sm_s,big_smoothing_window);
    s2 = sm_s-lf;
    % normalize signals
    s1 = (s1-mean(s1(:)))/std(s1(:));
    s2 = (s2-mean(s2(:)))/std(s2(:));  
        
    if length(s2)<=L_min, continue, end
    
    % one minute/step upsampling
    fac = round(handles.dt/(1/60));
    s1 = interp(s1,fac);
    s2 = interp(s2,fac);
    
    nucleus_size_interp = interp(nucleus_size,fac);
    nucleus_size_interp = nucleus_size_interp(2:numel(nucleus_size_interp)); % shift to compare with derivative
        
    r2 = -diff(s2); % derivative of nuc_size
    
    [pks_2,locs_2] = findpeaks(r2,'MinPeakWidth',fac,'MinPeakHeight',t_nucsize_drv_cond); % peaks of r2 derivative
            
     z2 = zeros(size(r2)); % z2 (minutes axis where ME are marked "1").. whatever - may be sizeof r2 or r3
     z2(locs_2) = 1; % assign them all as valid first, then analyze to set some of them to 0
     i2 = find(z2~=0);

        for kk=1:numel(i2)
            L2=i2(kk); % location of nuclear size derivative's jump
            if pks_2(kk)>t_nucsize_drv_uncond % unconditional acceptance
                L2_is_OK = true;
            else % maybe it will work with weak criterion
                L2_is_OK = pks_2(kk)>t_nucsize_drv_cond; % good starting point if satisfies weak criterion
                    if L2_is_OK %
                        bfre = max(1,L2-DT):L2;
                        aftr = L2:min(L2+DT,numel(nucleus_size_interp));
                        nuc_size_bfre = mean(nucleus_size_interp(bfre)); % max?
                        nuc_size_aftr = mean(nucleus_size_interp(aftr)); % min?
                        %
                        if nuc_size_bfre/nuc_size_aftr < t_nucsize_ratio % but unfortunately nuc_size differennce isn't big enough
                            L2_is_OK = false;
                        end
                    end
            end
            if ~L2_is_OK
                z2(L2) = 0; % set to 0 if satisfies neither unconditional, nor weak criterion
            end            
        end

      % delete twins from z2 (minutes axis where ME are marked "1")         
      if 0==sum(z2), continue, end
      %
      cnd = num2cell(find(z2==1)); % candidates
      exclusion = [];
      for k=2:numel(cnd)
          if cnd{k}-cnd{k-1}<dt_twin_max
              exclusion = [exclusion; k];
          end
      end
      cnd_corr = cnd;
      cnd_corr(exclusion)=[];       

% VISUALIZATION
%         z2 = zeros(size(z2));      
%         z2(cell2mat(cnd_corr))=1;                      
%         %
%         mes = find(z2~=0);
%         figure;
%         h2=subplot(2,1,1);
%         ff=1:length(s2);
%         plot(h2,ff,s1,'g:',ff,s2,'k.-','linewidth',2);
%         hold(h2,'on');
%             plot(h2,1:length(r2),r2,'m.-','linewidth',3); 
%             hold(h2,'on');
%                 plot(h2,mes,4*ones(size(mes)),'r*',locs_2,pks_2,'m*','linewidth',3);                 
%                 hold(h2,'off');                
%         grid(h2,'on');                
%             
%         h3=subplot(2,1,2);
%             me = [];
%             for k=1:numel(cnd_corr)
%                 within_track_ME_frame_index = max(1,round(cnd_corr{k}/fac));
%                 me = [me; within_track_ME_frame_index];
%             end        
%             plot(h3,1:length(nucleus_size),nucleus_size/mean(nucleus_size),'k.-',1:length(nucleus_size),FRET_ratio/mean(FRET_ratio),'c:',me,ones(size(me)),'r*','linewidth',2);
%             grid(h3,'on');
%            disp('');
% VISUALIZATION
        
        % calculate actual index in the track for every event and add to "mevents"
        for k=1:numel(cnd_corr)
            within_track_ME_frame_index = max(1,round(cnd_corr{k}/fac));
            mevents = [mevents; [track_ind within_track_ME_frame_index]];
        end    
end

size(mevents,1)

if strcmp('on',verbose) && ~isempty(hw), delete(hw), drawnow; end
warning('on');


function [tracks_out,track_indices,norm_FRET_ratio,norm_nuc_size,peak_shift,filenames_out] = ... 
                    get_mitotic_intervals_by_using_mitotic_events(handles,verbose)                
tracks_out = [];
track_indices = [];
norm_FRET_ratio = [];
norm_nuc_size = [];
peak_shift = [];
filenames_out = cell(0);

if isempty(handles.ME), return, end

big_smoothing_window    = handles.MED_big_smoothing_window; % 280; % min
small_smoothing_window  = handles.MED_small_smoothing_window; % 35; %min;

% convert to frames from minutes
dt_min = handles.dt*60; % [minutes]
big_smoothing_window  = round(big_smoothing_window/dt_min);
small_smoothing_window  = round(small_smoothing_window/dt_min);
    LLeft = round(handles.MI_LLeft/dt_min);
    LRight = round(handles.MI_LRight/dt_min);

if strcmp('on',verbose)
    hw = waitbar(0,'finding mitotic intervals.. please wait');
end

tracks = handles.ST_raw_data;
for k=1:size(handles.ME,1)
if strcmp('on',verbose) && ~isempty(hw), waitbar(k/size(handles.ME,1),hw); drawnow, end
    track_ind   = handles.ME(k,1);
    track = tracks{track_ind};
    FRET_ratio = squeeze(track(:,4));
    nucleus_size = squeeze(track(:,7));
    %
    sm_s = medfilt2(FRET_ratio,[small_smoothing_window 1]);
    [~,lf_FRET] = TD_high_pass_filter(sm_s,big_smoothing_window);
    s1 = sm_s-lf_FRET;
    % repeat for nuclei size
    sm_s = medfilt2(nucleus_size,[small_smoothing_window 1]);
    [~,lf] = TD_high_pass_filter(sm_s,big_smoothing_window);
    s2 = sm_s-lf;
    % normalize signals
    s1 = (s1-mean(s1(:)))/std(s1(:));
    s2 = (s2-mean(s2(:)))/std(s2(:));  

    mevent = handles.ME(k,2); % index of ME centre within track
    rb = mevent - LLeft;
    re = mevent + LRight;
    
    if rb>=1 && re<=length(s1)
                    tracks_out = [tracks_out; {track(rb:re,:)}];
                    track_indices = [ track_indices; track_ind];
                    norm_FRET_ratio = [norm_FRET_ratio; s1(rb:re)'];
                    norm_nuc_size = [norm_nuc_size; s2(rb:re)'];
                    peak_shift = [peak_shift; 0]; % not used
                    filenames_out = [filenames_out; handles.ST_raw_data_filenames{track_ind}];
    end    
end % loop over mitotic events
if strcmp('on',verbose) && ~isempty(hw), delete(hw), drawnow; end

%-------------------------------
function [features_lut,features_void,features_coeff] = set_features_to_data_correspondence(handles)
%%%%%%%
% handles.features = {'duration [h]', ...                   1
%                     'XY speed [um/min]', ...              2
%                     'XY directionality', ...              3
%                     '#neighbours', ...                    4
%                     'cell density', ...                   5
%                     'FRET ratio', ...                     6
%                     'FRET ratio variability', ...         7
%                     'Donor intensity',...                 8
%                     'Acc. intensity',...                  9
%                     'nucleus size [um^2]',...             10
%                     'D/A Pearson corr.',...               11
%                     't(start) [h]',...                    12
%                     'FRET ratio c-time [min]',...         13
%                     'FRET molar fraction', ...            14
%                     'Cell area (ref) [um^2]',...          15
%                     'Cell intensity (ref)',...            16
%                     'Nucleus intensity (ref)',...         17
%                     'Nuc/Cyt intensity (ref)',...         18
%                     'Nuc/Cell area (ref)',...             19
%                     'Nucleus solidity',...                20
%                     'Nucleus shape factor',...            21
%                     'Nucleus eccentricity',...            22
%                     'Nucleus orientation',...             23
features_coeff = ones(1,numel(handles.features));
features_coeff(5) = (handles.pixelsize)^(-2);
features_coeff(10) = (handles.pixelsize)^2;
features_coeff(15) = (handles.pixelsize)^2;
switch size(handles.ST_raw_data{1},2)
    case 11
        features_lut = [0 0 0 9 10 4 0 5 6 7 8 0 0 11 0 0 0 0 0 0 0 0 0];         % no 3rd channel, no nuc. shape
        features_void = 15:23;
    case 15
        features_lut = [0 0 0 9 10 4 0 5 6 7 8 0 0 11 12 13 14 15 0 0 0 0 0];     % with 3rd channel, no nuc shape
        features_void = 20:23;
    case 16
        features_lut = [0 0 0 9 10 4 0 5 6 7 8 0 0 11 0 0 0 0 0 12 13 14 15];     % no 3rd channel, with nuc shapes
        features_void = 15:19;
    case 20
        features_lut = [0 0 0 9 10 4 0 5 6 7 8 0 0 11 12 13 14 15 0 16 17 18 19]; % with 3rd channel, with nuc. shapes
        features_void = [];
end

% --------------------------------------------------------------------
function exportTrackData_Callback(hObject, eventdata, handles, should_filter)
% hObject    handle to exportTrackData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This function is designed to dump out all feature information for all
% fovs & tracks (over time). The current output is fov_track.csv

% Optional arg handling
if ~exist('should_filter', 'var')
    should_filter = false;
end

% Target save location
path = uigetdir(matlabroot, 'Export Directory');

% Quest dialog for output format
answer = questdlg('.mat or .csv output?', 'Output Format', '.mat', '.csv', '.mat');

% Header information
header = handles.features;
coeff = handles.features_coeff;
header = [{'frame_index'}, {'x'}, {'y'}, header(:)'];
coeff = [1, 1, 1, coeff(:)'];
features_lut = handles.features_lut;
features_lut(1, 1:3) = 1:3;
features_lut = nonzeros(features_lut);
present_headers = header(features_lut);
present_coeff = coeff(features_lut);

%WAAAHA :P
present_headers = strrep(present_headers, ' ', '_');
present_headers = strrep(present_headers, '[', '');
present_headers = strrep(present_headers, ']', '');
present_headers = strrep(present_headers, '/', '_per_');
present_headers = strrep(present_headers, '#', 'number_of_');
present_headers = strrep(present_headers, '.', '');
present_headers = strrep(present_headers, '^', '');
present_headers = strrep(present_headers, '(', '_');
present_headers = strrep(present_headers, ')', '');

% Add DT to our headers
present_headers = {present_headers{1}, 'dt', present_headers{2:end}, 'me_event_mask'};
present_coeff = [present_coeff(1), 1, present_coeff(2:end), 1];

% Add XY speed
add_speed = false;
if ~ismember(2, handles.features_void)
    add_speed = true;
    present_headers = [present_headers(:)', {'xy_speed'}];
    present_coeff = [present_coeff(:)', 1];
end

% Add nuclear cell area - not sure why this index is not in the
% features_lut?
add_nuc_area = false;
if ~ismember(19, handles.features_void)
    add_nuc_area = true;
    present_headers = [present_headers(:)', {'nuclear_area_ratio'}];
    present_coeff = [present_coeff(:)', 1];
end

% Track counter dict
track_counter = containers.Map;

% Loop through our tracks but retain source information
mat_output = struct();
count = numel(handles.ST_raw_data_filenames);
h1 = waitbar(0, "Exporting data - please wait");
for index = 1:count
    waitbar(index/count, h1, "Exporting data - please wait");
    track_source = handles.ST_raw_data_filenames{index, 1};
    [~, track_source, ~] = fileparts(track_source);
    
    % Track counter
    if ismember(track_source, keys(track_counter))
        track_counter(track_source) = track_counter(track_source) + 1;
    else
        track_counter(track_source) = 1;
    end
    
    % Get the data
    track_data = handles.raw_data{index, 1};    
    
    % Add dt
    track_data = [track_data(:, 1), track_data(:, 1) * handles.dt, track_data(:, 2:end)];
    
    % Add the MEs
    loc = find(handles.ME(:, 1)==index);
    binary_vector = zeros(size(track_data, 1), 1);
    if ~isnan(loc)
        binary_vector(handles.ME(loc, 2), 1) = 1;
    end
    track_data = [track_data(:, 1:end), binary_vector];
    
    % Optionally add speed
    if add_speed
        track_data = [track_data(:, 1:end), handles.velocity_t{index}];
    end
    
    % Optionally inject nuclear area ratio
    if add_nuc_area
        track_data = [track_data(:, 1:end), handles.nuc_cell_area_ratio_t{index}];
    end
    
    % Handle coefficients
    track_data = track_data .* present_coeff;
    
    % Optional filtering of outputs based on MEs
    if should_filter && sum(binary_vector) == 0
        continue
    end
    
    % Handle output format options
    track_name = strcat(track_source, '_track_', num2str(track_counter(track_source)));
    table = array2table(track_data, 'VariableNames', present_headers);
    switch answer
        
        % Save as mat
        case '.mat'
            track_name = genvarname(track_name);
            mat_output.(track_name) = table;
            
        % Save as csv
        case '.csv'
            % Save file
            filepath = strcat(path, filesep, track_name, '.csv');
            writetable(table, filepath);
    end
end
close(h1);

% If we want as a mat
if strcmp('.mat', answer)
    filepath = strcat(path, filesep, 'FRET_TrackPlotter_Output_', date, '.mat');
    save(filepath, 'mat_output');
end

% --------------------------------------------------------------------
function exportFiltered_Callback(hObject, eventdata, handles)
% hObject    handle to exportFiltered (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
exportTrackData_Callback(hObject, eventdata, handles, true);
