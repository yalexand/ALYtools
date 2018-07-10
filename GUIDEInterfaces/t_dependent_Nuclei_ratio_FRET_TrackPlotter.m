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

% Last Modified by GUIDE v2.5 09-Jul-2018 14:16:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
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
                    't(start) [h]'
                    };

handles.mask = [];    
handles.track_data = [];
                
if 0==nargin-3
    handles.raw_data = [];
    handles.dt = 1/12; % 5 minutes
    handles.pixelsize = 2; % microns
elseif 3==nargin-3
    handles.raw_data = varargin{1};
    handles.dt = varargin{2};
    handles.pixelsize = varargin{3};
    handles.track_data = calculate_track_data(hObject,handles);
end

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
       for k=1:12
            minmaxlimits(k+1,1)=min(squeeze(handles.track_data(:,k+2)));
            minmaxlimits(k+1,2)=max(squeeze(handles.track_data(:,k+2)));
       end
handles.mask = calculate_mask(hObject,handles);        
end

set(handles.filter_table, 'Data', minmaxlimits);
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
handles.track_data = calculate_track_data(hObject,handles);
minmaxlimits(1,1)=min(squeeze(handles.track_data(1,:)));
minmaxlimits(1,2)=max(squeeze(handles.track_data(2,:)));
       for k=1:12
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
handles.track_data = calculate_track_data(hObject,handles);
minmaxlimits(1,1)=min(squeeze(handles.track_data(1,:)));
minmaxlimits(1,2)=max(squeeze(handles.track_data(2,:)));
       for k=1:12
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

% --------------------------------------------------------------------
function load_trackmate_plus_data_Callback(hObject, eventdata, handles)
% hObject    handle to load_trackmate_plus_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [filename,pathname] = uigetfile({'*.mat','TrackMate+ Files'}, ...
                    'Select data file',pwd);
    if filename == 0, return, end;
    load([pathname filesep filename]);

    if ~exist('microns_per_pixel','var'), return, end
    
    handles.raw_data = tracks;
    handles.dt = dt;
    handles.pixelsize = microns_per_pixel;
    
    handles.track_data = calculate_track_data(hObject,handles);
       
    set(handles.figure1, 'Name', [handles.figureName ' : ' filename]);

    str = [{'time'} handles.features];
    minmaxlimits = zeros(numel(str),2);
    minmaxlimits(:,1)=-Inf;
    minmaxlimits(:,2)=Inf;
    if ~isempty(handles.track_data)
    minmaxlimits(1,1)=min(squeeze(handles.track_data(:,1)));
    minmaxlimits(1,2)=max(squeeze(handles.track_data(:,2)));
        for k=1:12
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
function track_data = calculate_track_data(hObject,handles)
D = handles.raw_data;
if isempty(D), return, end;

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
    mean_FRET_ratio = mean(FRET_ratio);
    mean_donor_intensity = mean(donor_intensity);
    mean_acceptor_intensity = mean(acceptor_intensity);
    mean_nucleus_size = mean(nucleus_size);
    mean_Pearson_correlation = mean(Pearson_correlation); 
    %
    track_data(k,6+2) = mean_FRET_ratio;
    %
    %smoothed_FRET_ratio = medfilt2(FRET_ratio,[5 1]);
    %track_data(k,7+2) = norm(smoothed_FRET_ratio - mean_FRET_ratio);
    %
%     figure(22);
%     if length(FRET_ratio)>100 && length(FRET_ratio)<600
%         s = FRET_ratio;
%         W = 40;
%         sm_s = medfilt2(FRET_ratio,[5 1]);
%         [~,lf] = TD_high_pass_filter(sm_s,W);        
%         subplot(2,1,1);
%         plot(1:length(s),s,'b.-',1:length(s),lf,'r-',1:length(s),sm_s,'k-');
%         grid on;
%         subplot(2,1,2);
%         plot(1:length(s),sm_s-lf,'b.-');
%         axis([1 600 -0.4 0.4]);
%         grid on;
%         disp(k);
%     end
        W = 40;
        sm_s = medfilt2(FRET_ratio,[5 1]);
        [~,lf] = TD_high_pass_filter(sm_s,W);        
        s = sm_s-lf;
        track_data(k,7+2) = sqrt(mean(s.*s));
    %
    
    %
    %quantify tracks XY trajectory
    x = squeeze(track(:,2));
    y = squeeze(track(:,3));
    z = zeros(size(x));
    [directionality,velocity,velocity_sd] = quantify_track(x',y',z','noZ');
    track_data(k,2+2) = velocity*handles.pixelsize/(handles.dt*60); % :0
    track_data(k,3+2) = directionality;
    %
    track_data(k,8+2) = mean_donor_intensity;
    track_data(k,9+2) = mean_acceptor_intensity;
    track_data(k,10+2) = mean_nucleus_size*(handles.pixelsize)^2;
    track_data(k,11+2) = mean_Pearson_correlation;
    track_data(k,12+2) = track_data(k,1); % start time
end

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

if strcmp(mode,'scatter')
    plot(handles.histo2_axes,x_data,y_data,'r.');
    grid(handles.histo2_axes,'on');
elseif strcmp(mode,'histo2')
    corr_map_W = 100;
    corrmap = correlation_map(x_data,y_data,corr_map_W);       
    AXES = handles.histo2_axes;
    imagesc(corrmap,'Parent',AXES);
    daspect(AXES,[1 1 1]);
    set(AXES, 'xticklabel', [], 'yticklabel', []);
end

% --------------------------------------------------------------------
function visualize_time_dependence(hObject,handles)
D = handles.track_data;
if isempty(D), return, end;

% need to use filtered data here
c_ind = get(handles.time_plot_colour_feature,'Value');
y_ind = get(handles.time_plot_Y_feature,'Value');
c_data = squeeze(D(:,c_ind+2));
y_data = squeeze(D(:,y_ind+2));
tb_data = squeeze(D(:,1));
te_data = squeeze(D(:,2));

mask = handles.mask;

% define colors
Ngrades = 256;
Colors = jet(Ngrades);
min_val = min(c_data);
max_val = max(c_data);
    
axes(handles.time_plot_axes);
for k = 1:numel(y_data)
    if 1==mask(k)
        index = max(1,round((c_data(k) - min_val)/(max_val - min_val)*Ngrades));
        plot(handles.time_plot_axes,[tb_data(k) te_data(k)],[y_data(k) y_data(k)],'Color',Colors(index,:));
        hold(handles.time_plot_axes,'on');
    end
end
plot(handles.time_plot_axes,tb_data(mask==1),y_data(mask==1),'k.');
hold(handles.time_plot_axes,'off');

c = colorbar(handles.time_plot_axes,'Ticks',linspace(min_val,max_val,10));
c.Label.String = handles.features(c_ind);

axis(handles.time_plot_axes,[min(tb_data) max(te_data) min(y_data) max(y_data)]);
xlabel(handles.time_plot_axes,'time [h]');
%str = get(handles.time_plot_Y_feature,'String')
%ylabel(handles.time_plot_axes,str{y_ind});
grid(handles.time_plot_axes,'on');

% --------------------------------------------------------------------
function mask = calculate_mask(hObject,handles)
    D = handles.track_data;
    if isempty(D), return, end;

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
