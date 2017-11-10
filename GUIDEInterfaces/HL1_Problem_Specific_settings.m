function varargout = HL1_Problem_Specific_settings(varargin)
% HL1_PROBLEM_SPECIFIC_SETTINGS MATLAB code for HL1_Problem_Specific_settings.fig
%      HL1_PROBLEM_SPECIFIC_SETTINGS, by itself, creates a new HL1_PROBLEM_SPECIFIC_SETTINGS or raises the existing
%      singleton*.
%
%      H = HL1_PROBLEM_SPECIFIC_SETTINGS returns the handle to a new HL1_PROBLEM_SPECIFIC_SETTINGS or the handle to
%      the existing singleton*.
%
%      HL1_PROBLEM_SPECIFIC_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HL1_PROBLEM_SPECIFIC_SETTINGS.M with the given input arguments.
%
%      HL1_PROBLEM_SPECIFIC_SETTINGS('Property','Value',...) creates a new HL1_PROBLEM_SPECIFIC_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HL1_Problem_Specific_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HL1_Problem_Specific_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HL1_Problem_Specific_settings

% Last Modified by GUIDE v2.5 05-Oct-2015 09:25:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HL1_Problem_Specific_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @HL1_Problem_Specific_settings_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    try
    gui_mainfcn(gui_State, varargin{:});
    catch
    end
end
% End initialization code - DO NOT EDIT




% --- Executes just before HL1_Problem_Specific_settings is made visible.
function HL1_Problem_Specific_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HL1_Problem_Specific_settings (see VARARGIN)

data_controller = varargin{1};
if isempty(data_controller.imgdata), 
    errordlg('no data loaded - can not continue'), 
    fh = ancestor(hObject,'figure');     
    delete(fh);
    return, 
end;

handles.data_controller = data_controller;
        % 
                handles.HL1_ref_channel = data_controller.HL1_ref_channel; 
                handles.HL1_TD_avreraging_size = data_controller.HL1_TD_avreraging_size;
                handles.HL1_TD_smoothing_size = data_controller.HL1_TD_smoothing_size;
                handles.HL1_df = data_controller.HL1_df;
                handles.HL1_min_std_T = data_controller.HL1_min_std_T; 
                handles.HL1_dynamic_amplitude_threshold = data_controller.HL1_dynamic_amplitude_threshold; 
                handles.HL1_calculate_phasemap = data_controller.HL1_calculate_phasemap;
                handles.HL1_binning_radius = data_controller.HL1_binning_radius;
                handles.HL1_invert_input = data_controller.HL1_invert_input;
                handles.HL1_MINPEAKDISTANCE = data_controller.HL1_MINPEAKDISTANCE;
                handles.HL1_MINPEAKHEIGHT_std_factor = data_controller.HL1_MINPEAKHEIGHT_std_factor;
                %
                handles.HL1_delineate_wave_front = data_controller.HL1_delineate_wave_front;
                handles.HL1_create_LF_movie = data_controller.HL1_create_LF_movie;
                handles.HL1_create_frame_indexed_excited_regions_movie = data_controller.HL1_create_frame_indexed_excited_regions_movie; %isochrones
                handles.H1_isochrones_map_step = data_controller.H1_isochrones_map_step;    
    
                
                handles.saved_HL1_ref_channel = data_controller.HL1_ref_channel; 
                handles.saved_HL1_TD_avreraging_size = data_controller.HL1_TD_avreraging_size;
                handles.saved_HL1_TD_smoothing_size = data_controller.HL1_TD_smoothing_size;
                handles.saved_HL1_df = data_controller.HL1_df;
                handles.saved_HL1_min_std_T = data_controller.HL1_min_std_T; 
                handles.saved_HL1_dynamic_amplitude_threshold = data_controller.HL1_dynamic_amplitude_threshold; 
                handles.saved_HL1_calculate_phasemap = data_controller.HL1_calculate_phasemap;
                handles.saved_HL1_binning_radius = data_controller.HL1_binning_radius;
                handles.saved_HL1_invert_input = data_controller.HL1_invert_input;
                handles.saved_HL1_MINPEAKDISTANCE = data_controller.HL1_MINPEAKDISTANCE;
                handles.saved_HL1_MINPEAKHEIGHT_std_factor = data_controller.HL1_MINPEAKHEIGHT_std_factor;
                %
                handles.saved_HL1_delineate_wave_front = data_controller.HL1_delineate_wave_front;
                handles.saved_HL1_create_LF_movie = data_controller.HL1_create_LF_movie;
                handles.saved_HL1_create_frame_indexed_excited_regions_movie = data_controller.HL1_create_frame_indexed_excited_regions_movie; %isochrones
                handles.saved_H1_isochrones_map_step = data_controller.H1_isochrones_map_step;    
                                                
                set(handles.edit35,'String',num2str(handles.HL1_ref_channel));
                set(handles.edit3,'String',num2str(handles.HL1_TD_avreraging_size)); 
                set(handles.edit27,'String',num2str(handles.HL1_df));
                set(handles.edit29,'String',num2str(handles.HL1_min_std_T));
                set(handles.edit33,'String',num2str(handles.HL1_dynamic_amplitude_threshold));
                set(handles.edit36,'String',num2str(handles.HL1_binning_radius));
                set(handles.checkbox1,'Value',handles.HL1_calculate_phasemap);
                set(handles.edit37,'String',num2str(handles.HL1_TD_smoothing_size));
                set(handles.checkbox2,'Value',handles.HL1_invert_input);
                set(handles.edit38,'String',num2str(handles.HL1_MINPEAKDISTANCE));
                set(handles.edit39,'String',num2str(handles.HL1_MINPEAKHEIGHT_std_factor));
                %
                set(handles.checkbox6,'Value',handles.HL1_delineate_wave_front);
                set(handles.checkbox5,'Value',handles.HL1_create_LF_movie);
                set(handles.checkbox4,'Value',handles.HL1_create_frame_indexed_excited_regions_movie);
                set(handles.edit42,'String',num2str(handles.H1_isochrones_map_step));
                                                                        
% Choose default command line output for HL1_Problem_Specific_settings
handles.output = hObject;

output_value = 1;

handles.output = output_value;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HL1_Problem_Specific_settings wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HL1_Problem_Specific_settings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

% varargout{1} = 1.;
% if isfield(handles,'output')
%     varargout{1} = handles.output;
% end;
varargout = [];
% ????????

function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
value = fix(str2double(get(hObject,'String')));
if ~isnan(value) && value >= 4
    handles.HL1_TD_avreraging_size = fix(value);
    guidata(hObject,handles);
else
    value = handles.HL1_TD_avreraging_size;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% OK
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data_controller = handles.data_controller;
                data_controller.HL1_ref_channel = handles.HL1_ref_channel;
                data_controller.HL1_TD_avreraging_size = handles.HL1_TD_avreraging_size;
                data_controller.HL1_TD_smoothing_size = handles.HL1_TD_smoothing_size;                
                data_controller.HL1_df = handles.HL1_df;
                data_controller.HL1_min_std_T = handles.HL1_min_std_T; 
                data_controller.HL1_dynamic_amplitude_threshold = handles.HL1_dynamic_amplitude_threshold; 
                data_controller.HL1_calculate_phasemap = handles.HL1_calculate_phasemap;
                data_controller.HL1_binning_radius = handles.HL1_binning_radius;
                data_controller.HL1_invert_input = handles.HL1_invert_input;
                data_controller.HL1_MINPEAKDISTANCE = handles.HL1_MINPEAKDISTANCE;
                data_controller.HL1_MINPEAKHEIGHT_std_factor = handles.HL1_MINPEAKHEIGHT_std_factor; 
                %
                data_controller.HL1_delineate_wave_front = handles.HL1_delineate_wave_front;
                data_controller.HL1_create_LF_movie = handles.HL1_create_LF_movie;
                data_controller.HL1_create_frame_indexed_excited_regions_movie = handles.HL1_create_frame_indexed_excited_regions_movie; %isochrones
                data_controller.H1_isochrones_map_step = handles.H1_isochrones_map_step;
                                
    fh = ancestor(hObject,'figure');     
    delete(fh);

% Cancel    
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    data_controller = handles.data_controller;
                data_controller.HL1_ref_channel = handles.saved_HL1_ref_channel;
                data_controller.HL1_TD_avreraging_size = handles.saved_HL1_TD_avreraging_size;
                data_controller.HL1_TD_smoothing_size = handles.saved_HL1_TD_smoothing_size;                
                data_controller.HL1_df = handles.saved_HL1_df;
                data_controller.HL1_min_std_T = handles.saved_HL1_min_std_T; 
                data_controller.HL1_dynamic_amplitude_threshold = handles.saved_HL1_dynamic_amplitude_threshold; 
                data_controller.HL1_calculate_phasemap = handles.saved_HL1_calculate_phasemap;
                data_controller.HL1_binning_radius = handles.saved_HL1_binning_radius;
                data_controller.HL1_invert_input = handles.saved_HL1_invert_input;
                data_controller.HL1_MINPEAKDISTANCE = handles.saved_HL1_MINPEAKDISTANCE;
                data_controller.HL1_MINPEAKHEIGHT_std_factor = handles.saved_HL1_MINPEAKHEIGHT_std_factor; 
                %
                data_controller.HL1_delineate_wave_front = handles.saved_HL1_delineate_wave_front;
                data_controller.HL1_create_LF_movie = handles.saved_HL1_create_LF_movie;
                data_controller.HL1_create_frame_indexed_excited_regions_movie = handles.saved_HL1_create_frame_indexed_excited_regions_movie; %isochrones
                data_controller.H1_isochrones_map_step = handles.saved_H1_isochrones_map_step;    
    fh = ancestor(hObject,'figure');     
    delete(fh);



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.HL1_df = fix(value);
    guidata(hObject,handles);
else
    value = handles.HL1_df;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit29_Callback(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit29 as text
%        str2double(get(hObject,'String')) returns contents of edit29 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.HL1_min_std_T = value;
    guidata(hObject,handles);
else
    value = handles.HL1_min_std_T;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit29_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit33_Callback(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit33 as text
%        str2double(get(hObject,'String')) returns contents of edit33 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 0
    handles.HL1_dynamic_amplitude_threshold = value;
    guidata(hObject,handles);
else
    value = handles.HL1_dynamic_amplitude_threshold;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit33_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit34 as text
%        str2double(get(hObject,'String')) returns contents of edit34 as a double
data_controller = handles.data_controller;
[~,~,~,sC,~]=size(data_controller.imgdata);
value = fix(str2double(get(hObject,'String')));
if ~isnan(value) && value >= 1 && value <= sC
    handles.HL1_ref_channel = fix(value);
    guidata(hObject,handles);
else
    value = handles.HL1_ref_channel;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
handles.HL1_calculate_phasemap = get(hObject,'Value');
guidata(hObject,handles);
uiresume(handles.figure1);



function edit36_Callback(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit36 as text
%        str2double(get(hObject,'String')) returns contents of edit36 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 0
    handles.HL1_binning_radius = value;
    guidata(hObject,handles);
else
    value = handles.HL1_binning_radius;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);



% --- Executes during object creation, after setting all properties.
function edit36_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit37_Callback(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit37 as text
%        str2double(get(hObject,'String')) returns contents of edit37 as a double
value = fix(str2double(get(hObject,'String')));
if ~isnan(value) && value >= 0
    handles.HL1_TD_smoothing_size = fix(value);
    guidata(hObject,handles);
else
    value = handles.HL1_TD_smoothing_size;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
handles.HL1_invert_input = get(hObject,'Value');
guidata(hObject,handles);
uiresume(handles.figure1);


function edit38_Callback(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit38 as text
%        str2double(get(hObject,'String')) returns contents of edit38 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value > 2
    handles.HL1_MINPEAKDISTANCE = fix(value);
    guidata(hObject,handles);
else
    value = handles.HL1_MINPEAKDISTANCE;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit39_Callback(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit39 as text
%        str2double(get(hObject,'String')) returns contents of edit39 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value > 0
    handles.HL1_MINPEAKHEIGHT_std_factor = value;
    guidata(hObject,handles);
else
    value = handles.HL1_MINPEAKHEIGHT_std_factor;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit39_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit39 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
handles.HL1_create_LF_movie = get(hObject,'Value');
guidata(hObject,handles);
uiresume(handles.figure1);


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6
handles.HL1_delineate_wave_front = get(hObject,'Value');
guidata(hObject,handles);
uiresume(handles.figure1);


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
handles.HL1_create_frame_indexed_excited_regions_movie = get(hObject,'Value');
guidata(hObject,handles);
uiresume(handles.figure1);



function edit42_Callback(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit42 as text
%        str2double(get(hObject,'String')) returns contents of edit42 as a double
data_controller = handles.data_controller;
[~,~,~,sC,~]=size(data_controller.imgdata);
value = fix(str2double(get(hObject,'String')));
if ~isnan(value) && value >= 1 && value <= 100
    handles.H1_isochrones_map_step = fix(value);
    guidata(hObject,handles);
else
    value = handles.H1_isochrones_map_step;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit42_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


