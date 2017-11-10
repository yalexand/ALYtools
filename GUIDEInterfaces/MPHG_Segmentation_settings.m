function varargout = MPHG_Segmentation_settings(varargin)
% MPHG_SEGMENTATION_SETTINGS MATLAB code for MPHG_Segmentation_settings.fig
%      MPHG_SEGMENTATION_SETTINGS, by itself, creates a new MPHG_SEGMENTATION_SETTINGS or raises the existing
%      singleton*.
%
%      H = MPHG_SEGMENTATION_SETTINGS returns the handle to a new MPHG_SEGMENTATION_SETTINGS or the handle to
%      the existing singleton*.
%
%      MPHG_SEGMENTATION_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MPHG_SEGMENTATION_SETTINGS.M with the given input arguments.
%
%      MPHG_SEGMENTATION_SETTINGS('Property','Value',...) creates a new MPHG_SEGMENTATION_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MPHG_Segmentation_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MPHG_Segmentation_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MPHG_Segmentation_settings

% Last Modified by GUIDE v2.5 27-Jan-2016 10:04:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MPHG_Segmentation_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @MPHG_Segmentation_settings_OutputFcn, ...
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


% --- Executes just before MPHG_Segmentation_settings is made visible.
function MPHG_Segmentation_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MPHG_Segmentation_settings (see VARARGIN)

% data_controller = varargin{1};
% handles.data_controller = data_controller;
% 
% %26,25,28,27,30,29

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_controller = varargin{1};
handles.data_controller = data_controller;

        handles.MPHG_nuc_scale = data_controller.MPHG_nuc_scale;
        handles.MPHG_nuc_rel_bg_scale = data_controller.MPHG_nuc_rel_bg_scale;
        handles.MPHG_nuc_threshold = data_controller.MPHG_nuc_threshold;
        handles.MPHG_nuc_smoothing_scale = data_controller.MPHG_nuc_smoothing_scale;        
        handles.MPHG_nuc_min_area = data_controller.MPHG_nuc_min_area;        
        handles.MPHG_cell_smoothing_radius = data_controller.MPHG_cell_smoothing_radius;        
        handles.MPHG_cell_rg_std_factor_int = data_controller.MPHG_cell_rg_std_factor_int; 
        handles.MPHG_cell_rg_std_factor_out = data_controller.MPHG_cell_rg_std_factor_out;        
        handles.MPHG_cell_overpeak_ratio = data_controller.MPHG_cell_overpeak_ratio;
        handles.MPHG_minimal_gran_size = data_controller.MPHG_minimal_gran_size;
        handles.MPHG_gran_overpeak_ratio = data_controller.MPHG_gran_overpeak_ratio;
         
        handles.saved_MPHG_nuc_scale = data_controller.MPHG_nuc_scale;
        handles.saved_MPHG_nuc_rel_bg_scale = data_controller.MPHG_nuc_rel_bg_scale;
        handles.saved_MPHG_nuc_threshold = data_controller.MPHG_nuc_threshold;
        handles.saved_MPHG_nuc_smoothing_scale = data_controller.MPHG_nuc_smoothing_scale;        
        handles.saved_MPHG_nuc_min_area = data_controller.MPHG_nuc_min_area;        
        handles.saved_MPHG_cell_smoothing_radius = data_controller.MPHG_cell_smoothing_radius;        
        handles.saved_MPHG_cell_rg_std_factor_int = data_controller.MPHG_cell_rg_std_factor_int; 
        handles.saved_MPHG_cell_rg_std_factor_out = data_controller.MPHG_cell_rg_std_factor_out;        
        handles.saved_MPHG_cell_overpeak_ratio = data_controller.MPHG_cell_overpeak_ratio;
        handles.saved_MPHG_minimal_gran_size = data_controller.MPHG_minimal_gran_size;
        handles.saved_MPHG_gran_overpeak_ratio = data_controller.MPHG_gran_overpeak_ratio;
                
        set(handles.edit26,'String',num2str(handles.MPHG_nuc_scale));
        set(handles.edit28,'String',num2str(handles.MPHG_nuc_rel_bg_scale));
        set(handles.edit30,'String',num2str(handles.MPHG_nuc_threshold));
        set(handles.edit33,'String',num2str(handles.MPHG_nuc_smoothing_scale));                
        set(handles.edit29,'String',num2str(handles.MPHG_nuc_min_area));                
        set(handles.edit35,'String',num2str(handles.MPHG_cell_smoothing_radius));
        set(handles.edit36,'String',num2str(handles.MPHG_cell_rg_std_factor_int));                
        set(handles.edit34,'String',num2str(handles.MPHG_cell_rg_std_factor_out));                
        set(handles.edit38,'String',num2str(handles.MPHG_cell_overpeak_ratio));
        set(handles.edit37,'String',num2str(handles.MPHG_minimal_gran_size));
        set(handles.edit39,'String',num2str(handles.MPHG_gran_overpeak_ratio));
                            
% Choose default command line output for MPHG_Segmentation_settings
handles.output = hObject;

output_value = 1;

handles.output = output_value;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MPHG_Segmentation_settings wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MPHG_Segmentation_settings_OutputFcn(hObject, eventdata, handles) 
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

% OK
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data_controller = handles.data_controller;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        data_controller.MPHG_nuc_scale = handles.MPHG_nuc_scale;
        data_controller.MPHG_nuc_rel_bg_scale = handles.MPHG_nuc_rel_bg_scale;
        data_controller.MPHG_nuc_threshold = handles.MPHG_nuc_threshold;
        data_controller.MPHG_nuc_smoothing_scale = handles.MPHG_nuc_smoothing_scale;        
        data_controller.MPHG_nuc_min_area = handles.MPHG_nuc_min_area;        
        data_controller.MPHG_cell_smoothing_radius = handles.MPHG_cell_smoothing_radius;        
        data_controller.MPHG_cell_rg_std_factor_int = handles.MPHG_cell_rg_std_factor_int; 
        data_controller.MPHG_cell_rg_std_factor_out = handles.MPHG_cell_rg_std_factor_out;        
        data_controller.MPHG_cell_overpeak_ratio = handles.MPHG_cell_overpeak_ratio; 
        data_controller.MPHG_minimal_gran_size = handles.MPHG_minimal_gran_size;
        data_controller.MPHG_gran_overpeak_ratio = handles.MPHG_gran_overpeak_ratio;
                          
    fh = ancestor(hObject,'figure');     
    delete(fh);

% Cancel    
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    data_controller = handles.data_controller;
        data_controller.MPHG_nuc_scale = handles.saved_MPHG_nuc_scale;
        data_controller.MPHG_nuc_rel_bg_scale = handles.saved_MPHG_nuc_rel_bg_scale;
        data_controller.MPHG_nuc_threshold = handles.saved_MPHG_nuc_threshold;
        data_controller.MPHG_nuc_smoothing_scale = handles.saved_MPHG_nuc_smoothing_scale;        
        data_controller.MPHG_nuc_min_area = handles.saved_MPHG_nuc_min_area;        
        data_controller.MPHG_cell_smoothing_radius = handles.saved_MPHG_cell_smoothing_radius;        
        data_controller.MPHG_cell_rg_std_factor_int = handles.saved_MPHG_cell_rg_std_factor_int; 
        data_controller.MPHG_cell_rg_std_factor_out = handles.saved_MPHG_cell_rg_std_factor_out;        
        data_controller.MPHG_cell_overpeak_ratio = handles.saved_MPHG_cell_overpeak_ratio;
        data_controller.MPHG_minimal_gran_size = handles.saved_MPHG_minimal_gran_size;
        data_controller.MPHG_gran_overpeak_ratio = handles.saved_MPHG_gran_overpeak_ratio;
        
    fh = ancestor(hObject,'figure');     
    delete(fh);

% Set default
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data_controller = handles.data_controller;
    fh = ancestor(hObject,'figure');     
    delete(fh);


function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value > 0
    handles.MPHG_nuc_scale = value;
    guidata(hObject,handles);
else
    value = handles.MPHG_nuc_scale;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data_controller = handles.data_controller;
        data_controller.MPHG_nuc_scale = handles.MPHG_nuc_scale;
        data_controller.MPHG_nuc_rel_bg_scale = handles.MPHG_nuc_rel_bg_scale;
        data_controller.MPHG_nuc_threshold = handles.MPHG_nuc_threshold;
        data_controller.MPHG_nuc_smoothing_scale = handles.MPHG_nuc_smoothing_scale;        
        data_controller.MPHG_nuc_min_area = handles.MPHG_nuc_min_area;        
        data_controller.MPHG_cell_smoothing_radius = handles.MPHG_cell_smoothing_radius;        
        data_controller.MPHG_cell_rg_std_factor_int = handles.MPHG_cell_rg_std_factor_int; 
        data_controller.MPHG_cell_rg_std_factor_out = handles.MPHG_cell_rg_std_factor_out;        
        data_controller.MPHG_cell_overpeak_ratio = handles.MPHG_cell_overpeak_ratio;
        data_controller.MPHG_minimal_gran_size = handles.MPHG_minimal_gran_size;
        data_controller.MPHG_gran_overpeak_ratio = handles.MPHG_gran_overpeak_ratio;        
        %
        send_adjustment_to_Icy = true;
        data_controller.do_MPHG_Segmentation(send_adjustment_to_Icy);


function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit28 as text
%        str2double(get(hObject,'String')) returns contents of edit28 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value > 1
    handles.MPHG_nuc_rel_bg_scale = value;
    guidata(hObject,handles);
else
    value = handles.MPHG_nuc_rel_bg_scale;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
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
if ~isnan(value) && value > 1
    handles.MPHG_nuc_min_area = value;
    guidata(hObject,handles);
else
    value = handles.MPHG_nuc_min_area;
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



function edit30_Callback(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit30 as text
%        str2double(get(hObject,'String')) returns contents of edit30 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value > 0
    handles.MPHG_nuc_threshold = value;
    guidata(hObject,handles);
else
    value = handles.MPHG_nuc_threshold;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit30 (see GCBO)
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
if ~isnan(value) && value > 0
    handles.MPHG_nuc_smoothing_scale = value;
    guidata(hObject,handles);
else
    value = handles.MPHG_nuc_smoothing_scale;
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



function edit34_Callback(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit34 as text
%        str2double(get(hObject,'String')) returns contents of edit34 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value > 0
    handles.MPHG_cell_rg_std_factor_out = value;
    guidata(hObject,handles);
else
    value = handles.MPHG_cell_rg_std_factor_out;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit34_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value > 0
    handles.MPHG_cell_smoothing_radius = value;
    guidata(hObject,handles);
else
    value = handles.MPHG_cell_smoothing_radius;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit36_Callback(hObject, eventdata, handles)
% hObject    handle to edit36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit36 as text
%        str2double(get(hObject,'String')) returns contents of edit36 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value > 0
    handles.MPHG_cell_rg_std_factor_int = value;
    guidata(hObject,handles);
else
    value = handles.MPHG_cell_rg_std_factor_int;
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
value = str2double(get(hObject,'String'));
if ~isnan(value) && value > 0
    handles.MPHG_minimal_gran_size = value;
    guidata(hObject,handles);
else
    value = handles.MPHG_minimal_gran_size;
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



function edit38_Callback(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit38 as text
%        str2double(get(hObject,'String')) returns contents of edit38 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 0
    handles.MPHG_cell_overpeak_ratio = value;
    guidata(hObject,handles);
else
    value = handles.MPHG_cell_overpeak_ratio;
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
if ~isnan(value) && value >= 0
    handles.MPHG_gran_overpeak_ratio = value;
    guidata(hObject,handles);
else
    value = handles.MPHG_gran_overpeak_ratio;
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
