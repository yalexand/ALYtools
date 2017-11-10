function varargout = NucCyt_Segmentation_settings(varargin)
% NUCCYT_SEGMENTATION_SETTINGS MATLAB code for NucCyt_Segmentation_settings.fig
%      NUCCYT_SEGMENTATION_SETTINGS, by itself, creates a new NUCCYT_SEGMENTATION_SETTINGS or raises the existing
%      singleton*.
%
%      H = NUCCYT_SEGMENTATION_SETTINGS returns the handle to a new NUCCYT_SEGMENTATION_SETTINGS or the handle to
%      the existing singleton*.
%
%      NUCCYT_SEGMENTATION_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NUCCYT_SEGMENTATION_SETTINGS.M with the given input arguments.
%
%      NUCCYT_SEGMENTATION_SETTINGS('Property','Value',...) creates a new NUCCYT_SEGMENTATION_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NucCyt_Segmentation_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NucCyt_Segmentation_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NucCyt_Segmentation_settings

% Last Modified by GUIDE v2.5 27-Nov-2015 11:05:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NucCyt_Segmentation_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @NucCyt_Segmentation_settings_OutputFcn, ...
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


% --- Executes just before NucCyt_Segmentation_settings is made visible.
function NucCyt_Segmentation_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NucCyt_Segmentation_settings (see VARARGIN)

% data_controller = varargin{1};
% handles.data_controller = data_controller;
% 
% %26,25,28,27,30,29

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_controller = varargin{1};
handles.data_controller = data_controller;

        handles.NC_nuc_scale = data_controller.NC_nuc_scale;
        handles.NC_nuc_rel_bg_scale = data_controller.NC_nuc_rel_bg_scale;
        handles.NC_nuc_threshold = data_controller.NC_nuc_threshold;
        handles.NC_nuc_smoothing = data_controller.NC_nuc_smoothing;        
        handles.NC_nuc_min_area = data_controller.NC_nuc_min_area;        
        handles.NC_cell_smoothing_radius = data_controller.NC_cell_smoothing_radius;        
        handles.NC_cell_rg_std_factor_int = data_controller.NC_cell_rg_std_factor_int; 
        handles.NC_cell_rg_std_factor_out = data_controller.NC_cell_rg_std_factor_out;        
        handles.NC_nuc_breacking_distmap_smoothing_scale = data_controller.NC_nuc_breacking_distmap_smoothing_scale;              
        handles.NC_cell_overpeak_ratio = data_controller.NC_cell_overpeak_ratio;              
         
        handles.saved_NC_nuc_scale = data_controller.NC_nuc_scale;
        handles.saved_NC_nuc_rel_bg_scale = data_controller.NC_nuc_rel_bg_scale;
        handles.saved_NC_nuc_threshold = data_controller.NC_nuc_threshold;
        handles.saved_NC_nuc_smoothing = data_controller.NC_nuc_smoothing;        
        handles.saved_NC_nuc_min_area = data_controller.NC_nuc_min_area;        
        handles.saved_NC_cell_smoothing_radius = data_controller.NC_cell_smoothing_radius;        
        handles.saved_NC_cell_rg_std_factor_int = data_controller.NC_cell_rg_std_factor_int; 
        handles.saved_NC_cell_rg_std_factor_out = data_controller.NC_cell_rg_std_factor_out;        
        handles.saved_NC_nuc_breacking_distmap_smoothing_scale = data_controller.NC_nuc_breacking_distmap_smoothing_scale;              
        handles.saved_NC_cell_overpeak_ratio = data_controller.NC_cell_overpeak_ratio;              
                
        set(handles.edit26,'String',num2str(handles.NC_nuc_scale));
        set(handles.edit28,'String',num2str(handles.NC_nuc_rel_bg_scale));
        set(handles.edit30,'String',num2str(handles.NC_nuc_threshold));
        set(handles.edit33,'String',num2str(handles.NC_nuc_smoothing));                
        set(handles.edit29,'String',num2str(handles.NC_nuc_min_area));                
        set(handles.edit35,'String',num2str(handles.NC_cell_smoothing_radius));
        set(handles.edit36,'String',num2str(handles.NC_cell_rg_std_factor_int));                
        set(handles.edit34,'String',num2str(handles.NC_cell_rg_std_factor_out));                
        set(handles.edit37,'String',num2str(handles.NC_nuc_breacking_distmap_smoothing_scale));
        set(handles.edit38,'String',num2str(handles.NC_cell_overpeak_ratio));
                            
% Choose default command line output for NucCyt_Segmentation_settings
handles.output = hObject;

output_value = 1;

handles.output = output_value;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NucCyt_Segmentation_settings wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NucCyt_Segmentation_settings_OutputFcn(hObject, eventdata, handles) 
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
        data_controller.NC_nuc_scale = handles.NC_nuc_scale;
        data_controller.NC_nuc_rel_bg_scale = handles.NC_nuc_rel_bg_scale;
        data_controller.NC_nuc_threshold = handles.NC_nuc_threshold;
        data_controller.NC_nuc_smoothing = handles.NC_nuc_smoothing;        
        data_controller.NC_nuc_min_area = handles.NC_nuc_min_area;        
        data_controller.NC_cell_smoothing_radius = handles.NC_cell_smoothing_radius;        
        data_controller.NC_cell_rg_std_factor_int = handles.NC_cell_rg_std_factor_int; 
        data_controller.NC_cell_rg_std_factor_out = handles.NC_cell_rg_std_factor_out;        
        data_controller.NC_nuc_breacking_distmap_smoothing_scale = handles.NC_nuc_breacking_distmap_smoothing_scale;              
        data_controller.NC_cell_overpeak_ratio = handles.NC_cell_overpeak_ratio;              
                  
    fh = ancestor(hObject,'figure');     
    delete(fh);

% Cancel    
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    data_controller = handles.data_controller;
        data_controller.NC_nuc_scale = handles.saved_NC_nuc_scale;
        data_controller.NC_nuc_rel_bg_scale = handles.saved_NC_nuc_rel_bg_scale;
        data_controller.NC_nuc_threshold = handles.saved_NC_nuc_threshold;
        data_controller.NC_nuc_smoothing = handles.saved_NC_nuc_smoothing;        
        data_controller.NC_nuc_min_area = handles.saved_NC_nuc_min_area;        
        data_controller.NC_cell_smoothing_radius = handles.saved_NC_cell_smoothing_radius;        
        data_controller.NC_cell_rg_std_factor_int = handles.saved_NC_cell_rg_std_factor_int; 
        data_controller.NC_cell_rg_std_factor_out = handles.saved_NC_cell_rg_std_factor_out;        
        data_controller.NC_nuc_breacking_distmap_smoothing_scale = handles.saved_NC_nuc_breacking_distmap_smoothing_scale;              
        data_controller.NC_cell_overpeak_ratio = handles.saved_NC_cell_overpeak_ratio;              
    fh = ancestor(hObject,'figure');     
    delete(fh);

% Set default
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data_controller = handles.data_controller;
%data_controller.Set_NucCyt_segmentation_default;
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
    handles.NC_nuc_scale = value;
    guidata(hObject,handles);
else
    value = handles.NC_nuc_scale;
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
        data_controller.NC_nuc_scale = handles.NC_nuc_scale;
        data_controller.NC_nuc_rel_bg_scale = handles.NC_nuc_rel_bg_scale;
        data_controller.NC_nuc_threshold = handles.NC_nuc_threshold;
        data_controller.NC_nuc_smoothing = handles.NC_nuc_smoothing;        
        data_controller.NC_nuc_min_area = handles.NC_nuc_min_area;        
        data_controller.NC_cell_smoothing_radius = handles.NC_cell_smoothing_radius;        
        data_controller.NC_cell_rg_std_factor_int = handles.NC_cell_rg_std_factor_int; 
        data_controller.NC_cell_rg_std_factor_out = handles.NC_cell_rg_std_factor_out;        
        data_controller.NC_nuc_breacking_distmap_smoothing_scale = handles.NC_nuc_breacking_distmap_smoothing_scale;                   
        data_controller.NC_cell_overpeak_ratio = handles.NC_cell_overpeak_ratio;                   
        %
        send_adjustment_to_Icy = true;
        data_controller.do_NC_Segmentation(send_adjustment_to_Icy);


function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit28 as text
%        str2double(get(hObject,'String')) returns contents of edit28 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value > 1
    handles.NC_nuc_rel_bg_scale = value;
    guidata(hObject,handles);
else
    value = handles.NC_nuc_rel_bg_scale;
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
    handles.NC_nuc_min_area = value;
    guidata(hObject,handles);
else
    value = handles.NC_nuc_min_area;
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
    handles.NC_nuc_threshold = value;
    guidata(hObject,handles);
else
    value = handles.NC_nuc_threshold;
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
    handles.NC_nuc_smoothing = value;
    guidata(hObject,handles);
else
    value = handles.NC_nuc_smoothing;
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
    handles.NC_cell_rg_std_factor_out = value;
    guidata(hObject,handles);
else
    value = handles.NC_cell_rg_std_factor_out;
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
    handles.NC_cell_smoothing_radius = value;
    guidata(hObject,handles);
else
    value = handles.NC_cell_smoothing_radius;
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
    handles.NC_cell_rg_std_factor_int = value;
    guidata(hObject,handles);
else
    value = handles.NC_cell_rg_std_factor_int;
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
    handles.NC_nuc_breacking_distmap_smoothing_scale = value;
    guidata(hObject,handles);
else
    value = handles.NC_nuc_breacking_distmap_smoothing_scale;
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
    handles.NC_cell_overpeak_ratio = value;
    guidata(hObject,handles);
else
    value = handles.NC_cell_overpeak_ratio;
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
