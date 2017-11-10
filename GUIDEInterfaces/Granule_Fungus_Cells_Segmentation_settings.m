function varargout = Granule_Fungus_Cells_Segmentation_settings(varargin)
% GRANULE_FUNGUS_CELLS_SEGMENTATION_SETTINGS MATLAB code for Granule_Fungus_Cells_Segmentation_settings.fig
%      GRANULE_FUNGUS_CELLS_SEGMENTATION_SETTINGS, by itself, creates a new GRANULE_FUNGUS_CELLS_SEGMENTATION_SETTINGS or raises the existing
%      singleton*.
%
%      H = GRANULE_FUNGUS_CELLS_SEGMENTATION_SETTINGS returns the handle to a new GRANULE_FUNGUS_CELLS_SEGMENTATION_SETTINGS or the handle to
%      the existing singleton*.
%
%      GRANULE_FUNGUS_CELLS_SEGMENTATION_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRANULE_FUNGUS_CELLS_SEGMENTATION_SETTINGS.M with the given input arguments.
%
%      GRANULE_FUNGUS_CELLS_SEGMENTATION_SETTINGS('Property','Value',...) creates a new GRANULE_FUNGUS_CELLS_SEGMENTATION_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Granule_Fungus_Cells_Segmentation_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Granule_Fungus_Cells_Segmentation_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Granule_Fungus_Cells_Segmentation_settings

% Last Modified by GUIDE v2.5 27-Jun-2015 14:16:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Granule_Fungus_Cells_Segmentation_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @Granule_Fungus_Cells_Segmentation_settings_OutputFcn, ...
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


% --- Executes just before Granule_Fungus_Cells_Segmentation_settings is made visible.
function Granule_Fungus_Cells_Segmentation_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Granule_Fungus_Cells_Segmentation_settings (see VARARGIN)

data_controller = varargin{1};
handles.data_controller = data_controller;

        % NB - SIZES IN MICRONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% granules
        handles.g_scale = data_controller.g_scale;
        handles.g_threshold = data_controller.g_threshold;
        handles.g_smoothing = data_controller.g_smoothing;
        handles.g_min_area = data_controller.g_min_area;
        handles.g_rel_bg_scale = data_controller.g_rel_bg_scale;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fungus
        handles.f_scale = data_controller.f_scale;
        handles.f_threshold = data_controller.f_threshold;
        handles.f_smoothing =data_controller.f_smoothing;
        handles.f_min_area = data_controller.f_min_area;
        handles.f_rel_bg_scale = data_controller.f_rel_bg_scale;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells
        handles.c_scale = data_controller.c_scale;
        handles.c_threshold = data_controller.c_threshold;
        handles.c_smoothing = data_controller.c_smoothing;
        handles.c_min_area = data_controller.c_min_area;
        handles.c_rel_bg_scale = data_controller.c_rel_bg_scale;
   
        %
        %%%%%%%%%%%%%%%%%% saved settings - to be restored on Cancel
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% granules
        handles.saved_g_scale = data_controller.g_scale;
        handles.saved_g_threshold = data_controller.g_threshold;
        handles.saved_g_smoothing = data_controller.g_smoothing;
        handles.saved_g_min_area = data_controller.g_min_area;
        handles.saved_g_rel_bg_scale = data_controller.g_rel_bg_scale;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fungus
        handles.saved_f_scale = data_controller.f_scale;
        handles.saved_f_threshold = data_controller.f_threshold;
        handles.saved_f_smoothing =data_controller.f_smoothing;
        handles.saved_f_min_area = data_controller.f_min_area;
        handles.saved_f_rel_bg_scale = data_controller.f_rel_bg_scale;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells
        handles.saved_c_scale = data_controller.c_scale;
        handles.saved_c_threshold = data_controller.c_threshold;
        handles.saved_c_smoothing = data_controller.c_smoothing;
        handles.saved_c_min_area = data_controller.c_min_area;
        handles.saved_c_rel_bg_scale = data_controller.c_rel_bg_scale;        
        %%%%%%%%%%%%%%%%%% saved settings - to be restored on Cancel
        
                
        set(handles.edit4,'String',num2str(handles.g_scale));        
        set(handles.edit3,'String',num2str(handles.g_rel_bg_scale));        
        set(handles.edit2,'String',num2str(handles.g_threshold));
        set(handles.edit5,'String',num2str(handles.g_smoothing));
        set(handles.edit6,'String',num2str(handles.g_min_area));
        
        set(handles.edit19,'String',num2str(handles.f_scale));        
        set(handles.edit18,'String',num2str(handles.f_rel_bg_scale));        
        set(handles.edit17,'String',num2str(handles.f_threshold));
        set(handles.edit20,'String',num2str(handles.f_smoothing));
        set(handles.edit21,'String',num2str(handles.f_min_area));

        set(handles.edit24,'String',num2str(handles.c_scale));        
        set(handles.edit23,'String',num2str(handles.c_rel_bg_scale));        
        set(handles.edit22,'String',num2str(handles.c_threshold));
        set(handles.edit25,'String',num2str(handles.c_smoothing));
        set(handles.edit26,'String',num2str(handles.c_min_area));
                        
% Choose default command line output for Granule_Fungus_Cells_Segmentation_settings
handles.output = hObject;

output_value = 1;

handles.output = output_value;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Granule_Fungus_Cells_Segmentation_settings wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Granule_Fungus_Cells_Segmentation_settings_OutputFcn(hObject, eventdata, handles) 
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

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.g_threshold = value;
    guidata(hObject,handles);
else
    value = handles.g_threshold;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.g_rel_bg_scale = value;
    guidata(hObject,handles);
else
    value = handles.g_rel_bg_scale;
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


function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.g_scale = value;
    guidata(hObject,handles);
else
    value = handles.g_scale;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.g_smoothing= value;
    guidata(hObject,handles);
else
    value = handles.g_smoothing;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.g_min_area = value;
    guidata(hObject,handles);
else
    value = handles.g_min_area;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double
% value = str2double(get(hObject,'String'));
% if ~isnan(value)
%     handles.close_radius = value;
%     guidata(hObject,handles);
% else
%     value = handles.close_radius;
%     set(hObject,'String',num2str(value));
% end
uiresume(handles.figure1);



% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% granules
         data_controller.g_scale = handles.g_scale;
         data_controller.g_threshold = handles.g_threshold;
         data_controller.g_smoothing = handles.g_smoothing;
         data_controller.g_min_area = handles.g_min_area;
         data_controller.g_rel_bg_scale = handles.g_rel_bg_scale;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fungus
         data_controller.f_scale = handles.f_scale; 
         data_controller.f_threshold = handles.f_threshold;
         data_controller.f_smoothing = handles.f_smoothing;
         data_controller.f_min_area = handles.f_min_area;
         data_controller.f_rel_bg_scale = handles.f_rel_bg_scale;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells
         data_controller.c_scale = handles.c_scale;
         data_controller.c_threshold = handles.c_threshold;
         data_controller.c_smoothing = handles.c_smoothing;
         data_controller.c_min_area = handles.c_min_area;
         data_controller.c_rel_bg_scale = handles.c_rel_bg_scale;
    fh = ancestor(hObject,'figure');     
    delete(fh);

% Cancel    
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    data_controller = handles.data_controller;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% granules
         data_controller.g_scale = handles.saved_g_scale;
         data_controller.g_threshold = handles.saved_g_threshold;
         data_controller.g_smoothing = handles.saved_g_smoothing;
         data_controller.g_min_area = handles.saved_g_min_area;
         data_controller.g_rel_bg_scale = handles.saved_g_rel_bg_scale;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fungus
         data_controller.f_scale = handles.saved_f_scale; 
         data_controller.f_threshold = handles.saved_f_threshold;
         data_controller.f_smoothing = handles.saved_f_smoothing;
         data_controller.f_min_area = handles.saved_f_min_area;
         data_controller.f_rel_bg_scale = handles.saved_f_rel_bg_scale;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells
         data_controller.c_scale = handles.saved_c_scale;
         data_controller.c_threshold = handles.saved_c_threshold;
         data_controller.c_smoothing = handles.saved_c_smoothing;
         data_controller.c_min_area = handles.saved_c_min_area;
         data_controller.c_rel_bg_scale = handles.saved_c_rel_bg_scale;
        
    fh = ancestor(hObject,'figure');     
    delete(fh);

% Set default
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data_controller = handles.data_controller;
data_controller.Set_Fungus_Dependent_Granule_Release_segmentation_default;
    fh = ancestor(hObject,'figure');     
    delete(fh);



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.f_threshold = value;
    guidata(hObject,handles);
else
    value = handles.f_threshold;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.f_rel_bg_scale = value;
    guidata(hObject,handles);
else
    value = handles.f_rel_bg_scale;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.f_scale = value;
    guidata(hObject,handles);
else
    value = handles.f_scale;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.f_smoothing = value;
    guidata(hObject,handles);
else
    value = handles.f_smoothing;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.f_min_area = value;
    guidata(hObject,handles);
else
    value = handles.f_min_area;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.c_threshold = value;
    guidata(hObject,handles);
else
    value = handles.c_threshold;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.c_rel_bg_scale = value;
    guidata(hObject,handles);
else
    value = handles.c_rel_bg_scale;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.c_scale = value;
    guidata(hObject,handles);
else
    value = handles.c_scale;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.c_smoothing = value;
    guidata(hObject,handles);
else
    value = handles.c_smoothing;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.c_min_area = value;
    guidata(hObject,handles);
else
    value = handles.c_min_area;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
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

% UPDATE...
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% granules
         data_controller.g_scale = handles.g_scale;
         data_controller.g_threshold = handles.g_threshold;
         data_controller.g_smoothing = handles.g_smoothing;
         data_controller.g_min_area = handles.g_min_area;
         data_controller.g_rel_bg_scale = handles.g_rel_bg_scale;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fungus
         data_controller.f_scale = handles.f_scale; 
         data_controller.f_threshold = handles.f_threshold;
         data_controller.f_smoothing = handles.f_smoothing;
         data_controller.f_min_area = handles.f_min_area;
         data_controller.f_rel_bg_scale = handles.f_rel_bg_scale;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells
         data_controller.c_scale = handles.c_scale;
         data_controller.c_threshold = handles.c_threshold;
         data_controller.c_smoothing = handles.c_smoothing;
         data_controller.c_min_area = handles.c_min_area;
         data_controller.c_rel_bg_scale = handles.c_rel_bg_scale;
% UPDATE
send_adjustment_to_Icy = true;
data_controller.do_Fungus_Dependent_Granule_Release_segmentation(send_adjustment_to_Icy);




