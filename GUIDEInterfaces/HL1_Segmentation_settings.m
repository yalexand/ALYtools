function varargout = HL1_Segmentation_settings(varargin)
% HL1_SEGMENTATION_SETTINGS MATLAB code for HL1_Segmentation_settings.fig
%      HL1_SEGMENTATION_SETTINGS, by itself, creates a new HL1_SEGMENTATION_SETTINGS or raises the existing
%      singleton*.
%
%      H = HL1_SEGMENTATION_SETTINGS returns the handle to a new HL1_SEGMENTATION_SETTINGS or the handle to
%      the existing singleton*.
%
%      HL1_SEGMENTATION_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HL1_SEGMENTATION_SETTINGS.M with the given input arguments.
%
%      HL1_SEGMENTATION_SETTINGS('Property','Value',...) creates a new HL1_SEGMENTATION_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HL1_Segmentation_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HL1_Segmentation_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HL1_Segmentation_settings

% Last Modified by GUIDE v2.5 16-Sep-2015 11:08:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HL1_Segmentation_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @HL1_Segmentation_settings_OutputFcn, ...
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


% --- Executes just before HL1_Segmentation_settings is made visible.
function HL1_Segmentation_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HL1_Segmentation_settings (see VARARGIN)

data_controller = varargin{1};
handles.data_controller = data_controller;

% %26,25,28,popupmenu1
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.HL1_sgm_minimal_radius = data_controller.HL1_sgm_minimal_radius;
        handles.HL1_sgm_stdT_threshold = data_controller.HL1_sgm_stdT_threshold;
        handles.HL1_sgm_I_threshold = data_controller.HL1_sgm_I_threshold;
        handles.HL1_sgm_mode = data_controller.HL1_sgm_mode;        
        
        handles.saved_HL1_sgm_minimal_radius = data_controller.HL1_sgm_minimal_radius;
        handles.saved_HL1_sgm_stdT_threshold = data_controller.HL1_sgm_stdT_threshold;
        handles.saved_HL1_sgm_I_threshold = data_controller.HL1_sgm_I_threshold;
        handles.saved_HL1_sgm_mode = data_controller.HL1_sgm_mode;        
                
        set(handles.edit26,'String',num2str(handles.HL1_sgm_minimal_radius));
        set(handles.edit25,'String',num2str(handles.HL1_sgm_I_threshold));
        set(handles.edit28,'String',num2str(handles.HL1_sgm_stdT_threshold));
        
        % NEED TO SET UP THE MODE!!!
        mode_index = 1;
        switch handles.HL1_sgm_mode
            case char(data_controller.HL1_sgm_modes{1})
                mode_index = 1;
            case char(data_controller.HL1_sgm_modes{2})
                mode_index = 2;
            case char(data_controller.HL1_sgm_modes{3})
                mode_index = 3;
            case char(data_controller.HL1_sgm_modes{4})
                mode_index = 4;                
        end
        set(handles.popupmenu1,'Value',mode_index);
                                        
% Choose default command line output for HL1_Segmentation_settings
handles.output = hObject;

output_value = 1;

handles.output = output_value;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HL1_Segmentation_settings wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HL1_Segmentation_settings_OutputFcn(hObject, eventdata, handles) 
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
        data_controller.HL1_sgm_minimal_radius = handles.HL1_sgm_minimal_radius;
        data_controller.HL1_sgm_stdT_threshold = handles.HL1_sgm_stdT_threshold;
        data_controller.HL1_sgm_I_threshold = handles.HL1_sgm_I_threshold;
        data_controller.HL1_sgm_mode = handles.HL1_sgm_mode;        

    fh = ancestor(hObject,'figure');     
    delete(fh);

% Cancel    
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
     data_controller = handles.data_controller;
        data_controller.HL1_sgm_minimal_radius = handles.saved_HL1_sgm_minimal_radius;
        data_controller.HL1_sgm_stdT_threshold = handles.saved_HL1_sgm_stdT_threshold;
        data_controller.HL1_sgm_I_threshold = handles.saved_HL1_sgm_I_threshold;
        data_controller.HL1_sgm_mode = handles.saved_HL1_sgm_mode;        
    fh = ancestor(hObject,'figure');     
    delete(fh);


function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.HL1_sgm_minimal_radius = value;
    guidata(hObject,handles);
else
    value = handles.HL1_sgm_minimal_radius;
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

function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.HL1_sgm_I_threshold = value;
    guidata(hObject,handles);
else
    value = handles.HL1_sgm_I_threshold;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
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
        data_controller.HL1_sgm_minimal_radius = handles.HL1_sgm_minimal_radius;
        data_controller.HL1_sgm_stdT_threshold = handles.HL1_sgm_stdT_threshold;
        data_controller.HL1_sgm_I_threshold = handles.HL1_sgm_I_threshold;
        data_controller.HL1_sgm_mode = handles.HL1_sgm_mode;        
        %
        send_adjustment_to_Icy = true;
        data_controller.do_HL1_Segmentation(send_adjustment_to_Icy);

function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit28 as text
%        str2double(get(hObject,'String')) returns contents of edit28 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.HL1_sgm_stdT_threshold = value;
    guidata(hObject,handles);
else
    value = handles.HL1_sgm_stdT_threshold;
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


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
data_controller = handles.data_controller;
handles.HL1_sgm_mode = char(data_controller.HL1_sgm_modes{get(hObject,'Value')});
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
