function varargout = PR_Segmentation_settings(varargin)
% PR_SEGMENTATION_SETTINGS MATLAB code for PR_Segmentation_settings.fig
%      PR_SEGMENTATION_SETTINGS, by itself, creates a new PR_SEGMENTATION_SETTINGS or raises the existing
%      singleton*.
%
%      H = PR_SEGMENTATION_SETTINGS returns the handle to a new PR_SEGMENTATION_SETTINGS or the handle to
%      the existing singleton*.
%
%      PR_SEGMENTATION_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PR_SEGMENTATION_SETTINGS.M with the given input arguments.
%
%      PR_SEGMENTATION_SETTINGS('Property','Value',...) creates a new PR_SEGMENTATION_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PR_Segmentation_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PR_Segmentation_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PR_Segmentation_settings

% Last Modified by GUIDE v2.5 07-Aug-2015 08:57:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PR_Segmentation_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @PR_Segmentation_settings_OutputFcn, ...
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


% --- Executes just before PR_Segmentation_settings is made visible.
function PR_Segmentation_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PR_Segmentation_settings (see VARARGIN)

% data_controller = varargin{1};
% handles.data_controller = data_controller;
% 
% %26,25,28,27,30,29

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_controller = varargin{1};
handles.data_controller = data_controller;

        handles.PR_K = data_controller.PR_K;
        handles.PR_S1 = data_controller.PR_S1;
        handles.PR_S2 = data_controller.PR_S2;
        handles.PR_a = data_controller.PR_a;        
        handles.PR_t = data_controller.PR_t;        
        handles.PR_min_size = data_controller.PR_min_size;        
        handles.PR_mode = data_controller.PR_mode;                
 
        handles.saved_PR_K = data_controller.PR_K;
        handles.saved_PR_S1 = data_controller.PR_S1;
        handles.saved_PR_S2 = data_controller.PR_S2;
        handles.saved_PR_a = data_controller.PR_a;        
        handles.saved_PR_t = data_controller.PR_t;        
        handles.saved_PR_min_size = data_controller.PR_min_size;        
        handles.saved_PR_mode = data_controller.PR_mode;                    
        
        set(handles.edit26,'String',num2str(handles.PR_S1));
        set(handles.edit25,'String',num2str(handles.PR_S2));
        set(handles.edit28,'String',num2str(handles.PR_K));
        set(handles.edit27,'String',num2str(handles.PR_a));                
        set(handles.edit30,'String',num2str(handles.PR_t));                
        set(handles.edit29,'String',num2str(handles.PR_min_size));
        
        if strcmp(handles.PR_mode,'Peak')
            set(handles.popupmenu1,'Value',1);
        else
            set(handles.popupmenu1,'Value',2);
        end
                   
% Choose default command line output for PR_Segmentation_settings
handles.output = hObject;

output_value = 1;

handles.output = output_value;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PR_Segmentation_settings wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PR_Segmentation_settings_OutputFcn(hObject, eventdata, handles) 
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
         data_controller.PR_K = handles.PR_K;
         data_controller.PR_S1 = handles.PR_S1;
         data_controller.PR_S2 = handles.PR_S2;
         data_controller.PR_a = handles.PR_a;        
         data_controller.PR_t = handles.PR_t;        
         data_controller.PR_min_size = handles.PR_min_size ;        
         data_controller.PR_mode = handles.PR_mode;                
                  
    fh = ancestor(hObject,'figure');     
    delete(fh);

% Cancel    
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    data_controller = handles.data_controller;
         data_controller.PR_K = handles.saved_PR_K;
         data_controller.PR_S1 = handles.saved_PR_S1;
         data_controller.PR_S2 = handles.saved_PR_S2;
         data_controller.PR_a = handles.saved_PR_a;        
         data_controller.PR_t = handles.saved_PR_t;        
         data_controller.PR_min_size = handles.saved_PR_min_size ;        
         data_controller.PR_mode = handles.saved_PR_mode;          
    fh = ancestor(hObject,'figure');     
    delete(fh);

% Set default
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data_controller = handles.data_controller;
data_controller.Set_PR_segmentation_default;
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
    handles.PR_S1 = value;
    guidata(hObject,handles);
else
    value = handles.PR_S1;
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
    handles.PR_S2 = value;
    guidata(hObject,handles);
else
    value = handles.PR_S2;
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
         data_controller.PR_K = handles.PR_K;
         data_controller.PR_S1 = handles.PR_S1;
         data_controller.PR_S2 = handles.PR_S2;
         data_controller.PR_a = handles.PR_a;        
         data_controller.PR_t = handles.PR_t;        
         data_controller.PR_min_size = handles.PR_min_size ;        
         data_controller.PR_mode = handles.PR_mode;                   
         %
         send_adjustment_to_Icy = true;
         if strcmp(data_controller.application,'per_image_TCSPC_FLIM') % parasite
             data_controller.do_per_image_TCSPC_FLIM_Segmentation(send_adjustment_to_Icy);
         else % normal
            data_controller.do_PR_Segmentation(send_adjustment_to_Icy);
         end

function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.PR_a = value;
    guidata(hObject,handles);
else
    value = handles.PR_a;
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



function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit28 as text
%        str2double(get(hObject,'String')) returns contents of edit28 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.PR_K = value;
    guidata(hObject,handles);
else
    value = handles.PR_K;
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
if ~isnan(value)
    handles.PR_min_size = value;
    guidata(hObject,handles);
else
    value = handles.PR_min_size;
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
if ~isnan(value)
    handles.PR_t = value;
    guidata(hObject,handles);
else
    value = handles.PR_t;
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

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
ind = get(hObject,'Value');
if ind == 1 
    handles.PR_mode = 'Peak';
else
    handles.PR_mode = 'Ridge';
end
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
