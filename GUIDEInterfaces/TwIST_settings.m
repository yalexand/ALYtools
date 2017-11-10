function varargout = TwIST_settings(varargin)
% TWIST_SETTINGS MATLAB code for TwIST_settings.fig
%      TWIST_SETTINGS, by itself, creates a new TWIST_SETTINGS or raises the existing
%      singleton*.
%
%      H = TWIST_SETTINGS returns the handle to a new TWIST_SETTINGS or the handle to
%      the existing singleton*.
%
%      TWIST_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TWIST_SETTINGS.M with the given input arguments.
%
%      TWIST_SETTINGS('Property','Value',...) creates a new TWIST_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TwIST_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TwIST_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TwIST_settings

% Last Modified by GUIDE v2.5 23-Jan-2015 16:16:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TwIST_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @TwIST_settings_OutputFcn, ...
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


% --- Executes just before TwIST_settings is made visible.
function TwIST_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TwIST_settings (see VARARGIN)

data_controller = varargin{1};
handles.data_controller = data_controller;

        handles.TwIST_TAU = data_controller.TwIST_TAU;
        handles.TwIST_LAMBDA = data_controller.TwIST_LAMBDA;
        handles.TwIST_ALPHA = data_controller.TwIST_ALPHA;
        handles.TwIST_BETA = data_controller.TwIST_BETA;
        handles.TwIST_STOPCRITERION = data_controller.TwIST_STOPCRITERION;
        handles.TwIST_TOLERANCEA = data_controller.TwIST_TOLERANCEA;
        handles.TwIST_TOLERANCED = data_controller.TwIST_TOLERANCED;
        handles.TwIST_DEBIAS = data_controller.TwIST_DEBIAS;
        handles.TwIST_MAXITERA = data_controller.TwIST_MAXITERA;
        handles.TwIST_MAXITERD = data_controller.TwIST_MAXITERD;
        handles.TwIST_MINITERA = data_controller.TwIST_MINITERA;
        handles.TwIST_MINITERD = data_controller.TwIST_MINITERD;
        handles.TwIST_INITIALIZATION = data_controller.TwIST_INITIALIZATION;
        handles.TwIST_MONOTONE = data_controller.TwIST_MONOTONE;
        handles.TwIST_SPARSE = data_controller.TwIST_SPARSE;
        handles.TwIST_VERBOSE = data_controller.TwIST_VERBOSE;

        set(handles.edit1,'String',num2str(handles.TwIST_TAU));        
        set(handles.edit2,'String',num2str(handles.TwIST_LAMBDA));        
        set(handles.edit3,'String',num2str(handles.TwIST_ALPHA));        
        set(handles.edit4,'String',num2str(handles.TwIST_BETA));        
        set(handles.edit5,'String',num2str(handles.TwIST_STOPCRITERION));        
        set(handles.edit6,'String',num2str(handles.TwIST_TOLERANCEA));        
        set(handles.edit7,'String',num2str(handles.TwIST_TOLERANCED));        
        set(handles.edit8,'String',num2str(handles.TwIST_DEBIAS));        
        set(handles.edit9,'String',num2str(handles.TwIST_MAXITERA));        
        set(handles.edit10,'String',num2str(handles.TwIST_MAXITERD));        
        set(handles.edit11,'String',num2str(handles.TwIST_MINITERA));        
        set(handles.edit12,'String',num2str(handles.TwIST_MINITERD));        
        set(handles.edit13,'String',num2str(handles.TwIST_INITIALIZATION));        
        set(handles.edit14,'String',num2str(handles.TwIST_MONOTONE));        
        set(handles.edit15,'String',num2str(handles.TwIST_SPARSE));        
        set(handles.edit16,'String',num2str(handles.TwIST_VERBOSE));                

% Choose default command line output for TwIST_settings
handles.output = hObject;

output_value = 1;

handles.output = output_value;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TwIST_settings wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TwIST_settings_OutputFcn(hObject, eventdata, handles) 
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

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.TwIST_TAU = value;
    guidata(hObject,handles);
else
    value = handles.TwIST_TAU;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.TwIST_LAMBDA = value;
    guidata(hObject,handles);
else
    value = handles.TwIST_LAMBDA;
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
    handles.TwIST_ALPHA = value;
    guidata(hObject,handles);
else
    value = handles.TwIST_ALPHA;
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
    handles.TwIST_BETA = value;
    guidata(hObject,handles);
else
    value = handles.TwIST_BETA;
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
    handles.TwIST_STOPCRITERION = value;
    guidata(hObject,handles);
else
    value = handles.TwIST_STOPCRITERION;
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
    handles.TwIST_TOLERANCEA = value;
    guidata(hObject,handles);
else
    value = handles.TwIST_TOLERANCEA;
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
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.TwIST_TOLERANCED = value;
    guidata(hObject,handles);
else
    value = handles.TwIST_TOLERANCED;
    set(hObject,'String',num2str(value));
end
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



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.TwIST_DEBIAS = value;
    guidata(hObject,handles);
else
    value = handles.TwIST_DEBIAS;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


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



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.TwIST_MAXITERA = value;
    guidata(hObject,handles);
else
    value = handles.TwIST_MAXITERA;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


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



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.TwIST_MAXITERD = value;
    guidata(hObject,handles);
else
    value = handles.TwIST_MAXITERD;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


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



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.TwIST_MINITERA = value;
    guidata(hObject,handles);
else
    value = handles.TwIST_MINITERD;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


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



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.TwIST_MINITERD = value;
    guidata(hObject,handles);
else
    value = handles.TwIST_MINITERD;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


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



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.TwIST_INITIALIZATION = value;
    guidata(hObject,handles);
else
    value = handles.TwIST_INITIALIZATION;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


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



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.TwIST_MONOTONE = value;
    guidata(hObject,handles);
else
    value = handles.TwIST_MONOTONE;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


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



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.TwIST_SPARSE = value;
    guidata(hObject,handles);
else
    value = handles.TwIST_SPARSE;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


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



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.TwIST_VERBOSE = value;
    guidata(hObject,handles);
else
    value = handles.TwIST_VERBOSE;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


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
        data_controller.TwIST_TAU = handles.TwIST_TAU;
        data_controller.TwIST_LAMBDA = handles.TwIST_LAMBDA;
        data_controller.TwIST_ALPHA = handles.TwIST_ALPHA;
        data_controller.TwIST_BETA = handles.TwIST_BETA;
        data_controller.TwIST_STOPCRITERION = handles.TwIST_STOPCRITERION;
        data_controller.TwIST_TOLERANCEA = handles.TwIST_TOLERANCEA;
        data_controller.TwIST_TOLERANCED = handles.TwIST_TOLERANCED;
        data_controller.TwIST_DEBIAS = handles.TwIST_DEBIAS;
        data_controller.TwIST_MAXITERA = handles.TwIST_MAXITERA;
        data_controller.TwIST_MAXITERD = handles.TwIST_MAXITERD;
        data_controller.TwIST_MINITERA = handles.TwIST_MINITERA;
        data_controller.TwIST_MINITERD = handles.TwIST_MINITERD;
        data_controller.TwIST_INITIALIZATION = handles.TwIST_INITIALIZATION;
        data_controller.TwIST_MONOTONE = handles.TwIST_MONOTONE;
        data_controller.TwIST_SPARSE = handles.TwIST_SPARSE;
        data_controller.TwIST_VERBOSE = handles.TwIST_VERBOSE;
    fh = ancestor(hObject,'figure');     
    delete(fh);

% Cancel    
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    fh = ancestor(hObject,'figure');     
    delete(fh);

% Set default
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data_controller = handles.data_controller;
data_controller.set_TwIST_settings_default;
    fh = ancestor(hObject,'figure');     
    delete(fh);


