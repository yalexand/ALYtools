function varargout = CIDR_Problem_Specific_settings(varargin)
% CIDR_PROBLEM_SPECIFIC_SETTINGS MATLAB code for CIDR_Problem_Specific_settings.fig
%      CIDR_PROBLEM_SPECIFIC_SETTINGS, by itself, creates a new CIDR_PROBLEM_SPECIFIC_SETTINGS or raises the existing
%      singleton*.
%
%      H = CIDR_PROBLEM_SPECIFIC_SETTINGS returns the handle to a new CIDR_PROBLEM_SPECIFIC_SETTINGS or the handle to
%      the existing singleton*.
%
%      CIDR_PROBLEM_SPECIFIC_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CIDR_PROBLEM_SPECIFIC_SETTINGS.M with the given input arguments.
%
%      CIDR_PROBLEM_SPECIFIC_SETTINGS('Property','Value',...) creates a new CIDR_PROBLEM_SPECIFIC_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CIDR_Problem_Specific_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CIDR_Problem_Specific_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CIDR_Problem_Specific_settings

% Last Modified by GUIDE v2.5 01-Aug-2015 19:09:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CIDR_Problem_Specific_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @CIDR_Problem_Specific_settings_OutputFcn, ...
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




% --- Executes just before CIDR_Problem_Specific_settings is made visible.
function CIDR_Problem_Specific_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CIDR_Problem_Specific_settings (see VARARGIN)

data_controller = varargin{1};
if isempty(data_controller.imgdata), 
    errordlg('no data loaded - can not continue'), 
    fh = ancestor(hObject,'figure');     
    delete(fh);
    return, 
end;

handles.data_controller = data_controller;
        % 
        handles.CIDR_ref_channel = data_controller.CIDR_ref_channel;        
        handles.CIDR_smoothing_scale = data_controller.CIDR_smoothing_scale;
        handles.CIDR_Nmax = data_controller.CIDR_Nmax;

        handles.saved_CIDR_ref_channel = data_controller.CIDR_ref_channel;        
        handles.saved_CIDR_smoothing_scale = data_controller.CIDR_smoothing_scale;
        handles.saved_CIDR_Nmax = data_controller.CIDR_Nmax;
                
        set(handles.edit3,'String',num2str(handles.CIDR_ref_channel)); 
        set(handles.edit27,'String',num2str(handles.CIDR_smoothing_scale)); 
        set(handles.edit29,'String',num2str(handles.CIDR_Nmax));        
                        
% Choose default command line output for CIDR_Problem_Specific_settings
handles.output = hObject;

output_value = 1;

handles.output = output_value;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CIDR_Problem_Specific_settings wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CIDR_Problem_Specific_settings_OutputFcn(hObject, eventdata, handles) 
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
data_controller = handles.data_controller;
[~,~,~,sC,~]=size(data_controller.imgdata);
value = fix(str2double(get(hObject,'String')));
if ~isnan(value) && value >= 1 && value <= sC
    handles.CIDR_ref_channel = fix(value);
    guidata(hObject,handles);
else
    value = handles.CIDR_ref_channel;
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
        data_controller.CIDR_ref_channel = handles.CIDR_ref_channel;        
        data_controller.CIDR_smoothing_scale = handles.CIDR_smoothing_scale;
        data_controller.CIDR_Nmax = handles.CIDR_Nmax;
   
    fh = ancestor(hObject,'figure');     
    delete(fh);

% Cancel    
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    data_controller = handles.data_controller;
        data_controller.CIDR_ref_channel = handles.saved_CIDR_ref_channel;        
        data_controller.CIDR_smoothing_scale = handles.saved_CIDR_smoothing_scale;
        data_controller.CIDR_Nmax = handles.saved_CIDR_Nmax;    
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
    handles.CIDR_smoothing_scale = value;
    guidata(hObject,handles);
else
    value = handles.CIDR_smoothing_scale;
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
    handles.CIDR_Nmax = value;
    guidata(hObject,handles);
else
    value = handles.CIDR_Nmax;
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
