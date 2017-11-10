function varargout = TTO_Problem_Specific_settings(varargin)
% TTO_PROBLEM_SPECIFIC_SETTINGS MATLAB code for TTO_Problem_Specific_settings.fig
%      TTO_PROBLEM_SPECIFIC_SETTINGS, by itself, creates a new TTO_PROBLEM_SPECIFIC_SETTINGS or raises the existing
%      singleton*.
%
%      H = TTO_PROBLEM_SPECIFIC_SETTINGS returns the handle to a new TTO_PROBLEM_SPECIFIC_SETTINGS or the handle to
%      the existing singleton*.
%
%      TTO_PROBLEM_SPECIFIC_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TTO_PROBLEM_SPECIFIC_SETTINGS.M with the given input arguments.
%
%      TTO_PROBLEM_SPECIFIC_SETTINGS('Property','Value',...) creates a new TTO_PROBLEM_SPECIFIC_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TTO_Problem_Specific_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TTO_Problem_Specific_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TTO_Problem_Specific_settings

% Last Modified by GUIDE v2.5 01-Aug-2015 11:27:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TTO_Problem_Specific_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @TTO_Problem_Specific_settings_OutputFcn, ...
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




% --- Executes just before TTO_Problem_Specific_settings is made visible.
function TTO_Problem_Specific_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TTO_Problem_Specific_settings (see VARARGIN)

data_controller = varargin{1};
if isempty(data_controller.imgdata), 
    errordlg('no data loaded - can not continue'), 
    fh = ancestor(hObject,'figure');     
    delete(fh);
    return, 
end;

handles.data_controller = data_controller;
        %   
        handles.TTO_ref_channel = data_controller.TTO_ref_channel;        
        handles.TTO_Lmin = data_controller.TTO_Lmin;
        handles.TTO_Lmax = data_controller.TTO_Lmax;
        handles.TTO_Nmax = data_controller.TTO_Nmax;        
        handles.TTO_fL = data_controller.TTO_fL;
        handles.TTO_fH = data_controller.TTO_fH;
        handles.TTO_df = data_controller.TTO_df;
        
        % to be used on Cancel
        handles.saved_TTO_ref_channel = data_controller.TTO_ref_channel;        
        handles.saved_TTO_Lmin = data_controller.TTO_Lmin;
        handles.saved_TTO_Lmax = data_controller.TTO_Lmax;
        handles.saved_TTO_Nmax = data_controller.TTO_Nmax;        
        handles.saved_TTO_fL = data_controller.TTO_fL;
        handles.saved_TTO_fH = data_controller.TTO_fH;
        handles.saved_TTO_df = data_controller.TTO_df;        
                
        set(handles.edit3,'String',num2str(handles.TTO_ref_channel)); 
        set(handles.edit27,'String',num2str(handles.TTO_Lmin)); 
        set(handles.edit28,'String',num2str(handles.TTO_Lmax)); 
        set(handles.edit29,'String',num2str(handles.TTO_Nmax));        
        set(handles.edit30,'String',num2str(handles.TTO_fL)); 
        set(handles.edit31,'String',num2str(handles.TTO_fH)); 
        set(handles.edit32,'String',num2str(handles.TTO_df)); 
                        
% Choose default command line output for TTO_Problem_Specific_settings
handles.output = hObject;

output_value = 1;

handles.output = output_value;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TTO_Problem_Specific_settings wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TTO_Problem_Specific_settings_OutputFcn(hObject, eventdata, handles) 
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

%         TTO_ref_channel = 1; % edit3
%         TTO_Lmin = 180; % edit27
%         TTO_Lmax = 300; % edit28
%         TTO_Nmax = 40000; % edit29
%         TTO_fL = 100; % edit30
%         TTO_fH = 6; % edit31
%         TTO_df = 2; % edit32

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
    handles.TTO_ref_channel = fix(value);
    guidata(hObject,handles);
else
    value = handles.TTO_ref_channel;
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
        data_controller.TTO_ref_channel = handles.TTO_ref_channel;        
        data_controller.TTO_Lmin = handles.TTO_Lmin;
        data_controller.TTO_Lmax = handles.TTO_Lmax;
        data_controller.TTO_Nmax = handles.TTO_Nmax;        
        data_controller.TTO_fL = handles.TTO_fL;
        data_controller.TTO_fH = handles.TTO_fH;
        data_controller.TTO_df = handles.TTO_df;
    fh = ancestor(hObject,'figure');     
    delete(fh);

% Cancel    
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    data_controller = handles.data_controller;
        data_controller.TTO_ref_channel = handles.saved_TTO_ref_channel;        
        data_controller.TTO_Lmin = handles.saved_TTO_Lmin;
        data_controller.TTO_Lmax = handles.saved_TTO_Lmax;
        data_controller.TTO_Nmax = handles.saved_TTO_Nmax;        
        data_controller.TTO_fL = handles.saved_TTO_fL;
        data_controller.TTO_fH = handles.saved_TTO_fH;
        data_controller.TTO_df = handles.saved_TTO_df;    
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
    handles.TTO_Lmin = value;
    guidata(hObject,handles);
else
    value = handles.TTO_Lmin;
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
    handles.TTO_Lmax = value;
    guidata(hObject,handles);
else
    value = handles.TTO_Lmax;
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
    handles.TTO_Nmax = value;
    guidata(hObject,handles);
else
    value = handles.TTO_Nmax;
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
    handles.TTO_fL = value;
    guidata(hObject,handles);
else
    value = handles.TTO_fL;
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


function edit31_Callback(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit31 as text
%        str2double(get(hObject,'String')) returns contents of edit31 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.TTO_fH = value;
    guidata(hObject,handles);
else
    value = handles.TTO_fH;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit31_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit32_Callback(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit32 as text
%        str2double(get(hObject,'String')) returns contents of edit32 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.TTO_df = value;
    guidata(hObject,handles);
else
    value = handles.TTO_df;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function edit32_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
