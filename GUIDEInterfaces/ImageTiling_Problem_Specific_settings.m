function varargout = ImageTiling_Problem_Specific_settings(varargin)
% IMAGETILING_PROBLEM_SPECIFIC_SETTINGS MATLAB code for ImageTiling_Problem_Specific_settings.fig
%      IMAGETILING_PROBLEM_SPECIFIC_SETTINGS, by itself, creates a new IMAGETILING_PROBLEM_SPECIFIC_SETTINGS or raises the existing
%      singleton*.
%
%      H = IMAGETILING_PROBLEM_SPECIFIC_SETTINGS returns the handle to a new IMAGETILING_PROBLEM_SPECIFIC_SETTINGS or the handle to
%      the existing singleton*.
%
%      IMAGETILING_PROBLEM_SPECIFIC_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGETILING_PROBLEM_SPECIFIC_SETTINGS.M with the given input arguments.
%
%      IMAGETILING_PROBLEM_SPECIFIC_SETTINGS('Property','Value',...) creates a new IMAGETILING_PROBLEM_SPECIFIC_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImageTiling_Problem_Specific_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImageTiling_Problem_Specific_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ImageTiling_Problem_Specific_settings

% Last Modified by GUIDE v2.5 15-Feb-2019 10:01:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImageTiling_Problem_Specific_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @ImageTiling_Problem_Specific_settings_OutputFcn, ...
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


% --- Executes just before ImageTiling_Problem_Specific_settings is made visible.
function ImageTiling_Problem_Specific_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImageTiling_Problem_Specific_settings (see VARARGIN)
data_controller = varargin{1};
if isempty(data_controller.imgdata) && isempty(data_controller.M_imgdata), 
    errordlg('no data loaded - can not continue'), 
    fh = ancestor(hObject,'figure');     
    delete(fh);
    return, 
end;    

handles.dc = data_controller;
        % 
        handles.ImageTiling_Ncols = handles.dc.ImageTiling_Ncols;
        handles.ImageTiling_Nrows = handles.dc.ImageTiling_Nrows;
        handles.ImageTiling_Ovlp_X = handles.dc.ImageTiling_Ovlp_X;
        handles.ImageTiling_Ovlp_Y = handles.dc.ImageTiling_Ovlp_Y;
        handles.ImageTiling_QT = handles.dc.ImageTiling_QT;
        handles.ImageTiling_mode = handles.dc.ImageTiling_mode;
                        
        set(handles.cols_edit,'String',num2str(handles.ImageTiling_Ncols)); 
        set(handles.rows_edit,'String',num2str(handles.ImageTiling_Nrows)); 
        set(handles.ovlp_x_edit,'String',num2str(handles.ImageTiling_Ovlp_X));
        set(handles.ovlp_y_edit,'String',num2str(handles.ImageTiling_Ovlp_Y));
        set(handles.QT_edit,'String',num2str(handles.ImageTiling_QT));
        
        mode = handles.dc.ImageTiling_mode;
        str = {'bleached_fluor','brightfield'};
        set(handles.modality,'String',str);
        if strcmp(char(str(1)),mode)
            set(handles.modality,'Value',1);
        else
            set(handles.modality,'Value',2);
        end
        handles.ImageTiling_mode = mode;
            
% Choose default command line output for ImageTiling_Problem_Specific_settings
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ImageTiling_Problem_Specific_settings wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ImageTiling_Problem_Specific_settings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function cols_edit_Callback(hObject, eventdata, handles)
value = str2double(get(hObject,'String'));
if isnumeric(value) && value >=1 
    handles.ImageTiling_Ncols = value;
    guidata(hObject,handles);
else
    value = handles.ImageTiling_Ncols;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function cols_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cols_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rows_edit_Callback(hObject, eventdata, handles)
value = str2double(get(hObject,'String'));
if isnumeric(value) && value>=1 
    handles.ImageTiling_Nrows = value;
    guidata(hObject,handles);
else
    value = handles.ImageTiling_Nrows;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function rows_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rows_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ovlp_x_edit_Callback(hObject, eventdata, handles)
value = str2double(get(hObject,'String'));
if isnumeric(value) && 0<value&&value<1 
    handles.ImageTiling_Ovlp_X = value;
    guidata(hObject,handles);
else
    value = handles.ImageTiling_Ovlp_X;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function ovlp_x_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ovlp_x_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ovlp_y_edit_Callback(hObject, eventdata, handles)
value = str2double(get(hObject,'String'));
if isnumeric(value) && 0<value&&value<1 
    handles.ImageTiling_Ovlp_Y = value;
    guidata(hObject,handles);
else
    value = handles.ImageTiling_Ovlp_Y;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function ovlp_y_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ovlp_y_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function QT_edit_Callback(hObject, eventdata, handles)
value = str2double(get(hObject,'String'));
if isnumeric(value) && value>=0 
    handles.ImageTiling_QT = value;
    guidata(hObject,handles);
else
    value = handles.ImageTiling_QT;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function QT_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to QT_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in modality.
function modality_Callback(hObject, eventdata, handles)
    index = get(hObject,'Value');
    str = get(hObject,'String');
    handles.ImageTiling_mode = str(index);
    guidata(hObject,handles);
    uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function modality_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in OK.
function OK_Callback(hObject, eventdata, handles)
dc = handles.dc;
        dc.ImageTiling_Ncols = handles.ImageTiling_Ncols;
        dc.ImageTiling_Nrows = handles.ImageTiling_Nrows;
        dc.ImageTiling_Ovlp_X = handles.ImageTiling_Ovlp_X;
        dc.ImageTiling_Ovlp_Y = handles.ImageTiling_Ovlp_Y;
        dc.ImageTiling_QT = handles.ImageTiling_QT;
        dc.ImageTiling_mode = handles.ImageTiling_mode;
fh = ancestor(hObject,'figure');     
    delete(fh);

% --- Executes on button press in Cancel.
function Cancel_Callback(hObject, eventdata, handles)
fh = ancestor(hObject,'figure');     
    delete(fh);
