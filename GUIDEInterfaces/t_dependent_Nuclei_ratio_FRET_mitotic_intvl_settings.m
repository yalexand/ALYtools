function varargout = t_dependent_Nuclei_ratio_FRET_mitotic_intvl_settings(varargin)
% T_DEPENDENT_NUCLEI_RATIO_FRET_MITOTIC_INTVL_SETTINGS MATLAB code for t_dependent_Nuclei_ratio_FRET_mitotic_intvl_settings.fig
%      T_DEPENDENT_NUCLEI_RATIO_FRET_MITOTIC_INTVL_SETTINGS, by itself, creates a new T_DEPENDENT_NUCLEI_RATIO_FRET_MITOTIC_INTVL_SETTINGS or raises the existing
%      singleton*.
%
%      H = T_DEPENDENT_NUCLEI_RATIO_FRET_MITOTIC_INTVL_SETTINGS returns the handle to a new T_DEPENDENT_NUCLEI_RATIO_FRET_MITOTIC_INTVL_SETTINGS or the handle to
%      the existing singleton*.
%
%      T_DEPENDENT_NUCLEI_RATIO_FRET_MITOTIC_INTVL_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in T_DEPENDENT_NUCLEI_RATIO_FRET_MITOTIC_INTVL_SETTINGS.M with the given input arguments.
%
%      T_DEPENDENT_NUCLEI_RATIO_FRET_MITOTIC_INTVL_SETTINGS('Property','Value',...) creates a new T_DEPENDENT_NUCLEI_RATIO_FRET_MITOTIC_INTVL_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before t_dependent_Nuclei_ratio_FRET_mitotic_intvl_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to t_dependent_Nuclei_ratio_FRET_mitotic_intvl_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help t_dependent_Nuclei_ratio_FRET_mitotic_intvl_settings

% Last Modified by GUIDE v2.5 23-Apr-2020 16:05:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @t_dependent_Nuclei_ratio_FRET_mitotic_intvl_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @t_dependent_Nuclei_ratio_FRET_mitotic_intvl_settings_OutputFcn, ...
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


% --- Executes just before t_dependent_Nuclei_ratio_FRET_mitotic_intvl_settings is made visible.
function t_dependent_Nuclei_ratio_FRET_mitotic_intvl_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to t_dependent_Nuclei_ratio_FRET_mitotic_intvl_settings (see VARARGIN)

handles.TrackPlotter_handles = varargin{1};

% handles.bckp_MI_LLeft = handles.TrackPlotter_handles.MI_LLeft;
% handles.bckp_MI_LRight = handles.TrackPlotter_handles.MI_LRight;
% handles.bckp_MI_peak_proximity_tol = handles.TrackPlotter_handles.MI_peak_proximity_tol;
% handles.bckp_MI_large_smoothing_window = handles.TrackPlotter_handles.MI_large_smoothing_window;
% handles.bckp_MI_small_smoothing_window = handles.TrackPlotter_handles.MI_small_smoothing_window;
% handles.bckp_MI_FRET_peak_tol = handles.TrackPlotter_handles.MI_FRET_peak_tol;
% handles.bckp_MI_nucsize_peak_tol = handles.TrackPlotter_handles.MI_nucsize_peak_tol;

handles.MI_LLeft = handles.TrackPlotter_handles.MI_LLeft;
handles.MI_LRight = handles.TrackPlotter_handles.MI_LRight;
handles.MI_peak_proximity_tol = handles.TrackPlotter_handles.MI_peak_proximity_tol;
handles.MI_large_smoothing_window = handles.TrackPlotter_handles.MI_large_smoothing_window;
handles.MI_small_smoothing_window = handles.TrackPlotter_handles.MI_small_smoothing_window;
handles.MI_FRET_peak_tol = handles.TrackPlotter_handles.MI_FRET_peak_tol;
handles.MI_nucsize_peak_tol = handles.TrackPlotter_handles.MI_nucsize_peak_tol;

set(handles.mi_left_margin,'String',num2str(handles.MI_LLeft));
set(handles.mi_right_margin,'String',num2str(handles.MI_LRight));
set(handles.mi_peak_proximity_tol,'String',num2str(handles.MI_peak_proximity_tol));
set(handles.mi_large_smoothing_window,'String',num2str(handles.MI_large_smoothing_window));
set(handles.mi_small_smoothing_window,'String',num2str(handles.MI_small_smoothing_window));
set(handles.mi_FRET_peak_tol,'String',num2str(handles.MI_FRET_peak_tol));
set(handles.mi_nucsize_peak_tol,'String',num2str(handles.MI_nucsize_peak_tol));

% Choose default command line output for t_dependent_Nuclei_ratio_FRET_mitotic_intvl_settings
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes t_dependent_Nuclei_ratio_FRET_mitotic_intvl_settings wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = t_dependent_Nuclei_ratio_FRET_mitotic_intvl_settings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function mi_left_margin_Callback(hObject, eventdata, handles)
    value = fix(str2double(get(hObject,'String')));
    if ~isnan(value) && value>30 && value<360
        handles.MI_LLeft = value;
        guidata(hObject,handles);
    else
        value = handles.MI_LLeft;
        set(hObject,'String',num2str(value));
    end
    uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function mi_left_margin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mi_left_margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mi_right_margin_Callback(hObject, eventdata, handles)
    value = fix(str2double(get(hObject,'String')));
    if ~isnan(value) && value>30 && value<360
        handles.MI_LRight = value;
        guidata(hObject,handles);
    else
        value = handles.LRight;
        set(hObject,'String',num2str(value));
    end
    uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function mi_right_margin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mi_right_margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mi_peak_proximity_tol_Callback(hObject, eventdata, handles)
    value = fix(str2double(get(hObject,'String')));
    if ~isnan(value) && value>1 && value<60
        handles.MI_peak_proximity_tol = value;
        guidata(hObject,handles);
    else
        value = handles.MI_peak_proximity_tol;
        set(hObject,'String',num2str(value));
    end
    uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function mi_peak_proximity_tol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mi_peak_proximity_tol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mi_large_smoothing_window_Callback(hObject, eventdata, handles)
    value = fix(str2double(get(hObject,'String')));
    if ~isnan(value) && value>60 && value<12000
        handles.MI_large_smoothing_window = value;
        guidata(hObject,handles);
    else
        value = handles.MI_large_smoothing_window;
        set(hObject,'String',num2str(value));
    end
    uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function mi_large_smoothing_window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mi_large_smoothing_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mi_small_smoothing_window_Callback(hObject, eventdata, handles)
    value = fix(str2double(get(hObject,'String')));
    if ~isnan(value) && value>10 && value<90
        handles.MI_small_smoothing_window = value;
        guidata(hObject,handles);
    else
        value = handles.MI_small_smoothing_window;
        set(hObject,'String',num2str(value));
    end
    uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function mi_small_smoothing_window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mi_small_smoothing_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mi_FRET_peak_tol_Callback(hObject, eventdata, handles)
    value = str2double(get(hObject,'String'));
    if ~isnan(value) && value>0 && value<2
        handles.MI_FRET_peak_tol = value;
        guidata(hObject,handles);
    else
        value = handles.MI_FRET_peak_tol;
        set(hObject,'String',num2str(value));
    end
    uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function mi_FRET_peak_tol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mi_FRET_peak_tol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mi_nucsize_peak_tol_Callback(hObject, eventdata, handles)
    value = str2double(get(hObject,'String'));
    if ~isnan(value) && value>0 && value<2
        handles.MI_nucsize_peak_tol = value;
        guidata(hObject,handles);
    else
        value = handles.MI_nucsize_peak_tol;
        set(hObject,'String',num2str(value));
    end
    uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function mi_nucsize_peak_tol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mi_nucsize_peak_tol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% saving results back to plotter to reuse other time
FIGS = findobj(0, 'type', 'figure');
for k=1:size(FIGS,1)
    h = FIGS(k,1);
    cur_handles = guidata(h);
    try
    if strcmp(cur_handles.timestamp,handles.TrackPlotter_handles.timestamp) % mine
        cur_handles.MI_LLeft = handles.MI_LLeft;
        cur_handles.MI_LRight = handles.MI_LRight;
        cur_handles.MI_peak_proximity_tol = handles.MI_peak_proximity_tol;
        cur_handles.MI_large_smoothing_window = handles.MI_large_smoothing_window;
        cur_handles.MI_small_smoothing_window = handles.MI_small_smoothing_window;
        cur_handles.MI_FRET_peak_tol = handles.MI_FRET_peak_tol;
        cur_handles.MI_nucsize_peak_tol = handles.MI_nucsize_peak_tol;
        guidata(h,cur_handles);        
    end
    catch
    end
end
    fh = ancestor(hObject,'figure');     
    delete(fh);

% --- Executes on button press in pushbutton_Cancel.
function pushbutton_Cancel_Callback(hObject, eventdata, handles)
    % gou out without saving
    fh = ancestor(hObject,'figure');     
    delete(fh);