function varargout = t_dependent_Nuclei_ratio_FRET_ME_MI_settings(varargin)
%T_DEPENDENT_NUCLEI_RATIO_FRET_ME_MI_SETTINGS MATLAB code file for t_dependent_Nuclei_ratio_FRET_ME_MI_settings.fig
%      T_DEPENDENT_NUCLEI_RATIO_FRET_ME_MI_SETTINGS, by itself, creates a new T_DEPENDENT_NUCLEI_RATIO_FRET_ME_MI_SETTINGS or raises the existing
%      singleton*.
%
%      H = T_DEPENDENT_NUCLEI_RATIO_FRET_ME_MI_SETTINGS returns the handle to a new T_DEPENDENT_NUCLEI_RATIO_FRET_ME_MI_SETTINGS or the handle to
%      the existing singleton*.
%
%      T_DEPENDENT_NUCLEI_RATIO_FRET_ME_MI_SETTINGS('Property','Value',...) creates a new T_DEPENDENT_NUCLEI_RATIO_FRET_ME_MI_SETTINGS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to t_dependent_Nuclei_ratio_FRET_ME_MI_settings_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      T_DEPENDENT_NUCLEI_RATIO_FRET_ME_MI_SETTINGS('CALLBACK') and T_DEPENDENT_NUCLEI_RATIO_FRET_ME_MI_SETTINGS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in T_DEPENDENT_NUCLEI_RATIO_FRET_ME_MI_SETTINGS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help t_dependent_Nuclei_ratio_FRET_ME_MI_settings

% Last Modified by GUIDE v2.5 25-Nov-2020 16:30:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @t_dependent_Nuclei_ratio_FRET_ME_MI_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @t_dependent_Nuclei_ratio_FRET_ME_MI_settings_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before t_dependent_Nuclei_ratio_FRET_ME_MI_settings is made visible.
function t_dependent_Nuclei_ratio_FRET_ME_MI_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

handles.TrackPlotter_handles = varargin{1};

handles.MI_LLeft_v                    = handles.TrackPlotter_handles.MI_LLeft;
handles.MI_LRight_v                   = handles.TrackPlotter_handles.MI_LRight;
handles.MED_t_nucsize_drv_uncond_v    = handles.TrackPlotter_handles.MED_t_nucsize_drv_uncond;
handles.MED_t_nucsize_drv_cond_v      = handles.TrackPlotter_handles.MED_t_nucsize_drv_cond;
handles.MED_dt_twin_max_v             = handles.TrackPlotter_handles.MED_dt_twin_max;
handles.MED_T_min_v                   = handles.TrackPlotter_handles.MED_T_min;
handles.MED_big_smoothing_window_v    = handles.TrackPlotter_handles.MED_big_smoothing_window;
handles.MED_small_smoothing_window_v  = handles.TrackPlotter_handles.MED_small_smoothing_window;
handles.MED_t_nucsize_ratio_v         = handles.TrackPlotter_handles.MED_t_nucsize_ratio;
handles.MED_DT_v                      = handles.TrackPlotter_handles.MED_DT;

set(handles.MI_LLeft,'String',num2str(handles.MI_LLeft_v));
set(handles.MI_LRight,'String',num2str(handles.MI_LRight_v));
set(handles.MED_t_nucsize_drv_uncond,'String',num2str(handles.MED_t_nucsize_drv_uncond_v));
set(handles.MED_t_nucsize_drv_cond,'String',num2str(handles.MED_t_nucsize_drv_cond_v));
set(handles.MED_dt_twin_max,'String',num2str(handles.MED_dt_twin_max_v));
set(handles.MED_T_min,'String',num2str(handles.MED_T_min_v));
set(handles.MED_big_smoothing_window,'String',num2str(handles.MED_big_smoothing_window_v));
set(handles.MED_small_smoothing_window,'String',num2str(handles.MED_small_smoothing_window_v));
set(handles.MED_t_nucsize_ratio,'String',num2str(handles.MED_t_nucsize_ratio_v));
set(handles.MED_DT,'String',num2str(handles.MED_DT_v));

% Choose default command line output for t_dependent_Nuclei_ratio_FRET_ME_MI_settings
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes t_dependent_Nuclei_ratio_FRET_ME_MI_settings wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = t_dependent_Nuclei_ratio_FRET_ME_MI_settings_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function MED_t_nucsize_drv_uncond_Callback(hObject, eventdata, handles)
    value = str2double(get(hObject,'String'));
    if ~isnan(value) && value>=0
        handles.MED_t_nucsize_drv_uncond_v = value;
        guidata(hObject,handles);
    else
        value = handles.MED_t_nucsize_drv_uncond_v;
    end
    set(hObject,'String',num2str(value));    
    uiresume(handles.figure1);



% --- Executes during object creation, after setting all properties.
function MED_t_nucsize_drv_uncond_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MED_t_nucsize_drv_uncond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MED_t_nucsize_drv_cond_Callback(hObject, eventdata, handles)
    value = str2double(get(hObject,'String'));
    if ~isnan(value) && value>0
        handles.MED_t_nucsize_drv_cond_v = value;
        guidata(hObject,handles);
    else
        value = handles.MED_t_nucsize_drv_cond_v;
    end
    set(hObject,'String',num2str(value));
    uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function MED_t_nucsize_drv_cond_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MED_t_nucsize_drv_cond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MED_small_smoothing_window_Callback(hObject, eventdata, handles)
    value = fix(str2double(get(hObject,'String')));
    if ~isnan(value) && value>10 && value<90
        handles.MED_small_smoothing_window_v = value;
        guidata(hObject,handles);
    else
        value = handles.MED_small_smoothing_window_v;
    end
    set(hObject,'String',num2str(value));
    uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function MED_small_smoothing_window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MED_small_smoothing_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MED_T_min_Callback(hObject, eventdata, handles)
    value = fix(str2double(get(hObject,'String')));
    if ~isnan(value) && value>1
        handles.MED_T_min_v = value;
        guidata(hObject,handles);
    else
        value = handles.MED_T_min_v;
    end
    set(hObject,'String',num2str(value));
    uiresume(handles.figure1);



% --- Executes during object creation, after setting all properties.
function MED_T_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MED_T_min (see GCBO)
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
        cur_handles.MI_LLeft = handles.MI_LLeft_v;
        cur_handles.MI_LRight= handles.MI_LRight_v;
        cur_handles.MED_t_nucsize_drv_uncond = handles.MED_t_nucsize_drv_uncond_v;
        cur_handles.MED_t_nucsize_drv_cond = handles.MED_t_nucsize_drv_cond_v;
        cur_handles.MED_dt_twin_max = handles.MED_dt_twin_max_v;
        cur_handles.MED_T_min = handles.MED_T_min_v;
        cur_handles.MED_big_smoothing_window = handles.MED_big_smoothing_window_v;
        cur_handles.MED_small_smoothing_window = handles.MED_small_smoothing_window_v;
        cur_handles.MED_t_nucsize_ratio = handles.MED_t_nucsize_ratio_v;
        cur_handles.MED_DT = handles.MED_DT_v;
        guidata(h,cur_handles);        
    end
    catch
    end
end
    fh = ancestor(hObject,'figure');     
    delete(fh);


% --- Executes on button press in pushbutton_Cancel.
function pushbutton_Cancel_Callback(hObject, eventdata, handles)
fh = ancestor(hObject,'figure');     
    delete(fh);



function MED_DT_Callback(hObject, eventdata, handles)
    value = fix(str2double(get(hObject,'String')));
    if ~isnan(value) && value>5 && value<90
        handles.MED_DT_v = value;
        guidata(hObject,handles);
    else
        value = handles.DT_v;
    end
    set(hObject,'String',num2str(value));    
    uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function MED_DT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MED_DT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MED_big_smoothing_window_Callback(hObject, eventdata, handles)
    value = fix(str2double(get(hObject,'String')));
    if ~isnan(value) && value>10 
        handles.MED_big_smoothing_window_v = value;
        guidata(hObject,handles);
    else
        value = handles.MED_big_smoothing_window_v;
    end
    set(hObject,'String',num2str(value));    
    uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function MED_big_smoothing_window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MED_big_smoothing_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MED_t_nucsize_ratio_Callback(hObject, eventdata, handles)
    value = str2double(get(hObject,'String'));
    if ~isnan(value) && value>1
        handles.MED_t_nucsize_ratio_v = value;
        guidata(hObject,handles);
    else
        value = handles.MED_t_nucsize_ratio_v;
    end
    set(hObject,'String',num2str(value));    
    uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function MED_t_nucsize_ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MED_t_nucsize_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MED_dt_twin_max_Callback(hObject, eventdata, handles)
    value = fix(str2double(get(hObject,'String')));
    if ~isnan(value) && value>=5
        handles.MED_dt_twin_max_v = value;
        guidata(hObject,handles);
    else
        value = handles.MED_dt_twin_max_v;
    end
    set(hObject,'String',num2str(value));    
    uiresume(handles.figure1);



% --- Executes during object creation, after setting all properties.
function MED_dt_twin_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MED_dt_twin_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MI_LLeft_Callback(hObject, eventdata, handles)
    value = fix(str2double(get(hObject,'String')));
    if ~isnan(value) && value>=5 && value<360
        handles.MI_LLeft_v = value;
        guidata(hObject,handles);
    else
        value = handles.MI_LLeft_v;
    end
    set(hObject,'String',num2str(value));    
    uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function MI_LLeft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MI_LLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function MI_LRight_Callback(hObject, eventdata, handles)
    value = fix(str2double(get(hObject,'String')));
    if ~isnan(value) && value>=5 && value<360
        handles.MI_LRight_v = value;
        guidata(hObject,handles);
    else
        value = handles.MI_LRight_v;
    end
    set(hObject,'String',num2str(value));    
    uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function MI_LRight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MI_LRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
