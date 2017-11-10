function varargout = ALYtools_General_Settings(varargin)
% ALYTOOLS_GENERAL_SETTINGS MATLAB code for ALYtools_General_Settings.fig
%      ALYTOOLS_GENERAL_SETTINGS, by itself, creates a new ALYTOOLS_GENERAL_SETTINGS or raises the existing
%      singleton*.
%
%      H = ALYTOOLS_GENERAL_SETTINGS returns the handle to a new ALYTOOLS_GENERAL_SETTINGS or the handle to
%      the existing singleton*.
%
%      ALYTOOLS_GENERAL_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ALYTOOLS_GENERAL_SETTINGS.M with the given input arguments.
%
%      ALYTOOLS_GENERAL_SETTINGS('Property','Value',...) creates a new ALYTOOLS_GENERAL_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ALYtools_General_Settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ALYtools_General_Settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ALYtools_General_Settings

% Last Modified by GUIDE v2.5 08-Nov-2017 17:34:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ALYtools_General_Settings_OpeningFcn, ...
                   'gui_OutputFcn',  @ALYtools_General_Settings_OutputFcn, ...
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


% --- Executes just before ALYtools_General_Settings is made visible.
function ALYtools_General_Settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ALYtools_General_Settings (see VARARGIN)

data_controller = varargin{1};
handles.data_controller = data_controller;
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.send_analysis_output_to_Icy = data_controller.send_analysis_output_to_Icy;
        handles.save_analysis_output_as_OMEtiff = data_controller.save_analysis_output_as_OMEtiff;
        handles.save_analysis_output_as_xls = data_controller.save_analysis_output_as_xls;
        
        set(handles.checkbox1,'Value',handles.send_analysis_output_to_Icy);
        set(handles.checkbox2,'Value',handles.save_analysis_output_as_OMEtiff);
        set(handles.checkbox3,'Value',handles.save_analysis_output_as_xls);
         
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ALYtools_General_Settings wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ALYtools_General_Settings_OutputFcn(hObject, eventdata, handles) 
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% granules
         data_controller.send_analysis_output_to_Icy = handles.send_analysis_output_to_Icy;
         data_controller.save_analysis_output_as_OMEtiff = handles.save_analysis_output_as_OMEtiff;
         data_controller.save_analysis_output_as_xls = handles.save_analysis_output_as_xls;         
         
    fh = ancestor(hObject,'figure');     
    delete(fh);

% Cancel    
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    %data_controller = handles.data_controller;
    fh = ancestor(hObject,'figure');     
    delete(fh);


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
value = get(hObject,'Value');
if ~isnan(value)
    handles.send_analysis_output_to_Icy = value;
    guidata(hObject,handles);
else
    value = handles.send_analysis_output_to_Icy;
    set(hObject,'Value',value);
end
uiresume(handles.figure1);


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
value = get(hObject,'Value');
if ~isnan(value)
    handles.save_analysis_output_as_OMEtiff = value;
    guidata(hObject,handles);
else
    value = handles.save_analysis_output_as_OMEtiff;
    set(hObject,'Value',value);
end
uiresume(handles.figure1);


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
value = get(hObject,'Value');
if ~isnan(value)
    handles.save_analysis_output_as_xls = value;
    guidata(hObject,handles);
else
    value = handles.save_analysis_output_as_xls;
    set(hObject,'Value',value);
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
