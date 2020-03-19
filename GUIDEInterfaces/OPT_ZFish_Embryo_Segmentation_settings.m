function varargout = OPT_ZFish_Embryo_Segmentation_settings(varargin)
% OPT_ZFISH_EMBRYO_SEGMENTATION_SETTINGS MATLAB code for OPT_ZFish_Embryo_Segmentation_settings.fig
%      OPT_ZFISH_EMBRYO_SEGMENTATION_SETTINGS, by itself, creates a new OPT_ZFISH_EMBRYO_SEGMENTATION_SETTINGS or raises the existing
%      singleton*.
%
%      H = OPT_ZFISH_EMBRYO_SEGMENTATION_SETTINGS returns the handle to a new OPT_ZFISH_EMBRYO_SEGMENTATION_SETTINGS or the handle to
%      the existing singleton*.
%
%      OPT_ZFISH_EMBRYO_SEGMENTATION_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPT_ZFISH_EMBRYO_SEGMENTATION_SETTINGS.M with the given input arguments.
%
%      OPT_ZFISH_EMBRYO_SEGMENTATION_SETTINGS('Property','Value',...) creates a new OPT_ZFISH_EMBRYO_SEGMENTATION_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OPT_ZFish_Embryo_Segmentation_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OPT_ZFish_Embryo_Segmentation_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OPT_ZFish_Embryo_Segmentation_settings

% Last Modified by GUIDE v2.5 19-Mar-2020 10:53:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OPT_ZFish_Embryo_Segmentation_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @OPT_ZFish_Embryo_Segmentation_settings_OutputFcn, ...
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


% --- Executes just before OPT_ZFish_Embryo_Segmentation_settings is made visible.
function OPT_ZFish_Embryo_Segmentation_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OPT_ZFish_Embryo_Segmentation_settings (see VARARGIN)

data_controller = varargin{1};
handles.data_controller = data_controller;

%26,25,28,27
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        handles.OPT_ZFish_Embryo_sgm_primary_scale = data_controller.OPT_ZFish_Embryo_sgm_primary_scale;
        handles.OPT_ZFish_Embryo_sgm_K21 = data_controller.OPT_ZFish_Embryo_sgm_K21;
        handles.OPT_ZFish_Embryo_sgm_K31 = data_controller.OPT_ZFish_Embryo_sgm_K31;
        handles.OPT_ZFish_Embryo_sgm_a1 = data_controller.OPT_ZFish_Embryo_sgm_a1; 
        handles.OPT_ZFish_Embryo_sgm_t = data_controller.OPT_ZFish_Embryo_sgm_t;
        handles.OPT_ZFish_Embryo_sgm_min_vol = data_controller.OPT_ZFish_Embryo_sgm_min_vol; 
        
        handles.saved_OPT_ZFish_Embryo_sgm_primary_scale = data_controller.OPT_ZFish_Embryo_sgm_primary_scale;
        handles.saved_OPT_ZFish_Embryo_sgm_K21 = data_controller.OPT_ZFish_Embryo_sgm_K21;
        handles.saved_OPT_ZFish_Embryo_sgm_K31 = data_controller.OPT_ZFish_Embryo_sgm_K31;
        handles.saved_OPT_ZFish_Embryo_sgm_a1 = data_controller.OPT_ZFish_Embryo_sgm_a1;
        handles.saved_OPT_ZFish_Embryo_sgm_t = data_controller.OPT_ZFish_Embryo_sgm_t;
        handles.saved_OPT_ZFish_Embryo_sgm_min_vol = data_controller.OPT_ZFish_Embryo_sgm_min_vol; 
                
        set(handles.edit26,'String',num2str(handles.OPT_ZFish_Embryo_sgm_primary_scale));
        set(handles.edit25,'String',num2str(handles.OPT_ZFish_Embryo_sgm_K21));
        set(handles.edit28,'String',num2str(handles.OPT_ZFish_Embryo_sgm_K31));
        set(handles.edit27,'String',num2str(handles.OPT_ZFish_Embryo_sgm_a1));
        set(handles.edit30,'String',num2str(handles.OPT_ZFish_Embryo_sgm_t));
        set(handles.edit29,'String',num2str(handles.OPT_ZFish_Embryo_sgm_min_vol));
                    
% Choose default command line output for OPT_ZFish_Embryo_Segmentation_settings
handles.output = hObject;

output_value = 1;

handles.output = output_value;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes OPT_ZFish_Embryo_Segmentation_settings wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = OPT_ZFish_Embryo_Segmentation_settings_OutputFcn(hObject, eventdata, handles) 
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
         data_controller.OPT_ZFish_Embryo_sgm_primary_scale = handles.OPT_ZFish_Embryo_sgm_primary_scale;
         data_controller.OPT_ZFish_Embryo_sgm_K21 = handles.OPT_ZFish_Embryo_sgm_K21;
         data_controller.OPT_ZFish_Embryo_sgm_K31 = handles.OPT_ZFish_Embryo_sgm_K31;
         data_controller.OPT_ZFish_Embryo_sgm_a1 = handles.OPT_ZFish_Embryo_sgm_a1;
         data_controller.OPT_ZFish_Embryo_sgm_t = handles.OPT_ZFish_Embryo_sgm_t;
         data_controller.OPT_ZFish_Embryo_sgm_min_vol = handles.OPT_ZFish_Embryo_sgm_min_vol; 

    fh = ancestor(hObject,'figure');     
    delete(fh);

% Cancel    
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    data_controller = handles.data_controller;
         data_controller.OPT_ZFish_Embryo_sgm_primary_scale = handles.saved_OPT_ZFish_Embryo_sgm_primary_scale;
         data_controller.OPT_ZFish_Embryo_sgm_K21 = handles.saved_OPT_ZFish_Embryo_sgm_K21;
         data_controller.OPT_ZFish_Embryo_sgm_K31 = handles.saved_OPT_ZFish_Embryo_sgm_K31;
         data_controller.OPT_ZFish_Embryo_sgm_a1 = handles.saved_OPT_ZFish_Embryo_sgm_a1;
         data_controller.OPT_ZFish_Embryo_sgm_t = handles.saved_OPT_ZFish_Embryo_sgm_t;
         data_controller.OPT_ZFish_Embryo_sgm_min_vol = handles.saved_OPT_ZFish_Embryo_sgm_min_vol;          
         
    fh = ancestor(hObject,'figure');     
    delete(fh);

% Set default
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data_controller = handles.data_controller;

         data_controller.OPT_ZFish_Embryo_sgm_primary_scale = 88;
         data_controller.OPT_ZFish_Embryo_sgm_K21 = 2;
         data_controller.OPT_ZFish_Embryo_sgm_K31 = 4;
         data_controller.OPT_ZFish_Embryo_sgm_a1 = 0.5;
         data_controller.OPT_ZFish_Embryo_sgm_t = 0.4;
         data_controller.OPT_ZFish_Embryo_sgm_min_vol = 681472000;
         
         handles.OPT_ZFish_Embryo_sgm_primary_scale = data_controller.OPT_ZFish_Embryo_sgm_primary_scale;
         handles.OPT_ZFish_Embryo_sgm_K21 = data_controller.OPT_ZFish_Embryo_sgm_K21;
         handles.OPT_ZFish_Embryo_sgm_K31 = data_controller.OPT_ZFish_Embryo_sgm_K31;
         handles.OPT_ZFish_Embryo_sgm_a1 = data_controller.OPT_ZFish_Embryo_sgm_a1; 
         handles.OPT_ZFish_Embryo_sgm_t = data_controller.OPT_ZFish_Embryo_sgm_t;
         handles.OPT_ZFish_Embryo_sgm_min_vol = data_controller.OPT_ZFish_Embryo_sgm_min_vol; 
         
         set(handles.edit26,'String',num2str(data_controller.OPT_ZFish_Embryo_sgm_primary_scale));
         set(handles.edit25,'String',num2str(data_controller.OPT_ZFish_Embryo_sgm_K21));
         set(handles.edit28,'String',num2str(data_controller.OPT_ZFish_Embryo_sgm_K31));
         set(handles.edit27,'String',num2str(data_controller.OPT_ZFish_Embryo_sgm_a1));
         set(handles.edit30,'String',num2str(data_controller.OPT_ZFish_Embryo_sgm_t));
         set(handles.edit29,'String',num2str(data_controller.OPT_ZFish_Embryo_sgm_min_vol));         

 guidata(hObject, handles);
% UIWAIT makes OPT_ZFish_Embryo_Segmentation_settings wait for user response (see UIRESUME)
uiwait(handles.figure1);       
        
%     fh = ancestor(hObject,'figure');     
%     delete(fh);

function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) || value<0 || value>1
    handles.OPT_ZFish_Embryo_sgm_primary_scale = value;
    guidata(hObject,handles);
else
    value = handles.OPT_ZFish_Embryo_sgm_primary_scale;
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
if ~isnan(value) || value<0 || value>1
    handles.OPT_ZFish_Embryo_sgm_K21 = value;
    guidata(hObject,handles);
else
    value = handles.OPT_ZFish_Embryo_sgm_K21;
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
         data_controller.OPT_ZFish_Embryo_sgm_primary_scale = handles.OPT_ZFish_Embryo_sgm_primary_scale;
         data_controller.OPT_ZFish_Embryo_sgm_K21 = handles.OPT_ZFish_Embryo_sgm_K21;
         data_controller.OPT_ZFish_Embryo_sgm_K31 = handles.OPT_ZFish_Embryo_sgm_K31;
         data_controller.OPT_ZFish_Embryo_sgm_a1 = handles.OPT_ZFish_Embryo_sgm_a1;  
         data_controller.OPT_ZFish_Embryo_sgm_t = handles.OPT_ZFish_Embryo_sgm_t;
         data_controller.OPT_ZFish_Embryo_sgm_min_vol = handles.OPT_ZFish_Embryo_sgm_min_vol;        
         %
         send_adjustment_to_Icy = true;
         data_controller.do_OPT_ZFish_Embryo_Segmentation(send_adjustment_to_Icy);


function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) || value<0 || value>1
    handles.OPT_ZFish_Embryo_sgm_a1 = value;
    guidata(hObject,handles);
else
    value = handles.OPT_ZFish_Embryo_sgm_a1;
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
if ~isnan(value) || value<1
    handles.OPT_ZFish_Embryo_sgm_K31 = value;
    guidata(hObject,handles);
else
    value = handles.OPT_ZFish_Embryo_sgm_K31;
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
value = str2double(get(hObject,'String'));
if ~isnan(value) || value<0
    handles.OPT_ZFish_Embryo_sgm_min_vol = value;
    guidata(hObject,handles);
else
    value = handles.OPT_ZFish_Embryo_sgm_min_vol;
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
value = str2double(get(hObject,'String'));
if ~isnan(value) || value<0
    handles.OPT_ZFish_Embryo_sgm_t = value;
    guidata(hObject,handles);
else
    value = handles.OPT_ZFish_Embryo_sgm_t;
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
