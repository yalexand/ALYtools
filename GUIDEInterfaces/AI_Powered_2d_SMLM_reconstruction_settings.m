function varargout = AI_Powered_2d_SMLM_reconstruction_settings(varargin)
% AI_POWERED_2D_SMLM_RECONSTRUCTION_SETTINGS MATLAB code for AI_Powered_2d_SMLM_reconstruction_settings.fig
%      AI_POWERED_2D_SMLM_RECONSTRUCTION_SETTINGS, by itself, creates a new AI_POWERED_2D_SMLM_RECONSTRUCTION_SETTINGS or raises the existing
%      singleton*.
%
%      H = AI_POWERED_2D_SMLM_RECONSTRUCTION_SETTINGS returns the handle to a new AI_POWERED_2D_SMLM_RECONSTRUCTION_SETTINGS or the handle to
%      the existing singleton*.
%
%      AI_POWERED_2D_SMLM_RECONSTRUCTION_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AI_POWERED_2D_SMLM_RECONSTRUCTION_SETTINGS.M with the given input arguments.
%
%      AI_POWERED_2D_SMLM_RECONSTRUCTION_SETTINGS('Property','Value',...) creates a new AI_POWERED_2D_SMLM_RECONSTRUCTION_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AI_Powered_2d_SMLM_reconstruction_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AI_Powered_2d_SMLM_reconstruction_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AI_Powered_2d_SMLM_reconstruction_settings

% Last Modified by GUIDE v2.5 17-May-2019 15:37:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AI_Powered_2d_SMLM_reconstruction_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @AI_Powered_2d_SMLM_reconstruction_settings_OutputFcn, ...
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


% --- Executes just before AI_Powered_2d_SMLM_reconstruction_settings is made visible.
function AI_Powered_2d_SMLM_reconstruction_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AI_Powered_2d_SMLM_reconstruction_settings (see VARARGIN)

% Choose default command line output for AI_Powered_2d_SMLM_reconstruction_settings
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AI_Powered_2d_SMLM_reconstruction_settings wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AI_Powered_2d_SMLM_reconstruction_settings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function pix_size_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pix_size_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pix_size_edit as text
%        str2double(get(hObject,'String')) returns contents of pix_size_edit as a double


% --- Executes during object creation, after setting all properties.
function pix_size_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pix_size_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function upscale_factor_edit_Callback(hObject, eventdata, handles)
% hObject    handle to upscale_factor_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upscale_factor_edit as text
%        str2double(get(hObject,'String')) returns contents of upscale_factor_edit as a double


% --- Executes during object creation, after setting all properties.
function upscale_factor_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upscale_factor_edit (see GCBO)
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



function network_edit_Callback(hObject, eventdata, handles)
% hObject    handle to network_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of network_edit as text
%        str2double(get(hObject,'String')) returns contents of network_edit as a double


% --- Executes during object creation, after setting all properties.
function network_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to network_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function extraction_scale_edit_Callback(hObject, eventdata, handles)
% hObject    handle to extraction_scale_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of extraction_scale_edit as text
%        str2double(get(hObject,'String')) returns contents of extraction_scale_edit as a double


% --- Executes during object creation, after setting all properties.
function extraction_scale_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extraction_scale_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function scale_ratio_edit_Callback(hObject, eventdata, handles)
% hObject    handle to scale_ratio_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scale_ratio_edit as text
%        str2double(get(hObject,'String')) returns contents of scale_ratio_edit as a double


% --- Executes during object creation, after setting all properties.
function scale_ratio_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scale_ratio_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_extraction_method_Callback(hObject, eventdata, handles)
% hObject    handle to edit_extraction_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_extraction_method as text
%        str2double(get(hObject,'String')) returns contents of edit_extraction_method as a double


% --- Executes during object creation, after setting all properties.
function edit_extraction_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_extraction_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function extraction_threshold_edit_Callback(hObject, eventdata, handles)
% hObject    handle to extraction_threshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of extraction_threshold_edit as text
%        str2double(get(hObject,'String')) returns contents of extraction_threshold_edit as a double


% --- Executes during object creation, after setting all properties.
function extraction_threshold_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extraction_threshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_merging_block_size_edit_Callback(hObject, eventdata, handles)
% hObject    handle to frame_merging_block_size_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_merging_block_size_edit as text
%        str2double(get(hObject,'String')) returns contents of frame_merging_block_size_edit as a double


% --- Executes during object creation, after setting all properties.
function frame_merging_block_size_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_merging_block_size_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function image_formation_method_edit_Callback(hObject, eventdata, handles)
% hObject    handle to image_formation_method_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of image_formation_method_edit as text
%        str2double(get(hObject,'String')) returns contents of image_formation_method_edit as a double


% --- Executes during object creation, after setting all properties.
function image_formation_method_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_formation_method_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function super_res_image_formation_scale_edit_Callback(hObject, eventdata, handles)
% hObject    handle to super_res_image_formation_scale_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of super_res_image_formation_scale_edit as text
%        str2double(get(hObject,'String')) returns contents of super_res_image_formation_scale_edit as a double


% --- Executes during object creation, after setting all properties.
function super_res_image_formation_scale_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to super_res_image_formation_scale_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function spurious_emitter_distance_edit_Callback(hObject, eventdata, handles)
% hObject    handle to spurious_emitter_distance_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spurious_emitter_distance_edit as text
%        str2double(get(hObject,'String')) returns contents of spurious_emitter_distance_edit as a double


% --- Executes during object creation, after setting all properties.
function spurious_emitter_distance_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spurious_emitter_distance_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in update_pushbutton.
function update_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to update_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in update_and_close_pushbutton.
function update_and_close_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to update_and_close_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cancel_pushbutton.
function cancel_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in network_path_pushbutton.
function network_path_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to network_path_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
