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

% Last Modified by GUIDE v2.5 18-May-2019 22:05:29

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

dc = varargin{1};
if isempty(dc.imgdata), 
    errordlg('no data loaded - can not continue'), 
    fh = ancestor(hObject,'figure');     
    delete(fh);
    return, 
end
handles.dc = dc;
        % 
handles.network =                       dc.AI_Powered_2D_SMLM_Reconstruction_network; % path to file
handles.upscale_factor =                dc.AI_Powered_2D_SMLM_Reconstruction_upscale_factor;
handles.vicinity_half_width =           dc.AI_Powered_2D_SMLM_Reconstruction_vicinity_half_width;
handles.pixel_size =                    dc.AI_Powered_2D_SMLM_Reconstruction_pixel_size; % nm/pixel
handles.extraction_scale =              dc.AI_Powered_2D_SMLM_Reconstruction_extraction_scale; %non-super res pixels            
handles.extraction_scale_ratio =        dc.AI_Powered_2D_SMLM_Reconstruction_extraction_scale_ratio; % ditto
handles.extraction_threshold =          dc.AI_Powered_2D_SMLM_Reconstruction_extraction_threshold; % std factor?
handles.extraction_method =             dc.AI_Powered_2D_SMLM_Reconstruction_extraction_method;
handles.max_distance_to_spurious_pixl = dc.AI_Powered_2D_SMLM_Reconstruction_max_distance_to_spurious_pixl; % nm
handles.time_dependent_block_size =     dc.AI_Powered_2D_SMLM_Reconstruction_time_dependent_block_size; % frames            
handles.image_formation_method =        dc.AI_Powered_2D_SMLM_Reconstruction_image_formation_method;
handles.image_formation_scale =         dc.AI_Powered_2D_SMLM_Reconstruction_image_formation_scale;
handles.NA =                            dc.AI_Powered_2D_SMLM_Reconstruction_NA;
handles.wavelength =                    dc.AI_Powered_2D_SMLM_Reconstruction_wavelength;

% back up initial values for "cancel"
handles.ini_network =                       dc.AI_Powered_2D_SMLM_Reconstruction_network; % path to file
handles.ini_upscale_factor =                dc.AI_Powered_2D_SMLM_Reconstruction_upscale_factor;
handles.ini_vicinity_half_width =           dc.AI_Powered_2D_SMLM_Reconstruction_vicinity_half_width;
handles.ini_pixel_size =                    dc.AI_Powered_2D_SMLM_Reconstruction_pixel_size; % nm/pixel
handles.ini_extraction_scale =              dc.AI_Powered_2D_SMLM_Reconstruction_extraction_scale; %non-super res pixels            
handles.ini_extraction_scale_ratio =        dc.AI_Powered_2D_SMLM_Reconstruction_extraction_scale_ratio; % ditto
handles.ini_extraction_threshold =          dc.AI_Powered_2D_SMLM_Reconstruction_extraction_threshold; % std factor?
handles.ini_extraction_method =             dc.AI_Powered_2D_SMLM_Reconstruction_extraction_method;
handles.ini_max_distance_to_spurious_pixl = dc.AI_Powered_2D_SMLM_Reconstruction_max_distance_to_spurious_pixl; % nm
handles.ini_time_dependent_block_size =     dc.AI_Powered_2D_SMLM_Reconstruction_time_dependent_block_size; % frames            
handles.ini_image_formation_method =        dc.AI_Powered_2D_SMLM_Reconstruction_image_formation_method;
handles.ini_image_formation_scale =         dc.AI_Powered_2D_SMLM_Reconstruction_image_formation_scale;
handles.ini_NA =                            dc.AI_Powered_2D_SMLM_Reconstruction_NA;
handles.ini_wavelength =                    dc.AI_Powered_2D_SMLM_Reconstruction_wavelength;
% back up initial values for "cancel"

set(handles.image_formation_method_popupmenu,'String',dc.AI_Powered_2D_SMLM_Reconstruction_image_formation_methods);
index = find(strcmp(dc.AI_Powered_2D_SMLM_Reconstruction_image_formation_methods,dc.AI_Powered_2D_SMLM_Reconstruction_image_formation_method));
set(handles.image_formation_method_popupmenu,'Value',index);

set(handles.extraction_method_popupmenu,'String',dc.AI_Powered_2D_SMLM_Reconstruction_extraction_methods); 
index = find(strcmp(dc.AI_Powered_2D_SMLM_Reconstruction_extraction_methods,dc.AI_Powered_2D_SMLM_Reconstruction_extraction_method));
set(handles.extraction_method_popupmenu,'Value',index);
                      
set(handles.network_edit,'String',handles.network);
set(handles.upscale_factor_edit,'String',num2str(handles.upscale_factor));
set(handles.vicinity_half_width_edit,'String',num2str(handles.vicinity_half_width));
set(handles.pix_size_edit,'String',num2str(handles.pixel_size)); 
set(handles.extraction_scale_edit,'String',num2str(handles.extraction_scale));
set(handles.scale_ratio_edit,'String',num2str(handles.extraction_scale_ratio));
set(handles.extraction_threshold_edit,'String',num2str(handles.extraction_threshold));
set(handles.spurious_emitter_distance_edit,'String',num2str(handles.max_distance_to_spurious_pixl));
set(handles.frame_merging_block_size_edit,'String',num2str(handles.time_dependent_block_size));
set(handles.super_res_image_formation_scale_edit,'String',num2str(handles.image_formation_scale)); 
set(handles.NA_edit,'String',num2str(handles.NA));
set(handles.wavelength_edit,'String',num2str(handles.wavelength));
      
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
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 1 && value <= 1000
    handles.pixel_size = value;
    guidata(hObject,handles);
else
    value = handles.pixel_size;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

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
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 2 && value <= 10
    handles.upscale_factor = value;
    guidata(hObject,handles);
else
    value = handles.upscale_factor;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

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



function NA_edit_Callback(hObject, eventdata, handles)
% hObject    handle to NA_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NA_edit as text
%        str2double(get(hObject,'String')) returns contents of NA_edit as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= .1 && value <= 1.9
    handles.NA = value;
    guidata(hObject,handles);
else
    value = handles.NA;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function NA_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NA_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vicinity_half_width_edit_Callback(hObject, eventdata, handles)
% hObject    handle to vicinity_half_width_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vicinity_half_width_edit as text
%        str2double(get(hObject,'String')) returns contents of vicinity_half_width_edit as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 4 && value <= 20
    handles.vicinity_half_width = 2*ceil(value/2)+1;
    guidata(hObject,handles);
else
    value = handles.vicinity_half_width;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function vicinity_half_width_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vicinity_half_width_edit (see GCBO)
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
%  NO VALIDITY CHECK
set(handles.network_edit,'String',get(hObject,'String')); % NO VALIDITY CHECK
%  NO VALIDITY CHECK

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
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= .5 && value <= 15
    handles.extraction_scale = value;
    guidata(hObject,handles);
else
    value = handles.extraction_scale;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

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
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 1.2 && value <= 10
    handles.extraction_scale_ratio = value;
    guidata(hObject,handles);
else
    value = handles.extraction_scale_ratio;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

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



function extraction_method_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to extraction_method_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of extraction_method_popupmenu as text
%        str2double(get(hObject,'String')) returns contents of extraction_method_popupmenu as a double
index = get(hObject,'Value');
str = get(hObject,'String');
handles.extraction_method = str(index);

% --- Executes during object creation, after setting all properties.
function extraction_method_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extraction_method_popupmenu (see GCBO)
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
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 0 
    handles.extraction_threshold = value;
    guidata(hObject,handles);
else
    value = handles.extraction_threshold;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

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
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 10 && value <= 600
    handles.time_dependent_block_size = value;
    guidata(hObject,handles);
else
    value = handles.time_dependent_block_size;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

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



function image_formation_method_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to image_formation_method_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of image_formation_method_popupmenu as text
%        str2double(get(hObject,'String')) returns contents of image_formation_method_popupmenu as a double
index = get(hObject,'Value');
str = get(hObject,'String');
handles.image_formation_method = str(index);

% --- Executes during object creation, after setting all properties.
function image_formation_method_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_formation_method_popupmenu (see GCBO)
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
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 0 && value <= 50
    handles.image_formation_scale = value;
    guidata(hObject,handles);
else
    value = handles.image_formation_scale;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

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
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 1 && value <= 5000
    handles.max_distance_to_spurious_pixl = value;
    guidata(hObject,handles);
else
    value = handles.max_distance_to_spurious_pixl;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

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

function update_handles(handles)
dc = handles.dc;
dc.AI_Powered_2D_SMLM_Reconstruction_network = handles.network;
dc.AI_Powered_2D_SMLM_Reconstruction_upscale_factor = handles.upscale_factor;
dc.AI_Powered_2D_SMLM_Reconstruction_vicinity_half_width = handles.vicinity_half_width;
dc.AI_Powered_2D_SMLM_Reconstruction_pixel_size = handles.pixel_size;
dc.AI_Powered_2D_SMLM_Reconstruction_extraction_scale = handles.extraction_scale;  
dc.AI_Powered_2D_SMLM_Reconstruction_extraction_scale_ratio = handles.extraction_scale_ratio;
dc.AI_Powered_2D_SMLM_Reconstruction_extraction_threshold = handles.extraction_threshold;
dc.AI_Powered_2D_SMLM_Reconstruction_extraction_method = handles.extraction_method;
dc.AI_Powered_2D_SMLM_Reconstruction_max_distance_to_spurious_pixl = handles.max_distance_to_spurious_pixl;
dc.AI_Powered_2D_SMLM_Reconstruction_time_dependent_block_size = handles.time_dependent_block_size;
dc.AI_Powered_2D_SMLM_Reconstruction_image_formation_method = handles.image_formation_method;
dc.AI_Powered_2D_SMLM_Reconstruction_image_formation_scale = handles.image_formation_scale;
dc.AI_Powered_2D_SMLM_Reconstruction_NA = handles.NA;
dc.AI_Powered_2D_SMLM_Reconstruction_wavelength = handles.wavelength;

% --- Executes on button press in update_pushbutton.
function update_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to update_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_handles(handles);

% --- Executes on button press in update_and_close_pushbutton.
function update_and_close_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to update_and_close_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_handles(handles);
    fh = ancestor(hObject,'figure');     
    delete(fh);

% --- Executes on button press in cancel_pushbutton.
function cancel_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%needs to restore initial values
dc = handles.dc;
dc.AI_Powered_2D_SMLM_Reconstruction_network = handles.ini_network;
dc.AI_Powered_2D_SMLM_Reconstruction_upscale_factor = handles.ini_upscale_factor;
dc.AI_Powered_2D_SMLM_Reconstruction_vicinity_half_width = handles.ini_vicinity_half_width;
dc.AI_Powered_2D_SMLM_Reconstruction_pixel_size = handles.ini_pixel_size;
dc.AI_Powered_2D_SMLM_Reconstruction_extraction_scale = handles.ini_extraction_scale;  
dc.AI_Powered_2D_SMLM_Reconstruction_extraction_scale_ratio = handles.ini_extraction_scale_ratio;
dc.AI_Powered_2D_SMLM_Reconstruction_extraction_threshold = handles.ini_extraction_threshold;
dc.AI_Powered_2D_SMLM_Reconstruction_extraction_method = handles.ini_extraction_method;
dc.AI_Powered_2D_SMLM_Reconstruction_max_distance_to_spurious_pixl = handles.ini_max_distance_to_spurious_pixl;
dc.AI_Powered_2D_SMLM_Reconstruction_time_dependent_block_size = handles.ini_time_dependent_block_size;
dc.AI_Powered_2D_SMLM_Reconstruction_image_formation_method = handles.ini_image_formation_method;
dc.AI_Powered_2D_SMLM_Reconstruction_image_formation_scale = handles.ini_image_formation_scale;
dc.AI_Powered_2D_SMLM_Reconstruction_NA = handles.ini_NA;
dc.AI_Powered_2D_SMLM_Reconstruction_wavelength = handles.ini_wavelength;

    fh = ancestor(hObject,'figure');     
    delete(fh);

% --- Executes on button press in network_path_pushbutton.
function network_path_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to network_path_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function wavelength_edit_Callback(hObject, eventdata, handles)
% hObject    handle to wavelength_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wavelength_edit as text
%        str2double(get(hObject,'String')) returns contents of wavelength_edit as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 200 && value <= 800
    handles.wavelength = value;
    guidata(hObject,handles);
else
    value = handles.wavelength;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function wavelength_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wavelength_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
