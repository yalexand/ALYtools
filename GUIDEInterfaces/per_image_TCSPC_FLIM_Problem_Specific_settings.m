function varargout = per_image_TCSPC_FLIM_Problem_Specific_settings(varargin)
% PER_IMAGE_TCSPC_FLIM_PROBLEM_SPECIFIC_SETTINGS MATLAB code for per_image_TCSPC_FLIM_Problem_Specific_settings.fig
%      PER_IMAGE_TCSPC_FLIM_PROBLEM_SPECIFIC_SETTINGS, by itself, creates a new PER_IMAGE_TCSPC_FLIM_PROBLEM_SPECIFIC_SETTINGS or raises the existing
%      singleton*.
%
%      H = PER_IMAGE_TCSPC_FLIM_PROBLEM_SPECIFIC_SETTINGS returns the handle to a new PER_IMAGE_TCSPC_FLIM_PROBLEM_SPECIFIC_SETTINGS or the handle to
%      the existing singleton*.
%
%      PER_IMAGE_TCSPC_FLIM_PROBLEM_SPECIFIC_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PER_IMAGE_TCSPC_FLIM_PROBLEM_SPECIFIC_SETTINGS.M with the given input arguments.
%
%      PER_IMAGE_TCSPC_FLIM_PROBLEM_SPECIFIC_SETTINGS('Property','Value',...) creates a new PER_IMAGE_TCSPC_FLIM_PROBLEM_SPECIFIC_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before per_image_TCSPC_FLIM_Problem_Specific_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to per_image_TCSPC_FLIM_Problem_Specific_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help per_image_TCSPC_FLIM_Problem_Specific_settings

% Last Modified by GUIDE v2.5 05-Jun-2018 12:04:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @per_image_TCSPC_FLIM_Problem_Specific_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @per_image_TCSPC_FLIM_Problem_Specific_settings_OutputFcn, ...
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




% --- Executes just before per_image_TCSPC_FLIM_Problem_Specific_settings is made visible.
function per_image_TCSPC_FLIM_Problem_Specific_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to per_image_TCSPC_FLIM_Problem_Specific_settings (see VARARGIN)

dc = varargin{1};
if isempty(dc.M_imgdata), 
    errordlg('no data loaded - can not continue'), 
    fh = ancestor(hObject,'figure');     
    delete(fh);
    return, 
end;

if strcmp(dc.problem,'per_image_TCSPC_FLIM_PHASOR')
    set(handles.text29,'Visible','off');
    set(handles.fitting_model,'Visible','off');
    set(handles.fixed_tauD_checkbox,'Visible','off');
    set(handles.fixed_tauD,'Visible','off');    
end

            handles.dc = dc;         
%             per_image_TCSPC_FLIM_t = [];
%             per_image_TCSPC_FLIM_irf = [];
            handles.per_image_TCSPC_FLIM_irf_shift = dc.per_image_TCSPC_FLIM_irf_shift;
            handles.per_image_TCSPC_FLIM_irf_background = dc.per_image_TCSPC_FLIM_irf_background;
            handles.per_image_TCSPC_FLIM_rep_rate = dc.per_image_TCSPC_FLIM_rep_rate;
            handles.per_image_TCSPC_FLIM_irf_filename = dc.per_image_TCSPC_FLIM_irf_filename;
            %
            handles.per_image_TCSPC_FLIM_fit_model = dc.per_image_TCSPC_FLIM_fit_model;
            handles.per_image_TCSPC_FLIM_Tmin = dc.per_image_TCSPC_FLIM_Tmin;
            handles.per_image_TCSPC_FLIM_Tmax = dc.per_image_TCSPC_FLIM_Tmax;
            handles.per_image_TCSPC_FLIM_background_value = dc.per_image_TCSPC_FLIM_background_value;
            handles.per_image_TCSPC_FLIM_saturation_value = dc.per_image_TCSPC_FLIM_saturation_value;
                        
            handles.per_image_TCSPC_FLIM_irf_filename = dc.per_image_TCSPC_FLIM_irf_filename;
            handles.per_image_TCSPC_FLIM_tvb_filename = dc.per_image_TCSPC_FLIM_tvb_filename;
            handles.per_image_TCSPC_FLIM_tvb_scaling = dc.per_image_TCSPC_FLIM_tvb_scaling;
            
            handles.per_image_TCSPC_FLIM_fixed_tauD = dc.per_image_TCSPC_FLIM_fixed_tauD;
            handles.per_image_TCSPC_FLIM_conv_irf_pp_69_70 = dc.per_image_TCSPC_FLIM_conv_irf_pp_69_70;

            handles.per_image_TCSPC_FLIM_Ref_lifetime = dc.per_image_TCSPC_FLIM_Ref_lifetime;
            handles.per_image_TCSPC_FLIM_averaging_sigma = dc.per_image_TCSPC_FLIM_averaging_sigma;
                        
            set(handles.irf_shift,'String',num2str(handles.per_image_TCSPC_FLIM_irf_shift));
            set(handles.IRF_background,'String',num2str(handles.per_image_TCSPC_FLIM_irf_background)); 
            set(handles.fitting_model,'String',dc.global_fit_models');
            %
            %needs to set up model choice..
            handles.per_image_TCSPC_FLIM_fit_model = dc.per_image_TCSPC_FLIM_fit_model;
            index = find(strcmp(dc.global_fit_models,dc.per_image_TCSPC_FLIM_fit_model));                        
            set(handles.fitting_model,'Value',index);
            %
            set(handles.rep_rate,'String',num2str(handles.per_image_TCSPC_FLIM_rep_rate));
            set(handles.dead_zone_left,'String',num2str(handles.per_image_TCSPC_FLIM_Tmin));
            set(handles.dead_zone_right,'String',num2str(handles.per_image_TCSPC_FLIM_Tmax));
            set(handles.background_value,'String',num2str(handles.per_image_TCSPC_FLIM_background_value));
            set(handles.tvb_scaling,'String',num2str(handles.per_image_TCSPC_FLIM_tvb_scaling));
            set(handles.saturation_value,'String',num2str(handles.per_image_TCSPC_FLIM_saturation_value));
            set(handles.IRF_reference_lifetime,'String',num2str(handles.per_image_TCSPC_FLIM_Ref_lifetime));
            set(handles.averaging_sigma,'String',num2str(handles.per_image_TCSPC_FLIM_averaging_sigma));
            %
            if ~isempty(handles.per_image_TCSPC_FLIM_irf_filename)
                [~,IRF_NAME,IRF_EXT] = fileparts(handles.per_image_TCSPC_FLIM_irf_filename);
                set(handles.irf_name_text,'String',[IRF_NAME IRF_EXT]);
                guidata(hObject, handles);                
                uiresume(handles.figure1);                
            end
            %
            if ~isempty(handles.per_image_TCSPC_FLIM_tvb_filename)
                [~,TVB_NAME,TVB_EXT] = fileparts(handles.per_image_TCSPC_FLIM_tvb_filename);                
                set(handles.tvb_name_text,'String',[TVB_NAME TVB_EXT]);
                guidata(hObject, handles);                
                uiresume(handles.figure1);
            end
            
    if 0 == handles.per_image_TCSPC_FLIM_fixed_tauD
        set(handles.fixed_tauD,'Enable','off');
        set(handles.fixed_tauD,'String','...');
        set(handles.fixed_tauD_checkbox,'Value',0);
    else
        set(handles.fixed_tauD,'Enable','on');
        set(handles.fixed_tauD,'String',num2str(handles.per_image_TCSPC_FLIM_fixed_tauD));
        set(handles.fixed_tauD_checkbox,'Value',1);        
    end
            
    if 0 == handles.per_image_TCSPC_FLIM_conv_irf_pp_69_70
        set(handles.conv_pp_69_70,'Value',0);
    else
        set(handles.conv_pp_69_70,'Value',1);        
    end
                           
% Choose default command line output for per_image_TCSPC_FLIM_Problem_Specific_settings
handles.output = hObject;

output_value = 1;

handles.output = output_value;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes per_image_TCSPC_FLIM_Problem_Specific_settings wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = per_image_TCSPC_FLIM_Problem_Specific_settings_OutputFcn(hObject, eventdata, handles) 
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

% function edit3_Callback(hObject, eventdata, handles)
% % hObject    handle to edit3 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edit3 as text
% %        str2double(get(hObject,'String')) returns contents of edit3 as a double
% data_controller = handles.data_controller;
% [~,~,~,sC,~]=size(data_controller.imgdata);
% value = fix(str2double(get(hObject,'String')));
% if ~isnan(value) && value >= 1 && value <= sC
%     handles.NC_chNuc = fix(value);
%     guidata(hObject,handles);
% else
%     value = handles.NC_chNuc;
%     set(hObject,'String',num2str(value));
% end
% uiresume(handles.figure1);
% 
% 
% % --- Executes during object creation, after setting all properties.
% function edit3_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit3 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% OK
% --- Executes on button press in OK.
function OK_Callback(hObject, eventdata, handles)
% hObject    handle to OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dc = handles.dc;
        
            dc.per_image_TCSPC_FLIM_irf_shift = handles.per_image_TCSPC_FLIM_irf_shift;
            dc.per_image_TCSPC_FLIM_irf_background = handles.per_image_TCSPC_FLIM_irf_background;
            dc.per_image_TCSPC_FLIM_rep_rate = handles.per_image_TCSPC_FLIM_rep_rate;
            %
            dc.per_image_TCSPC_FLIM_fit_model = handles.per_image_TCSPC_FLIM_fit_model;
            dc.per_image_TCSPC_FLIM_Tmin = handles.per_image_TCSPC_FLIM_Tmin;
            dc.per_image_TCSPC_FLIM_Tmax = handles.per_image_TCSPC_FLIM_Tmax;
            dc.per_image_TCSPC_FLIM_background_value = handles.per_image_TCSPC_FLIM_background_value;
            dc.per_image_TCSPC_FLIM_saturation_value = handles.per_image_TCSPC_FLIM_saturation_value;
            dc.per_image_TCSPC_FLIM_fit_model = handles.per_image_TCSPC_FLIM_fit_model;
            
            dc.per_image_TCSPC_FLIM_irf_filename = handles.per_image_TCSPC_FLIM_irf_filename;
            dc.per_image_TCSPC_FLIM_tvb_filename = handles.per_image_TCSPC_FLIM_tvb_filename;
            dc.per_image_TCSPC_FLIM_tvb_scaling = handles.per_image_TCSPC_FLIM_tvb_scaling;
            
            dc.per_image_TCSPC_FLIM_fixed_tauD = handles.per_image_TCSPC_FLIM_fixed_tauD;
            dc.per_image_TCSPC_FLIM_conv_irf_pp_69_70 = handles.per_image_TCSPC_FLIM_conv_irf_pp_69_70;
            
            dc.per_image_TCSPC_FLIM_Ref_lifetime = handles.per_image_TCSPC_FLIM_Ref_lifetime;
            dc.per_image_TCSPC_FLIM_averaging_sigma = handles.per_image_TCSPC_FLIM_averaging_sigma;
            
    fh = ancestor(hObject,'figure');     
    delete(fh);

% Cancel    
% --- Executes on button press in Cancel.
function Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        
    fh = ancestor(hObject,'figure');     
    delete(fh);



function irf_shift_Callback(hObject, eventdata, handles)
% hObject    handle to irf_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of irf_shift as text
%        str2double(get(hObject,'String')) returns contents of irf_shift as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.per_image_TCSPC_FLIM_irf_shift = value;
    guidata(hObject,handles);
else
    value = handles.per_image_TCSPC_FLIM_irf_shift;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function irf_shift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to irf_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fitting_model_Callback(hObject, eventdata, handles)
% hObject    handle to fitting_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fitting_model as text
%        str2double(get(hObject,'String')) returns contents of fitting_model as a double
% value = str2double(get(hObject,'String'));
% if ~isnan(value) && value <= 1. && value >= 0.
%     handles.NC_bckg_subtraction_proportion = value;
%     guidata(hObject,handles);
% else
%     value = handles.NC_bckg_subtraction_proportion;
%     set(hObject,'String',num2str(value));
% end

index = get(hObject,'Value');
str = get(hObject,'String');
handles.per_image_TCSPC_FLIM_fit_model = str(index);

guidata(hObject,handles);

uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function fitting_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fitting_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rep_rate_Callback(hObject, eventdata, handles)
% hObject    handle to rep_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rep_rate as text
%        str2double(get(hObject,'String')) returns contents of rep_rate as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 1 
    handles.per_image_TCSPC_FLIM_rep_rate = value;
    guidata(hObject,handles);
else
    value = handles.per_image_TCSPC_FLIM_rep_rate;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function rep_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rep_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in irf_load.
function irf_load_Callback(hObject, eventdata, handles)
% hObject    handle to irf_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.dc;
dc.load_irf();
if ~isempty(dc.per_image_TCSPC_FLIM_irf_filename)
    handles.per_image_TCSPC_FLIM_irf_filename = dc.per_image_TCSPC_FLIM_irf_filename;
    [~,IRF_NAME,IRF_EXT] = fileparts(handles.per_image_TCSPC_FLIM_irf_filename);
    set(handles.irf_name_text,'String',[IRF_NAME IRF_EXT]);    
end
guidata(hObject,handles);
uiresume(handles.figure1);


function dead_zone_left_Callback(hObject, eventdata, handles)
% hObject    handle to dead_zone_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dead_zone_left as text
%        str2double(get(hObject,'String')) returns contents of dead_zone_left as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 1 
    handles.per_image_TCSPC_FLIM_Tmin = value;
    guidata(hObject,handles);
else
    value = handles.per_image_TCSPC_FLIM_Tmin;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function dead_zone_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dead_zone_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dead_zone_right_Callback(hObject, eventdata, handles)
% hObject    handle to dead_zone_right_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dead_zone_right_text as text
%        str2double(get(hObject,'String')) returns contents of dead_zone_right_text as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 1 
    handles.per_image_TCSPC_FLIM_Tmax = value;
    guidata(hObject,handles);
else
    value = handles.per_image_TCSPC_FLIM_Tmax;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function dead_zone_right_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dead_zone_right_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function background_value_Callback(hObject, eventdata, handles)
% hObject    handle to background_value_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of background_value_text as text
%        str2double(get(hObject,'String')) returns contents of background_value_text as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 0 
    handles.per_image_TCSPC_FLIM_background_value = value;
    guidata(hObject,handles);
else
    value = handles.per_image_TCSPC_FLIM_background_value;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function background_value_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to background_value_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tvb_load.
function tvb_load_Callback(hObject, eventdata, handles)
% hObject    handle to tvb_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.dc;
dc.load_tvb();
if ~isempty(dc.per_image_TCSPC_FLIM_tvb_filename)
    handles.per_image_TCSPC_FLIM_tvb_filename = dc.per_image_TCSPC_FLIM_tvb_filename;        
    [~,TVB_NAME,TVB_EXT] = fileparts(handles.per_image_TCSPC_FLIM_tvb_filename);                
    set(handles.tvb_name_text,'String',[TVB_NAME TVB_EXT]);    
end
guidata(hObject,handles);
uiresume(handles.figure1);





% --- Executes on button press in apply.
function apply_Callback(hObject, eventdata, handles)
% hObject    handle to apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.dc;
        
            dc.per_image_TCSPC_FLIM_irf_shift = handles.per_image_TCSPC_FLIM_irf_shift;
            dc.per_image_TCSPC_FLIM_irf_background = handles.per_image_TCSPC_FLIM_irf_background;      
            dc.per_image_TCSPC_FLIM_rep_rate = handles.per_image_TCSPC_FLIM_rep_rate;
            %
            dc.per_image_TCSPC_FLIM_fit_model = handles.per_image_TCSPC_FLIM_fit_model;
            dc.per_image_TCSPC_FLIM_Tmin = handles.per_image_TCSPC_FLIM_Tmin;
            dc.per_image_TCSPC_FLIM_Tmax = handles.per_image_TCSPC_FLIM_Tmax;
            dc.per_image_TCSPC_FLIM_background_value = handles.per_image_TCSPC_FLIM_background_value;
            dc.per_image_TCSPC_FLIM_saturation_value = handles.per_image_TCSPC_FLIM_saturation_value;
            dc.per_image_TCSPC_FLIM_fit_model = handles.per_image_TCSPC_FLIM_fit_model;
            
            dc.per_image_TCSPC_FLIM_irf_filename = handles.per_image_TCSPC_FLIM_irf_filename;
            dc.per_image_TCSPC_FLIM_tvb_filename = handles.per_image_TCSPC_FLIM_tvb_filename;
            dc.per_image_TCSPC_FLIM_tvb_scaling = handles.per_image_TCSPC_FLIM_tvb_scaling;
            
            dc.per_image_TCSPC_FLIM_fixed_tauD = handles.per_image_TCSPC_FLIM_fixed_tauD;
            dc.per_image_TCSPC_FLIM_conv_irf_pp_69_70 = handles.per_image_TCSPC_FLIM_conv_irf_pp_69_70;
            
            dc.per_image_TCSPC_FLIM_Ref_lifetime = handles.per_image_TCSPC_FLIM_Ref_lifetime;
            dc.per_image_TCSPC_FLIM_averaging_sigma = handles.per_image_TCSPC_FLIM_averaging_sigma;

% --- Executes on button press in fixed_tauD_checkbox.
function fixed_tauD_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to fixed_tauD_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fixed_tauD_checkbox

toggle_state = get(hObject,'Value');

if ~toggle_state
    handles.per_image_TCSPC_FLIM_fixed_tauD = 0;
    set(handles.fixed_tauD,'Enable','off');
    set(handles.fixed_tauD,'String','...');
else
    set(handles.fixed_tauD,'Enable','on');
end

guidata(hObject,handles);
%uiresume(handles.figure1);


function fixed_tauD_Callback(hObject, eventdata, handles)
% hObject    handle to fixed_tauD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fixed_tauD as text
%        str2double(get(hObject,'String')) returns contents of fixed_tauD as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value > 0 
    handles.per_image_TCSPC_FLIM_fixed_tauD = value;
    guidata(hObject,handles);
else
    value = handles.per_image_TCSPC_FLIM_fixed_tauD;
    if 0~=value
        set(hObject,'String',num2str(value));
    else
        set(hObject,'Enable','off');
        set(hObject,'String','...'); 
        set(handles.fixed_tauD_checkbox,'Value',0);
    end
    guidata(hObject,handles);
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function fixed_tauD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixed_tauD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function saturation_value_Callback(hObject, eventdata, handles)
% hObject    handle to saturation_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of saturation_value as text
%        str2double(get(hObject,'String')) returns contents of saturation_value as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 0 
    handles.per_image_TCSPC_FLIM_saturation_value = value;
    guidata(hObject,handles);
else
    value = handles.per_image_TCSPC_FLIM_saturation_value;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function saturation_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saturation_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in conv_pp_69_70.
function conv_pp_69_70_Callback(hObject, eventdata, handles)
% hObject    handle to conv_pp_69_70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of conv_pp_69_70
toggle_state = get(hObject,'Value');
if ~toggle_state
    handles.per_image_TCSPC_FLIM_conv_irf_pp_69_70 = 0;
else
    handles.per_image_TCSPC_FLIM_conv_irf_pp_69_70 = 1;
end
guidata(hObject,handles);
%uiresume(handles.figure1);


function IRF_background_Callback(hObject, eventdata, handles)
% hObject    handle to IRF_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IRF_background as text
%        str2double(get(hObject,'String')) returns contents of IRF_background as a double
value = str2double(get(hObject,'String'));
if ~isnan(value)
    handles.per_image_TCSPC_FLIM_irf_background = value;
    guidata(hObject,handles);
else
    value = handles.per_image_TCSPC_FLIM_irf_background;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function IRF_background_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IRF_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tvb_scaling_Callback(hObject, eventdata, handles)
% hObject    handle to tvb_scaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tvb_scaling as text
%        str2double(get(hObject,'String')) returns contents of tvb_scaling as a double

% value = str2double(get(hObject,'String'));
value = eval(get(hObject,'String'));
if ~isnan(value) && value >= 0 
    handles.per_image_TCSPC_FLIM_tvb_scaling = value;
    guidata(hObject,handles);
else
    value = handles.per_image_TCSPC_FLIM_tvb_scaling;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function tvb_scaling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tvb_scaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IRF_reference_lifetime_Callback(hObject, eventdata, handles)
% hObject    handle to IRF_reference_lifetime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IRF_reference_lifetime as text
%        str2double(get(hObject,'String')) returns contents of IRF_reference_lifetime as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 0
    handles.per_image_TCSPC_FLIM_Ref_lifetime = value;
    guidata(hObject,handles);
else
    value = handles.per_image_TCSPC_FLIM_Ref_lifetime;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function IRF_reference_lifetime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IRF_reference_lifetime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function averaging_sigma_Callback(hObject, eventdata, handles)
% hObject    handle to averaging_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of averaging_sigma as text
%        str2double(get(hObject,'String')) returns contents of averaging_sigma as a double
value = str2double(get(hObject,'String'));
if ~isnan(value) && value >= 0 && value <= 21 
    handles.per_image_TCSPC_FLIM_averaging_sigma = value;
    guidata(hObject,handles);
else
    value = handles.per_image_TCSPC_FLIM_averaging_sigma;
    set(hObject,'String',num2str(value));
end
uiresume(handles.figure1);

% --- Executes during object creation, after setting all properties.
function averaging_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to averaging_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
