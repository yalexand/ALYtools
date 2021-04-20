function varargout = MicroscopyImageFormatter(varargin)
% MICROSCOPYIMAGEFORMATTER MATLAB code for MicroscopyImageFormatter.fig
%      MICROSCOPYIMAGEFORMATTER, by itself, creates a new MICROSCOPYIMAGEFORMATTER or raises the existing
%      singleton*.
%
%      H = MICROSCOPYIMAGEFORMATTER returns the handle to a new MICROSCOPYIMAGEFORMATTER or the handle to
%      the existing singleton*.
%
%      MICROSCOPYIMAGEFORMATTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MICROSCOPYIMAGEFORMATTER.M with the given input arguments.
%
%      MICROSCOPYIMAGEFORMATTER('Property','Value',...) creates a new MICROSCOPYIMAGEFORMATTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MicroscopyImageFormatter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MicroscopyImageFormatter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MicroscopyImageFormatter

% Last Modified by GUIDE v2.5 29-Mar-2021 11:17:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MicroscopyImageFormatter_OpeningFcn, ...
                   'gui_OutputFcn',  @MicroscopyImageFormatter_OutputFcn, ...
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


% --- Executes just before MicroscopyImageFormatter is made visible.
function MicroscopyImageFormatter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MicroscopyImageFormatter (see VARARGIN)

set(handles.image_type,'String',{'Generic','Optosplit 2 channels','Optosplit 3 channels'});

handles.Optosplit_registration_roix = [];
handles.Optosplit_registration_roiy1 = [];
handles.Optosplit_registration_roiy2 = [];
handles.Optosplit_registration_tform = [];
handles.Optosplit_registration_droi_x = [];
handles.Optosplit_registration_droi_y = [];
handles.Optosplit_registration_a3 = 270; %degrees
handles.Optosplit_registration_f3 = 1/1.8; 
handles.Optosplit_registration_tform3 = [];
handles.Optosplit_registration_roi3x = [];
handles.Optosplit_registration_roi3y = [];

    set(handles.src_dir,'String',['c:' filesep]);
    set(handles.dst_dir,'String',['c:' filesep]);
    set(handles.ref_image_file,'String','xxx');
    % values
    handles.umppix = 0.650;
    handles.offset = 0; % !!!
    handles.downsample = 1;
    handles.min_per_frame = 5;
    src_channels = '12345';
    dst_channels = '12345';
    set(handles.setup_Optosplit_registration,'Enable','Off');
    set(handles.setup_Optosplit_registration,'Visible','Off');
    handles.polynom_order = 9;

set(handles.umppix_edit,'String',num2str(handles.umppix));
set(handles.offset_edit,'String',num2str(handles.offset));
set(handles.downsample_edit,'String',num2str(handles.downsample));
set(handles.min_per_frame_edit,'String',num2str(handles.min_per_frame));
set(handles.src_channels,'String',src_channels);
set(handles.dst_channels,'String',dst_channels);
set(handles.t_dep_fitting_poly_order,'String',num2str(handles.polynom_order));
set(handles.clean_reference,'Value',1);

set(handles.image_raw,'XTick',[]);
set(handles.image_raw,'YTick',[]);
set(handles.image_corrected,'XTick',[]);
set(handles.image_corrected,'YTick',[]);

set(handles.show_channel,'String',{'1','2','3','4','5','All'});

set(handles.compensate_CMHF,'Value',1);
handles.f_CMHF_t = ones(100000,1); % :) should be enough..

handles.ref_img = [];
handles.raw_img = [];
handles.corrected_img = [];

handles.raw_image_filename = [];

handles.p_xy = cell(0);
handles.f_t = cell(0);
handles.Eb = cell(0);

handles.model_icon = cell(4,1);
[~,~,handles.model_icon{1}] = bfopen_v([pwd filesep 'GUIDEInterfaces' filesep 'Formatter_model1.tif']);
[~,~,handles.model_icon{2}] = bfopen_v([pwd filesep 'GUIDEInterfaces' filesep 'Formatter_model2.tif']);
[~,~,handles.model_icon{3}] = bfopen_v([pwd filesep 'GUIDEInterfaces' filesep 'Formatter_model3.tif']);
[~,~,handles.model_icon{4}] = bfopen_v([pwd filesep 'GUIDEInterfaces' filesep 'Formatter_model4.tif']);

imshow(handles.model_icon{1}, 'Parent', handles.icon_1);
imshow(handles.model_icon{1}, 'Parent', handles.icon_2);
imshow(handles.model_icon{1}, 'Parent', handles.icon_3);
imshow(handles.model_icon{1}, 'Parent', handles.icon_4);
imshow(handles.model_icon{1}, 'Parent', handles.icon_5);

bckg_color = get(handles.figure1,'Color');
    R = bckg_color(1)*ones(size(handles.model_icon{1}));
    G = bckg_color(2)*ones(size(handles.model_icon{1}));
    B = bckg_color(3)*ones(size(handles.model_icon{1}));    
handles.dummy = cat(3,R,G,B);

handles.W = []; % cross-talk matrix
handles.g = []; % g factor rendering intensities in the same units
set(handles.cross,'Visible','off');
set(handles.talk,'Visible','off');
set(handles.correction,'Visible','off');
set(handles.unset_spectral_cross_talk_correction,'Enable','off');
set(handles.unset_spectral_cross_talk_correction,'Visible','off');

% Choose default command line output for MicroscopyImageFormatter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MicroscopyImageFormatter wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MicroscopyImageFormatter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function umppix_edit_Callback(hObject, eventdata, handles)
    v = str2double(get(hObject,'String'));
    if ~isempty(v) && v>0 
        handles.umppix = v;        
    else
        set(hObject,'String',num2str(handles.umppix));
    end
    guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function umppix_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to umppix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function min_per_frame_edit_Callback(hObject, eventdata, handles)
    v = str2double(get(hObject,'String'));
    if ~isempty(v) && v>0 
        handles.min_per_frame = v;        
    else
        set(hObject,'String',num2str(handles.min_per_frame));
    end
    guidata(hObject,handles);

        
% --- Executes during object creation, after setting all properties.
function min_per_frame_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_per_frame_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function downsample_edit_Callback(hObject, eventdata, handles)
    v = str2double(get(hObject,'String'));
    if ~isempty(v) && v>0 
        handles.downsample = v;        
    else
        set(hObject,'String',num2str(handles.downsample));
    end
    guidata(hObject,handles);
    
    

% --- Executes during object creation, after setting all properties.
function downsample_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to downsample_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function src_dir_Callback(hObject, eventdata, handles)
% hObject    handle to src_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of src_dir as text
%        str2double(get(hObject,'String')) returns contents of src_dir as a double


% --- Executes during object creation, after setting all properties.
function src_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to src_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dst_dir_Callback(hObject, eventdata, handles)
% hObject    handle to dst_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dst_dir as text
%        str2double(get(hObject,'String')) returns contents of dst_dir as a double


% --- Executes during object creation, after setting all properties.
function dst_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dst_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ref_image_file_Callback(hObject, eventdata, handles)
% hObject    handle to ref_image_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ref_image_file as text
%        str2double(get(hObject,'String')) returns contents of ref_image_file as a double


% --- Executes during object creation, after setting all properties.
function ref_image_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ref_image_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function src_channels_Callback(hObject, eventdata, handles)
s = get(hObject,'String');
if ~ismember(s,{'1','12','123','1234','12345'})
    set(hObject,'String','xxxxx');
end
guidata(hObject,handles);
    


% --- Executes during object creation, after setting all properties.
function src_channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to src_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dst_channels_Callback(hObject, eventdata, handles)
dst_channels = get(hObject,'String');
%
L_src = length(get(handles.src_channels,'String')); % may be <= 5
L_dst = length(dst_channels);
%
if L_dst<=L_src
    dst_channels_OK = true;    
    for k=1:L_dst
        ch_k = str2num(dst_channels(k));
        if isempty(ch_k) || ch_k>L_src
            dst_channels_OK = false;
            break;
        end
    end
else
    dst_channels_OK = false;
end
if ~dst_channels_OK
    set(hObject,'String','xxxxx');
else
            % clean images
            handles.ref_img = [];
            handles.raw_img = [];
            handles.corrected_img = [];
            %
            cla(handles.image_raw,'reset');
            cla(handles.image_corrected,'reset');    
end
setup_model_controls_visibility(handles);
guidata(hObject,handles);
    

% --- Executes during object creation, after setting all properties.
function dst_channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dst_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in model_1.
function model_1_Callback(hObject, eventdata, handles)
    imshow(handles.model_icon{get(hObject,'Value')}, 'Parent', handles.icon_1);
    guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function model_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in model_2.
function model_2_Callback(hObject, eventdata, handles)
    imshow(handles.model_icon{get(hObject,'Value')}, 'Parent', handles.icon_2);
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function model_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in model_3.
function model_3_Callback(hObject, eventdata, handles)
    imshow(handles.model_icon{get(hObject,'Value')}, 'Parent', handles.icon_3);
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function model_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in model_4.
function model_4_Callback(hObject, eventdata, handles)
    imshow(handles.model_icon{get(hObject,'Value')}, 'Parent', handles.icon_4);
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function model_4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in model_5.
function model_5_Callback(hObject, eventdata, handles)
    imshow(handles.model_icon{get(hObject,'Value')}, 'Parent', handles.icon_5);
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function model_5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function offset_edit_Callback(hObject, eventdata, handles)
    v = str2double(get(hObject,'String'));
    if ~isempty(v) && v>0 
        handles.offset = v;        
    else
        set(hObject,'String',num2str(handles.offset));
    end
    guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function offset_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offset_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clean_reference.
function clean_reference_Callback(hObject, eventdata, handles)
% hObject    handle to clean_reference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of clean_reference


% --- Executes on button press in load_image.
function load_image_Callback(hObject, eventdata, handles)

            [filename,pathname] = uigetfile({'*.tif;*.tiff','Image Files'}, ...
                'Select image file',get(handles.src_dir,'String'));
            if filename == 0, return, end

            full_path_to_file =  [pathname filesep filename];           

            v = load_microscopy_image(handles,full_path_to_file);
            if isempty(v), return, end
            %
            handles.raw_image_filename = filename;
            %
            handles.raw_img = v;
            guidata(hObject,handles);          
            %
            [~,FNAME,FEXT] = fileparts(full_path_to_file);
            show_image(handles,'raw_img','image_raw',[FNAME FEXT]);
            show_image(handles,'raw_img','image_corrected',[]);                
        
% --- Executes on button press in recalculate_corrections.
function recalculate_corrections_Callback(hObject, eventdata, handles)
% applies pre-calculated corrections to the loaded image
%
if isempty(handles.raw_img), return, end

[SX,SY,n_channels,~,st] = size(handles.raw_img);
handles.corrected_img = zeros(size(handles.raw_img));
for c = 1:n_channels
    
    s = get(handles.model_1,'String'); % all the same
    switch c
        case 1
            model = s{get(handles.model_1,'Value')};
        case 2
            model = s{get(handles.model_2,'Value')};            
        case 3
            model = s{get(handles.model_3,'Value')};            
        case 4
            model = s{get(handles.model_4,'Value')};            
        case 5            
            model = s{get(handles.model_5,'Value')};            
    end %switch
    %
    Eb = handles.Eb{c};
    p_xy  = handles.p_xy{c};
    f_t = handles.f_t{c};
            %
            if 1==get(handles.compensate_CMHF,'Value')
                f_CMHF_t = handles.f_CMHF_t;
            else
                f_CMHF_t = ones(size(handles.f_CMHF_t));
            end            
            %
    hw = waitbar(0,['introducing corrections to channel ' num2str(c)]);
    for f = 1:st
        I = squeeze(handles.raw_img(:,:,c,1,f)) - handles.offset;
        I = I/f_CMHF_t(f);
        switch model
            case 'additive (1)'         % f(t) acts only on background
                EO = I - Eb*f_t(f)*p_xy;
            case 'multiplicative (1)'   
                EO = I./p_xy - Eb*f_t(f);                
            case 'additive (2)'         % f(t) acts on everyhting
                EO = I/f_t(f) - Eb*p_xy;                
            case 'multiplicative (2)'   
                EO = I./( f_t(f)*p_xy ) - Eb;
        end
        EO(EO<0)=0;
        handles.corrected_img(:,:,c,1,f) = EO;
        if ~isempty(hw), waitbar(f/st,hw); end
    end
    if ~isempty(hw), delete(hw), drawnow; end
end
%
if ~isempty(handles.W)
    handles.corrected_img = introduce_cross_talk_corrections(handles.corrected_img,handles);
end
%
guidata(hObject,handles);

show_image(handles,'corrected_img','image_corrected',[]);


function frame_to_show_Callback(hObject, eventdata, handles)
        show_image(handles,'raw_img','image_raw',[]);
        show_image(handles,'corrected_img','image_corrected',[]);

% --- Executes during object creation, after setting all properties.
function frame_to_show_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_to_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in show_channel.
function show_channel_Callback(hObject, eventdata, handles)
        show_image(handles,'raw_img','image_raw',[]);
        show_image(handles,'corrected_img','image_corrected',[]);

% --- Executes during object creation, after setting all properties.
function show_channel_CreateFcn(hObject, eventdata, ~)
% hObject    handle to show_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t_dep_fitting_poly_order_Callback(hObject, eventdata, handles)
    v = int64(str2double(get(hObject,'String')));
    if ~isempty(v) && isinteger(v) && v>=1 
        handles.polynom_order = v;        
    else
        set(hObject,'String',num2str(handles.polynom_order));
    end
    guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function t_dep_fitting_poly_order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_dep_fitting_poly_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in process_images.
function process_images_Callback(hObject, eventdata, handles)

[filenames,pathname] = uigetfile('*.tif','Select image files',get(handles.src_dir,'String'),'MultiSelect','on');                
if isempty(filenames), return, end       
if isnumeric(filenames) && 0==filenames, return, end

if ischar(filenames)
    filenames = {filenames};
end
for k=1:numel(filenames)

    try
        handles.raw_img = load_microscopy_image(handles,[get(handles.src_dir,'String') filesep filenames{k}]);
    catch
        disp(['error when loading ' filenames{k}]);
        continue;
    end
    
        [SX,SY,n_channels,~,st] = size(handles.raw_img);
        handles.corrected_img = zeros(size(handles.raw_img));
        for c = 1:n_channels
            %            
            s = get(handles.model_1,'String'); % all the same
            switch c
                case 1
                    model = s{get(handles.model_1,'Value')};
                case 2
                    model = s{get(handles.model_2,'Value')};            
                case 3
                    model = s{get(handles.model_3,'Value')};            
                case 4
                    model = s{get(handles.model_4,'Value')};            
                case 5            
                    model = s{get(handles.model_5,'Value')};            
            end %switch
            %
            Eb = handles.Eb{c};
            p_xy  = handles.p_xy{c};
            f_t = handles.f_t{c};
            %
            if 1==get(handles.compensate_CMHF,'Value')
                f_CMHF_t = handles.f_CMHF_t;
            else
                f_CMHF_t = ones(size(handles.f_CMHF_t));
            end            
            %
            hw = waitbar(0,['introducing corrections to channel ' num2str(c)]);
            for f = 1:st
                I = squeeze(handles.raw_img(:,:,c,1,f)) - handles.offset;
                I = I/f_CMHF_t(f);
                switch model
                    case 'additive (1)'         % f(t) acts only on background
                        EO = I - Eb*f_t(f)*p_xy;
                    case 'multiplicative (1)'   
                        EO = I./p_xy - Eb*f_t(f);                
                    case 'additive (2)'         % f(t) acts on everyhting
                        EO = I/f_t(f) - Eb*p_xy;                
                    case 'multiplicative (2)'   
                        EO = I./( f_t(f)*p_xy ) - Eb;
                end
                EO(EO<0)=0;
                handles.corrected_img(:,:,c,1,f) = EO;
                if ~isempty(hw), waitbar(f/st,hw); end
            end
            if ~isempty(hw), delete(hw), drawnow; end
        end
        %            
        if ~isempty(handles.W)
           handles.corrected_img = introduce_cross_talk_corrections(handles.corrected_img,handles);
        end        
        %
        guidata(hObject,handles);
        
        show_image(handles,'corrected_img','image_corrected',filenames{k});
        show_image(handles,'raw_img','image_raw',[]);
        
        guidata(hObject,handles);        
        
        savename = filenames{k};
        if ~contains(savename,'ome')
            savename = strrep(savename,'.tif','.ome.tif'); % aaaaa....
        end
        bfsave(uint16(handles.corrected_img),[get(handles.dst_dir,'String') filesep savename],'Compression','LZW','BigTiff', true,'dimensionOrder','XYCZT');
end

%--------------------------------------------------------------
function v = load_microscopy_image(handles,full_path_to_file)
    s = get(handles.image_type,'String');
    switch char(s(get(handles.image_type,'Value')))
        case 'Optosplit 2 channels'                      
            v = load_Optosplit_image(handles,full_path_to_file);
        case 'Optosplit 3 channels'
            v = load_Optosplit_image_3(handles,full_path_to_file);
        case 'Generic'
            v = load_Generic_image(handles,full_path_to_file);
    end
    setup_model_controls_visibility(handles);


% --- Executes on button press in load_ref.
function load_ref_Callback(hObject, eventdata, handles)

    full_path_to_file = get(handles.ref_image_file,'String');

        v = load_microscopy_image(handles,full_path_to_file);
        if isempty(v), return, end
        handles.ref_img = single(v);
        handles.raw_img = single(v);
        guidata(hObject,handles);
        %
        [~,FNAME,FEXT] = fileparts(full_path_to_file);
        show_image(handles,'raw_img','image_raw',[FNAME FEXT]);
        show_image(handles,'raw_img','image_corrected',[]);
        
        get_correction_functions_from_ref(hObject,handles);

%--------------------------------------------------------------    
function show_image(handles,what_to_show,where_to_show,TITLE)    

    image = eval(['handles.' what_to_show]);
    ax = eval(['handles.' where_to_show]);
    
    if isempty(image), return, end
    
        frame = str2num(get(handles.frame_to_show,'String'));
        if isempty(frame) || frame<=0 || frame > size(image,5)
            frame = 1;
            set(handles.frame_to_show,'String',num2str(frame))
        end
                        
        s = get(handles.show_channel,'String');
        ind = get(handles.show_channel,'Value');
        current_channel_to_show = str2num(s{ind});
        %
        % safety
        n_channels = size(image,3);
        if current_channel_to_show > n_channels
            set(handles.show_channel,'Value',1);
            set(handles.show_channel,'String',s(1:n_channels));
            current_channel_to_show = 1;
        end
        %        
            img = single(image(:,:,current_channel_to_show,1,frame));
            %
             t = quantile(img(:),0.999);
             img(img>t) = t;
            %
            imshow(uint8(map(img,0,255)), 'Parent', ax);

            if ~isempty(TITLE)
                title(ax,TITLE);
            end
    
%--------------------------------------------------------------
function v = load_Optosplit_image(handles,full_path_to_image)

v = [];

roix = handles.Optosplit_registration_roix;
roiy1 = handles.Optosplit_registration_roiy1;
roiy2 = handles.Optosplit_registration_roiy2;
tform = handles.Optosplit_registration_tform;
droi_x = handles.Optosplit_registration_droi_x;
droi_y = handles.Optosplit_registration_droi_y;

if isempty(roix)
    disp('cannot load Optosplit image - registration not set up');
    return;
end

    [~,~,I] = bfopen_v(full_path_to_image);
    [sx,sy,sc,sz,st] = size(I);
    %
        s = get(handles.dst_channels,'String');
            nch1 = int64(str2num(s(1)));
            nch2 = int64(str2num(s(2)));
    %
    if length(s)>2
        s = s(1:2);
        set(handles.dst_channels,'String',s);
    end
    %
    for f=1:st
        u = single(squeeze(I(:,:,1,1,f)));
        u1 = u(roix,roiy1);
        u2 = u(roix,roiy2);
        
        u2_reg = imwarp(u2,tform,'OutputView',imref2d(size(u1)));
        
        u1_ = u1(droi_x,droi_y);
        u2_ = u2_reg(droi_x,droi_y);
        if isempty(v)
            v = zeros(size(u1_,1),size(u1_,2),2,1,st);
        end
        v(:,:,nch1,1,f) = u1_;
        v(:,:,nch2,1,f) = u2_;            
    end
    %

% --- Executes on button press in set_src_dir.
function set_src_dir_Callback(hObject, eventdata, handles)
directoryname = uigetdir(pwd,'Pick SRC Directory');
if 0==directoryname, return, end
if isfolder(directoryname) 
    set(handles.src_dir,'String',directoryname);
    guidata(hObject,handles);
end

% --- Executes on button press in set_dst_dir.
function set_dst_dir_Callback(hObject, eventdata, handles)
directoryname = uigetdir(pwd,'Pick DST Directory');
if 0==directoryname, return, end
if isfolder(directoryname) 
    set(handles.dst_dir,'String',directoryname);
    guidata(hObject,handles);
end

% --- Executes on button press in set_ref_image_file.
function set_ref_image_file_Callback(hObject, eventdata, handles)
    [filename,pathname] = uigetfile({'*.tif;*.tiff','Image Files'}, ...
    'Select image file',get(handles.src_dir,'String'));
    if filename == 0, return, end
    full_path_to_example =  [pathname filesep filename];           
    %
    set(handles.ref_image_file,'String',full_path_to_example);
    guidata(hObject,handles);

    
% --- Executes on button press in send_raw_to_Icy.
function send_raw_to_Icy_Callback(hObject, eventdata, handles)
    send_image_to_Icy(handles.raw_img,'raw');

% --- Executes on button press in send_corrected_to_Icy.
function send_corrected_to_Icy_Callback(hObject, eventdata, handles)
    send_image_to_Icy(handles.corrected_img,'corrected');
        

%-------------------------------------------------------------- 
    function send_image_to_Icy(img,title)    
    if isempty(img),return, end
    if ~ischar(title)
        title = '';
    end
    try
        icy_imshow(uint16(img),title);
    catch
        disp('cannot send whole image to Icy, send first 20 frames');
        img = img(:,:,:,:,1:20);
        icy_imshow(uint16(img),title);
    end
        
    
%-------------------------------------------------------------- 
function get_correction_functions_from_ref(hObject,handles)

offset = handles.offset;
polynom_order  = handles.polynom_order;

[SX,SY,n_channels,~,st] = size(handles.ref_img);

colors = {'r','g','b','k','c'};
AXES = handles.time_dependence_corr;
reset(AXES);
LEGEND = cell(0);

for channel = 1:n_channels

    if size(handles.ref_img,5)==1 % specific case when no t-dependence is available
        handles.f_t{channel} = ones(100000,1); % :) should be enough..      
        %handles.p_xy{channel} = ones(size(handles.ref_img(),1),size(handles.ref_img(),2));
        %handles.Eb{channel} = 0;                        
            u = handles.ref_img(:,:,channel,1,1);
            % calculate normalized profile
            prof = u - offset;    
            % smooth
                smooth_scale = 10;
                prof = imopen(prof,strel('disk',smooth_scale,0));
            if get(handles.clean_reference,'Value')
                prof = gsderiv(prof,smooth_scale,0);
            end
            [xmax,ymax] = find(prof==max(prof(:)));
            prof = prof/prof(xmax(1),ymax(1));
            handles.p_xy{channel} = prof;
            
            % calculate Eb
            ref = squeeze(handles.ref_img(:,:,channel,1,1));
            ds = 10;
            rx = (xmax-ds):(xmax+ds);
            ry = (ymax-ds):(ymax+ds);    
            rx(rx<1)=[];    
            ry(ry<1)=[];
            rx(rx>size(ref,1))=[];
            ry(ry>size(ref,2))=[];            
            %
            sample = ref(rx,ry,1);
            Eb = mean(sample(:)) - offset;                                
            handles.Eb{channel} = Eb;
    else
    
    ref = squeeze(handles.ref_img(:,:,channel,1,:));
    
    intensity = zeros(st,1);

    w = 40; % exclude noisy pixels?
    rx = w:SX-w;
    ry = w:SY-w;

    parfor f=1:st
        intensity(f) = mean(ref(rx,ry,f),'All') - offset;
    end

    frms = (1:st)';
    p = polyfit(frms,intensity,polynom_order);
    intensity_fit = polyval(p,frms);
    f_t = intensity_fit/intensity_fit(1);

    taxis = (frms-1)*handles.min_per_frame/60;    
    semilogy(AXES,taxis,intensity,[colors{channel} '.-'],taxis,intensity_fit,'k:','linewidth',3);
    hold(AXES,'on');
    LEGEND = [LEGEND num2str(channel) 'fit'];

    % calculate normalized profile
    prof = squeeze(sum(ref(:,:,1:st),3)) - offset;    
    % smooth
    smooth_scale = 3;    
    prof = imopen(prof,strel('disk',smooth_scale,0));
    if get(handles.clean_reference,'Value')
        prof = gsderiv(prof,smooth_scale,0);
    end

    [xmax,ymax] = find(prof==max(prof(:)));
    prof = prof/prof(xmax(1),ymax(1));
    %icy_imshow(prof,num2str(channel));
    
    % calculate Eb    
    ds = 10;
    rx = (xmax-ds):(xmax+ds);
    ry = (ymax-ds):(ymax+ds);    
    rx(rx<1)=[];    
    ry(ry<1)=[];           
    %
    sample = ref(rx,ry,1:3);
    Eb = mean(sample(:)) - offset;
           
    handles.f_t{channel} = f_t;
    handles.p_xy{channel} = prof;
    handles.Eb{channel} = Eb;
    
    end
        
end

% define CMHF correction (Common Multiplicative HF)
    w = 40; 
    rx = w:SX-w;
    ry = w:SY-w;
        ref = squeeze(sum(handles.ref_img,3)) - n_channels*offset;
        parfor f=1:st
            intensity(f) = mean(ref(rx,ry,f),'All');
        end
    frms = (1:st)';
    p = polyfit(frms,intensity,polynom_order);
    intensity_fit = polyval(p,frms);
    handles.f_CMHF_t = intensity./intensity_fit; % one needs to divide by this after offset subtraction, to introduce correction    
    %
%     taxis = (frms-1)*handles.min_per_frame/60;
%     figure(22);semilogy(taxis,intensity,'b.-',taxis,intensity_fit,'r.-');grid(gca,'on'); 
%     figure(23);semilogy(taxis,intensity./intensity_fit,'k.-');grid(gca,'on'); 
      %

guidata(hObject,handles);

hold(AXES,'off');
    xlabel(AXES,'time [h]','fontsize',8);
    ylabel(AXES,'offset-subtracted mean ref. intensity','fontsize',8);
    grid(AXES,'on');
legend(AXES,LEGEND);


% --- Executes on selection change in image_type.
function image_type_Callback(hObject, eventdata, handles)
    if ismember(get(hObject,'Value'),[2 3]) % Optosplit
        flag = 'On';
    else
        flag = 'Off';
    end
    set(handles.setup_Optosplit_registration,'Enable',flag);
    set(handles.setup_Optosplit_registration,'Visible',flag);    
    guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function image_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------
function v = load_Generic_image(handles,full_path_to_file)
v = [];

if ~isfile(full_path_to_file), return, end

s = get(handles.src_channels,'String');
n_src_channels = numel(s);
src_channels = zeros(n_src_channels,1);
for k=1:n_src_channels
    src_channels(k) = str2num(s(k));
end
s = get(handles.dst_channels,'String');
n_dst_channels = numel(s);
dst_channels = zeros(n_dst_channels,1);
for k=1:n_dst_channels
    dst_channels(k) = str2num(s(k));
end

[~,~,I] = bfopen_v(full_path_to_file);
[sx_src,sy_src,~,st_src] = size(I);

sx = sx_src;
sy = sy_src;
f = handles.downsample;
if 1~=f
    u = squeeze(single(I(:,:,1,1)));
    [sx,sy] = size(imresize(u,1/f));    
end
%
st = st_src/n_src_channels;
%
v = zeros(sx,sy,n_dst_channels,1,st);

for k=1:st
    i_b = n_src_channels*(k-1)+1;
    i_e = i_b + n_src_channels-1;
    slice = I(:,:,1,i_b:i_e);
    for m=1:n_dst_channels
            u = squeeze(slice(:,:,dst_channels(m)));
        if 1~=f
            u = imresize(u,1/f);
        end
        v(:,:,m,1,k) = u;
    end
end

s = cell(0);
for k=1:n_dst_channels, s = [s num2str(k)]; end
set(handles.show_channel,'String',s);

%--------------------------------------------------------------
function [roix,roiy1,roiy2,tform,droi_x,droi_y] = get_Optosplit_registration_parameters(full_path_to_file)

% this is 1 channel, multi-frame image, where 2 channels are imaged in
% every frame side to side
[~,~,I] = bfopen_v(full_path_to_file);

[SX,SY,~,~,st] = size(I);

u = squeeze(sum(single(I),5));

D = 10;
YC = round(SY/2);
roiy1 = D:(YC-D);
roiy2 = (YC+D):(SY-D);
roix = D:(SX-D);
u1 = u(roix,roiy1);
    deltax = 0; % to debug 
    deltay = 0;
u2 = u(roix+deltax,roiy2+deltay);

u1 = u1/mean(u1(:));
u2 = u2/mean(u2(:));

     fixed = u1;
     warped = u2;
     
     fixed_plus = fixed + gsderiv(fixed,2,0);
     warped_plus = warped + gsderiv(warped,2,0);
     
     [optimizer, metric] = imregconfig('multimodal');
     optimizer.MaximumIterations = 600;
     
     translation_only = false;
     if translation_only
         tform = imregtform(warped_plus,fixed_plus,'translation', optimizer, metric);
         % can visualize at this stage
         % registered = imwarp(warped,tform,'OutputView',imref2d(size(fixed)));    

         dy = round(tform.T(3,1));
         dx = round(tform.T(3,2));
         disp([dx dy]);

        % design of diminished ROIs
        roi_x = 1:size(u1,1);
        droi_x = roi_x;
        if dx>=0
            droi_x = (dx+2):roi_x(end);
        elseif dx<0
            droi_x = 1:(roi_x(end)-abs(dx)-2);
        end
        %
        roi_y = 1:size(u1,2);
        droi_y = roi_y;
        if dy>=0
            droi_y = (dy+2):roi_y(end);
        elseif dy<0
            droi_y = 1:(roi_y(end)-abs(dy)-2);
        end
     else
         tform = imregtform(warped_plus,fixed_plus,'rigid', optimizer, metric);         
         % can visualize at this stage
         registered = imwarp(warped,tform,'OutputView',imref2d(size(fixed)));    
         %
         tform.T
         %
        % finding diminished ROIs
         z=registered>0;
         S = FindLargestSquares(z);
         [C, H, W, M] = FindLargestRectangles(z,[0 0 1]);
         % icy_imshow(z+M);
         [~, pos] = max(C(:));
         [r, c] = ind2sub(size(S), pos);
         droi_y = c:W(r,c);
         droi_x = r:H(r,c);         
     end                 
%   
% visualize
%     fixed = fixed(droi_x,droi_y);
%     warped = warped(droi_x,droi_y);
%     registered = registered(droi_x,droi_y);    
%     iv = zeros(size(fixed,1),size(fixed,2),3,1,1);            
%     iv(:,:,1,1,1) = fixed;
%     iv(:,:,2,1,1) = warped;
%     iv(:,:,3,1,1) = registered;        
%     icy_imshow(iv);

%     v = [];
%     for f=1:st
%         u = squeeze(I(:,:,1,1,f));
%         u1 = u(roix,roiy1);
%         u2 = u(roix,roiy2);
%         
%         u2_reg = imwarp(u2,tform,'OutputView',imref2d(size(u1)));
%         
%         u1_ = u1(droi_x,droi_y);
%         u2_ = u2_reg(droi_x,droi_y);
%         if isempty(v)
%             v = zeros(size(u1_,1),size(u1_,2),2,1,st);
%         end
%         v(:,:,1,1,f) = u1_;
%         v(:,:,2,1,f) = u2_;            
%     end
%     %icy_imshow(uint16(v),['dx = ' num2str(dx) ', dy = ' num2str(dy)]);
%     icy_imshow(uint16(v));

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function setup_Optosplit_registration_Callback(hObject, eventdata, handles)

            [filename,pathname] = uigetfile({'*.tif;*.tiff','Image Files'}, ...
                'Select image file',get(handles.src_dir,'String'));
            if filename == 0, return, end
            full_path_to_example =  [pathname filesep filename];           
    %
    image_type = get(handles.image_type,'Value');
    switch image_type
        case 2 % 2-channel Optosplit
            [handles.Optosplit_registration_roix, ...
            handles.Optosplit_registration_roiy1, ...
            handles.Optosplit_registration_roiy2, ...
            handles.Optosplit_registration_tform, ...
            handles.Optosplit_registration_droi_x, ...
            handles.Optosplit_registration_droi_y] = get_Optosplit_registration_parameters(full_path_to_example);
            %
            set(handles.show_channel,'String',{'1','2'});
            set(handles.src_channels,'String','12');
            set(handles.dst_channels,'String','12');    
        case 3 % 3-channel Optosplit
            [handles.Optosplit_registration_roix, ...
            handles.Optosplit_registration_roiy1, ...
            handles.Optosplit_registration_roiy2, ...
            handles.Optosplit_registration_tform, ...
            handles.Optosplit_registration_droi_x, ...
            handles.Optosplit_registration_droi_y, ...
            handles.Optosplit_registration_tform3, ...
            handles.Optosplit_registration_roi3x, ...
            handles.Optosplit_registration_roi3y] = get_Optosplit_3_registration_parameters(full_path_to_example,handles);
            %
            set(handles.show_channel,'String',{'1','2','3'});
            set(handles.src_channels,'String','123');
            set(handles.dst_channels,'String','123');                
    end
    guidata(hObject,handles);
           
% --------------------------------------------------------------------
function derive_corrections_from_data_Callback(hObject, eventdata, handles)

[filenames,pathname] = uigetfile('*.tif','Select image files',get(handles.src_dir,'String'),'MultiSelect','on');                
if isempty(filenames), return, end       
if isnumeric(filenames) && 0==filenames, return, end

%
if ischar(filenames)
    filenames = {filenames};
end

offset = handles.offset;

% % WHY parfor IS IT SO UNSTABLE? 
img_acc = cell(numel(filenames),1);
%parfor k=1:numel(filenames)
for k=1:numel(filenames)
    raw_img = [];
    try
        raw_img = load_microscopy_image(handles,[get(handles.src_dir,'String') filesep filenames{k}]);
    catch
        disp(['error when loading ' filenames{k}]);
        continue;
    end
    img_acc{k} = raw_img;
    [k numel(filenames)]
end

img_data = [];
for k=1:numel(filenames)
    img_data = cat(5,img_data,single(img_acc{k}));
    size(img_data)
end

[sx,sy,sc,~,~] = size(img_data);
%
img_data = img_data - offset;
% ??
t = quantile(img_data(:),0.99);
img_data(img_data>t)=t;

ref = zeros(sx,sy,sc,1,1);
%
for c=1:sc
    ref_c = zeros(sx,sy);
    for x=1:sx
        disp(x);
        parfor y=1:sy
            s = img_data(x,y,c,1,:);
            t = quantile(s,0.01); % take minimum value
            ref_c(x,y)=t;
        end
    end
    ref(:,:,c) = ref_c;
end
%
% a bit of correction. as we need rather average, not minimum
% we presume Poissonic noise, i.e. that [ref_measured = ref_true - sqrt(ref_true)] and then find ref_true
ref = 1/2*( 1 + 2*ref + sqrt(4*ref + 1) );
clear('img_data');
%
% create dimensionless p_xy correction images for channels
xmax = zeros(sc,1);
ymax = zeros(sc,1);
for channel = 1:sc
            prof = ref(:,:,channel,1,1);
            % calculate normalized profile
            % smooth            
                smooth_scale = 10;
                prof = imopen(prof,strel('disk',smooth_scale,0));
            if get(handles.clean_reference,'Value')                
                prof = gsderiv(prof,smooth_scale,0);
            end
            [xmax(channel),ymax(channel)] = find(prof==max(prof(:)));
            prof = prof/prof(xmax(1),ymax(1));            
%                 icy_imshow(handles.p_xy{channel},['reference ' num2str(channel)]);
%                 icy_imshow(prof,['derived from data ' num2str(channel)]);            
            handles.p_xy{channel} = prof;
end
%
% work on Eb and temporal dependencies
f_data = zeros(size(img_acc,1),sc,size(img_acc{1},5));
st = size(img_acc{1},5);

    for k=1:size(img_acc,1)  
        k
        for c=1:sc
            %            
            u = img_acc{k};            
            ds = 10;
            rx = (xmax(c)-ds):(xmax(c)+ds);
            ry = (ymax(c)-ds):(ymax(c)+ds);    
            rx(rx<1)=[];
            ry(ry<1)=[];
            rx(rx>size(u,1))=[];
            ry(ry>size(u,2))=[];
            parfor f=1:st
                u_f = u(:,:,c,1,f);
                sample = u_f(rx,ry,1);
                f_data(k,c,f) = mean(sample(:)) - offset;
                f
            end
        end
    end

f_data = squeeze(min(f_data,[],1));
%
colors = {'r','g','b','k','c'};
AXES = handles.time_dependence_corr;
reset(AXES);
LEGEND = cell(0);

polynom_order  = handles.polynom_order;
%
for c=1:sc
    handles.Eb{c} = f_data(c,1);
    %
    intensity = f_data(c,:)';
    frms = (1:st)';
    p = polyfit(frms,intensity,polynom_order);
    intensity_fit = polyval(p,frms);
    f_t = intensity_fit/intensity_fit(1);

    taxis = (frms-1)*handles.min_per_frame/60;    
    semilogy(AXES,taxis,intensity,[colors{c} '.-'],taxis,intensity_fit,'k:','linewidth',3);
    hold(AXES,'on');
    LEGEND = [LEGEND num2str(channel) 'fit'];
    %
    handles.f_t{c} = f_t;
end

guidata(hObject,handles);

hold(AXES,'off');
    xlabel(AXES,'time [h]','fontsize',8);
    ylabel(AXES,'offset-subtracted mean ref. intensity','fontsize',8);
    grid(AXES,'on');
legend(AXES,LEGEND);


% --------------------------------------------------------------------
function calculate_intensity_histograms_for_models_Callback(hObject, eventdata, handles)
            %
            if isempty(handles.raw_img),return,end
            %
            [sx,sy,sc,~,st] = size(handles.raw_img);
            %
            S = round(handles.umppix*8/handles.downsample);
            K = 2.5;
            t = 0.2;
            %
            in = handles.raw_img;
            out = in;
            %
            hw = waitbar(0,'segmenting frames.. ');
            for c=1:sc
                parfor f=1:st
                    u = squeeze(in(:,:,c,1,f));
                    nth = nonlinear_tophat(u,S,K)-1;
                    nth(nth<t)=0;
                    nth = bwmorph(nth,'clean');
                   out(:,:,c,1,f) = nth;
                end
                if ~isempty(hw), waitbar(c/sc,hw); end
            end
            if ~isempty(hw), delete(hw), drawnow; end            
            %icy_imshow(out);
            %
            macc = cell(sc,4);
 
            offset = handles.offset;
            
            if 1==get(handles.compensate_CMHF,'Value')
                f_CMHF_t = handles.f_CMHF_t;
            else
                f_CMHF_t = ones(size(handles.f_CMHF_t));
            end
            
                hw = waitbar(0,'gathering data for correction models.. ');
                for m=1:4                     
                    corrimg = zeros(size(handles.raw_img));                    
                    for c=1:sc 
                            %
                            Eb = handles.Eb{c};
                            p_xy  = handles.p_xy{c};
                            f_t = handles.f_t{c};
                            %
                            parfor f = 1:st
                                I = squeeze(in(:,:,c,1,f)) - offset;
                                I = I/f_CMHF_t(f);
                                EO = [];
                                switch m
                                    case 1 % 'additive (1)'         % f(t) acts only on background
                                        EO = I - Eb*f_t(f)*p_xy;
                                    case 2 % 'multiplicative (1)'   
                                        EO = I./p_xy - Eb*f_t(f);                
                                    case 3 % 'additive (2)'         % f(t) acts on everyhting
                                        EO = I/f_t(f) - Eb*p_xy;                
                                    case 4 % 'multiplicative (2)'   
                                        EO = I./( f_t(f)*p_xy ) - Eb;
                                end
                                EO(EO<0)=0;
                                corrimg(:,:,c,1,f) = EO;
                            end
                    end 
                    %
                    if ~isempty(handles.W)
                        corrimg = introduce_cross_talk_corrections(corrimg,handles);
                    end
                    %                                        
                    parfor c=1:sc                        
                        for f=1:st
                            vals = corrimg(:,:,c,1,f);
                            labs = out(:,:,c,1,f);
                            %
                            % over pixels
                             sample = vals(labs>0);
                             macc{c,m} = [macc{c,m}; sample];
                            % over pixels
                            %
                            % over objects
%                             labs(labs>0)=1;
%                             labs = bwlabel(labs);
%                             for L=1:max(labs(:))    
%                                 s = vals(labs==L);
%                                 macc{c,m} = [macc{c,m}; mean(s(:))];
%                             end
                            % over objects
                        end                        
                    end
                if ~isempty(hw), waitbar(m/4,hw); end                    
                end
                if ~isempty(hw), delete(hw), drawnow; end
                %

colors = {'r','g','b','k','c'};
                for c=1:sc
                    figure('units','normalized','outerposition',[0 0 1 1],'name',[handles.raw_image_filename ', channel: ' num2str(c)]);
                    ax=gca;
                    [N,EDGES] = histcounts(macc{c,1},'Normalization','pdf');
                    plot(ax,EDGES(1:numel(N)),N,[colors{1} '.-'],'linewidth',2);                                                                               
                    hold(ax,'on')
                    [N,EDGES] = histcounts(macc{c,2},'Normalization','pdf');
                    plot(ax,EDGES(1:numel(N)),N,[colors{2} '.-'],'linewidth',2);                                                                               
                    hold(ax,'on')
                    [N,EDGES] = histcounts(macc{c,3},'Normalization','pdf');
                    plot(ax,EDGES(1:numel(N)),N,[colors{3} '.-'],'linewidth',2);                                                                               
                    hold(ax,'on')
                    [N,EDGES] = histcounts(macc{c,4},'Normalization','pdf');
                    plot(ax,EDGES(1:numel(N)),N,[colors{4} '.-'],'linewidth',2);                                                                               
                    hold(ax,'off'); 
                    grid(ax,'on');
                    xlabel(ax,['intensity, channel ' num2str(c)]);
                    ylabel(ax,'pdf');
                    legend(ax,{ ...
                        ['additive (1)' ' kurtosis ' num2str(kurtosis(macc{c,1}-3)) ' std ' num2str(std(macc{c,1}))], ...
                        ['multiplicative (1)' ' kurtosis ' num2str(kurtosis(macc{c,2}-3)) ' std ' num2str(std(macc{c,2}))], ...
                        ['additive (2)' ' kurtosis ' num2str(kurtosis(macc{c,3}-3)) ' std ' num2str(std(macc{c,3}))], ...
                        ['multiplicative (2)' ' kurtosis ' num2str(kurtosis(macc{c,4}-3)) ' std ' num2str(std(macc{c,4}))]});
                end


% --------------------------------------------------------------------
function save_settings_mat_Callback(hObject, eventdata, handles)

            start_dir = get(handles.src_dir,'String');
            [fname, fpath] = uiputfile('*.mat','Save Settings as..',[start_dir filesep 'MicroscopyImageFormatter_settings']);
            if fpath == 0; return; end
            filespec = fullfile(fpath,fname);

            saved_handles.dst_channels = get(handles.dst_channels,'String');
            saved_handles.src_channels = get(handles.src_channels,'String');
            saved_handles.image_type = get(handles.image_type,'Value');
            %
            saved_handles.model_1 = get(handles.model_1,'Value');
            saved_handles.model_2 = get(handles.model_2,'Value');
            saved_handles.model_3 = get(handles.model_3,'Value');
            saved_handles.model_4 = get(handles.model_4,'Value');
            saved_handles.model_5 = get(handles.model_5,'Value');
            
            saved_handles.src_dir = get(handles.src_dir,'String');
            saved_handles.dst_dir = get(handles.dst_dir,'String');
            saved_handles.ref_image_file = get(handles.ref_image_file,'String');
            
            saved_handles.frame_to_show = get(handles.frame_to_show,'Value');
            saved_handles.show_channel = get(handles.show_channel,'Value');
            saved_handles.clean_reference = get(handles.clean_reference,'Value');

            saved_handles.Optosplit_registration_roix = handles.Optosplit_registration_roix;
            saved_handles.Optosplit_registration_roiy1 = handles.Optosplit_registration_roiy1;
            saved_handles.Optosplit_registration_roiy2 = handles.Optosplit_registration_roiy2;
            saved_handles.Optosplit_registration_tform = handles.Optosplit_registration_tform;
            saved_handles.Optosplit_registration_droi_x = handles.Optosplit_registration_droi_x;
            saved_handles.Optosplit_registration_droi_y = handles.Optosplit_registration_droi_y;
            saved_handles.umppix = handles.umppix;
            saved_handles.offset = handles.offset;
            saved_handles.downsample = handles.downsample;
            saved_handles.min_per_frame = handles.min_per_frame;
            saved_handles.polynom_order = handles.polynom_order;
            saved_handles.raw_image_filename = handles.raw_image_filename;
            saved_handles.p_xy = handles.p_xy;
            saved_handles.f_t = handles.f_t;
            saved_handles.Eb = handles.Eb;
            
            saved_handles.compensate_CMHF = get(handles.compensate_CMHF,'Value');
            saved_handles.f_CMHF_t = handles.f_CMHF_t;

            saved_handles.W = handles.W;
            saved_handles.g = handles.g;
            
            saved_handles.Optosplit_registration_a3 = handles.Optosplit_registration_a3;
            saved_handles.Optosplit_registration_f3 = handles.Optosplit_registration_f3; 
            saved_handles.Optosplit_registration_tform3 = handles.Optosplit_registration_tform3;
            saved_handles.Optosplit_registration_roi3x = handles.Optosplit_registration_roi3x;
            saved_handles.Optosplit_registration_roi3y = handles.Optosplit_registration_roi3y;
                                                            
            save(filespec,'saved_handles');
            
% --------------------------------------------------------------------
function load_settings_mat_Callback(hObject, eventdata, handles)
            
           [filename,pathname] = uigetfile({'*.mat','mat files'}, ...
                'Select settings file',get(handles.src_dir,'String'));
            if filename == 0, return, end                        
try            
            load([pathname filesep filename]);

            set(handles.dst_channels,'String',saved_handles.dst_channels);
            set(handles.src_channels,'String',saved_handles.src_channels);
            
            set(handles.image_type,'Value',saved_handles.image_type);

            set(handles.model_1,'Value',saved_handles.model_1);
            set(handles.model_2,'Value',saved_handles.model_2);
            set(handles.model_3,'Value',saved_handles.model_3);            
            set(handles.model_4,'Value',saved_handles.model_4);
            set(handles.model_5,'Value',saved_handles.model_5); 
            model_1_Callback(handles.model_1, eventdata, handles);
            model_2_Callback(handles.model_2, eventdata, handles);
            model_3_Callback(handles.model_3, eventdata, handles);
            model_4_Callback(handles.model_4, eventdata, handles);
            model_5_Callback(handles.model_5, eventdata, handles);            
            %
            set(handles.src_dir,'String',saved_handles.src_dir);
            set(handles.dst_dir,'String',saved_handles.dst_dir);
            set(handles.ref_image_file,'String',saved_handles.ref_image_file);
            
            set(handles.frame_to_show,'Value',saved_handles.frame_to_show);
            set(handles.show_channel,'Value',saved_handles.show_channel);
            set(handles.clean_reference,'Value',saved_handles.clean_reference);           

            handles.Optosplit_registration_roix = saved_handles.Optosplit_registration_roix;
            handles.Optosplit_registration_roiy1 = saved_handles.Optosplit_registration_roiy1;
            handles.Optosplit_registration_roiy2 = saved_handles.Optosplit_registration_roiy2;
            handles.Optosplit_registration_tform = saved_handles.Optosplit_registration_tform;
            handles.Optosplit_registration_droi_x = saved_handles.Optosplit_registration_droi_x;
            handles.Optosplit_registration_droi_y = saved_handles.Optosplit_registration_droi_y;
                handles.umppix = saved_handles.umppix;
                handles.offset = saved_handles.offset;
                handles.downsample = saved_handles.downsample;
                handles.min_per_frame = saved_handles.min_per_frame;
                handles.polynom_order = saved_handles.polynom_order;
            handles.raw_image_filename = saved_handles.raw_image_filename;
            handles.p_xy = saved_handles.p_xy;
            handles.f_t = saved_handles.f_t;
            handles.Eb = saved_handles.Eb;
                        
            set(handles.umppix_edit,'String',num2str(handles.umppix));
            set(handles.offset_edit,'String',num2str(handles.offset));
            set(handles.downsample_edit,'String',num2str(handles.downsample));
            set(handles.min_per_frame_edit,'String',num2str(handles.min_per_frame));
            set(handles.t_dep_fitting_poly_order,'String',num2str(handles.polynom_order));
            
            if ismember(get(handles.image_type,'Value'),[2 3]) % Optosplit
                flag = 'On';
            else
                flag = 'Off';                
            end
                set(handles.setup_Optosplit_registration,'Enable',flag);
                set(handles.setup_Optosplit_registration,'Visible',flag);

            show_corrections_temporal_dependencies(handles); 
            %
            % clean images
            handles.ref_img = [];
            handles.raw_img = [];
            handles.corrected_img = [];
            %
            cla(handles.image_raw,'reset');
            cla(handles.image_corrected,'reset');
                set(handles.image_raw,'XTick',[]);
                set(handles.image_raw,'YTick',[]);
                set(handles.image_corrected,'XTick',[]);
                set(handles.image_corrected,'YTick',[]);            

            set(handles.compensate_CMHF,'Value',saved_handles.compensate_CMHF);
            handles.f_CMHF_t = saved_handles.f_CMHF_t;                
                
            setup_model_controls_visibility(handles);
                
            if isfield(saved_handles,'W') && ~isempty(saved_handles.W)
                handles.W = saved_handles.W;
                handles.g = saved_handles.g;
                cross_talk_correction_flag = 'on';
            else
                handles.W = [];
                handles.g = [];
                cross_talk_correction_flag = 'off';                
            end
                set(handles.cross,'Visible',cross_talk_correction_flag);
                set(handles.talk,'Visible',cross_talk_correction_flag);
                set(handles.correction,'Visible',cross_talk_correction_flag);
                set(handles.unset_spectral_cross_talk_correction,'Enable',cross_talk_correction_flag);
                set(handles.unset_spectral_cross_talk_correction,'Visible',cross_talk_correction_flag);
                        
            guidata(hObject,handles);
            
            handles.Optosplit_registration_a3 = saved_handles.Optosplit_registration_a3;
            handles.Optosplit_registration_f3 = saved_handles.Optosplit_registration_f3; 
            handles.Optosplit_registration_tform3 = saved_handles.Optosplit_registration_tform3;
            handles.Optosplit_registration_roi3x = saved_handles.Optosplit_registration_roi3x;
            handles.Optosplit_registration_roi3y = saved_handles.Optosplit_registration_roi3y;
            
            guidata(hObject,handles);
                                  
catch
     disp('error when loading mat setups');
end

% --------------------------------------------------------------------
function show_corrections_temporal_dependencies(handles)            

if isempty(handles.f_t), return, end

colors = {'r','g','b','k','c'};
AXES = handles.time_dependence_corr;
reset(AXES);

n_channels = numel(handles.Eb);
LEGEND = cell(n_channels,1);

for channel = 1:n_channels

    intensity_fit = handles.f_t{channel};
    head_val = handles.Eb{channel};
    %
    intensity = intensity_fit*head_val;
    %
    st = length(intensity_fit);
    frms = (1:st)';        
    taxis = (frms-1)*handles.min_per_frame/60;    
    semilogy(AXES,taxis,intensity,[colors{channel} '.-'],'linewidth',3);
    hold(AXES,'on');
    LEGEND{channel} = num2str(channel);
           
end

hold(AXES,'off');
    xlabel(AXES,'time [h]','fontsize',8);
    ylabel(AXES,'offset-subtracted mean ref. intensity','fontsize',8);
    grid(AXES,'on');
legend(AXES,LEGEND);


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function setup_model_controls_visibility(handles)

    text_item = {'text12','text13','text14','text15','text16'};
    popup_item = {'model_1','model_2','model_3','model_4','model_5'};
    icon_item = {'icon_1','icon_2','icon_3','icon_4','icon_5'};
    %     
    L = length(get(handles.dst_channels,'String'));
    %
    for m=1:5
            flag = 'On';
        if m>L
            flag = 'Off';
        end
        set(eval(['handles.' text_item{m}]),'Visible',flag);
        set(eval(['handles.' popup_item{m}]),'Visible',flag);
        set(eval(['handles.' popup_item{m}]),'Enable',flag);
        if strcmp('Off',flag)
            imshow(handles.dummy,'Parent',eval(['handles.' icon_item{m}]));
        else
            imshow(handles.model_icon{1},'Parent',eval(['handles.' icon_item{m}]));
        end
    end

% --- Executes on button press in compensate_CMHF.
function compensate_CMHF_Callback(hObject, eventdata, handles)
% hObject    handle to compensate_CMHF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of compensate_CMHF


% --------------------------------------------------------------------
function show_CMHF_correction_on_REF_Callback(hObject, eventdata, handles)

if isempty(handles.ref_img), return, end

offset = handles.offset;

[SX,SY,n_channels,~,st] = size(handles.ref_img);

colors = {'r','g','b','k','c'};
figure;
AXES = gca; 
reset(AXES);
LEGEND = cell(0);

    frms = (1:st)';
    taxis = (frms-1)*handles.min_per_frame/60;    

for channel = 1:n_channels
        
    ref = squeeze(handles.ref_img(:,:,channel,1,:));
    
    intensity = zeros(st,1);

    w = 40; % exclude noisy pixels?
    rx = w:SX-w;
    ry = w:SY-w;

    parfor f=1:st
        intensity(f) = mean(ref(rx,ry,f),'All') - offset;
    end
   
    f_CMHF_t = handles.f_CMHF_t;  
    intensity_CMHF_corrected = intensity./f_CMHF_t;

    semilogy(AXES,taxis,intensity,[colors{channel} '.-'],taxis,intensity(1)*handles.f_t{channel},'k:',taxis,intensity_CMHF_corrected,'m:','linewidth',2);
    hold(AXES,'on');
    LEGEND = [LEGEND num2str(channel) [num2str(channel) ' fit'] [num2str(channel) ' CMHF corrected']];
                   
end

hold(AXES,'off');
    xlabel(AXES,'time [h]','fontsize',8);
    ylabel(AXES,'offset-subtracted mean ref. intensity','fontsize',8);
    grid(AXES,'on');
legend(AXES,LEGEND);


% --------------------------------------------------------------------
function setup_spectral_cross_talk_correction_Callback(hObject, eventdata, handles)

           [filename,pathname] = uigetfile({'*.mat','mat files'}, ...
                'Select spectral cross talk matrix mat file',get(handles.src_dir,'String'));
            if filename == 0, return, end 
try
    load([pathname filesep filename]);
    if ~( exist('W','var') && isnumeric(W) && size(W,1)==size(W,2) && size(W,1)==length(get(handles.dst_channels,'String')) )
        disp('wrong or incompatible object was loaded for cross talk matrix, cannot set it up');
        return, 
    end
    
    if ~( exist('g','var') && isnumeric(g) && length(g)==length(get(handles.dst_channels,'String')) )
        disp('wrong or incompatible object was loaded for g factor, cannot set it up');
        return, 
    end
    if ~iscolumn(g), g=g'; end
    
    if sum(W(:)>=0) ~= size(W,1)*size(W,1)
        disp('entries in the cross talk matrix must be non-negative, cannot set it up');
        return,
    end
    
    % this is done immediately before usage (but, good to remember)
    %     for k=1:size(W,1)
    %         W(:,k) = W(:,k)/norm(W(:,k),1);
    %     end
    
    handles.W = W;
    handles.g = g;
    
    set(handles.cross,'Visible','on');
    set(handles.talk,'Visible','on');
    set(handles.correction,'Visible','on');
    set(handles.unset_spectral_cross_talk_correction,'Enable','on');
    set(handles.unset_spectral_cross_talk_correction,'Visible','on');
    %
    guidata(hObject,handles);
catch 
    disp('error when trying to setup cross-talk matrix');
end

% --------------------------------------------------------------------
function unset_spectral_cross_talk_correction_Callback(hObject, eventdata, handles)
handles.W = []; % cross-talk matrix
handles.g = []; % g factors for channels
set(handles.cross,'Visible','off');
set(handles.talk,'Visible','off');
set(handles.correction,'Visible','off');
set(handles.unset_spectral_cross_talk_correction,'Enable','off');
set(handles.unset_spectral_cross_talk_correction,'Visible','off');
guidata(hObject,handles);

% --------------------------------------------------------------------
function img_corr = introduce_cross_talk_corrections(img,handles)
    W = handles.W;
    g = handles.g;
    img_corr = img;
    if isempty(W), return, end
    %
    [sx,sy,sc,~,st] = size(img);
    %
    % columns of W must sum up to 1
    for c=1:sc
        W(:,c) = W(:,c)/norm(W(:,c),1);
    end
    %
    Winv = eye(sc)/W;
    
    tic
    parfor x=1:sx
        for y=1:sy
            for f=1:st
                v_corr = Winv*( squeeze(img(x,y,:,1,f)).*g ) ;
                v_corr(v_corr<0)=0;
                img_corr(x,y,:,1,f) = v_corr;
            end
        end
        x
    end
    toc/60

% this is about 1.5 times slower :(
%     tic
%     img_corr2 = img;
%     for row=1:sc
%         img_row = zeros(sx,sy,1,1,st);
%         parfor col=1:sc
%             img_row = img_row + img(:,:,col,1,:)*g(col)*Winv(row,col);
%         end
%         img_corr2(:,:,row,1,:) = img_row;
%     end
%     toc/60
    
    
%--------------------------------------------------------------
function [roix,roiy1,roiy2,tform,droi_x,droi_y,tform3,roi3x,roi3y] = get_Optosplit_3_registration_parameters(full_path_to_file,handles)

roix=[];roiy1=[];roiy2=[];tform=[];droi_x=[];droi_y=[];tform3=[];roi3x=[];roi3y=[];
 
[~,~,Ifull] = bfopen_v(full_path_to_file);
[SX,SY,~,~,~] = size(Ifull);

I3 = squeeze(Ifull(:,:,1,1));
I = squeeze(Ifull(:,:,1,2));

[SX,SY] = size(I);

u = squeeze(sum(single(I),5));

D = 10;
YC = round(SY/2);
roiy1 = D:(YC-D);
roiy2 = (YC+D):(SY-D);
roix = D:(SX-D);
u1 = u(roix,roiy1);
    deltax = 0; % to debug 
    deltay = 0;
u2 = u(roix+deltax,roiy2+deltay);

u1 = u1/mean(u1(:));
u2 = u2/mean(u2(:));

     fixed = u1+1;
     warped = u2+1;
     
     fixed_plus = fixed + gsderiv(fixed,2,0);
     warped_plus = warped + gsderiv(warped,2,0);
     
     [optimizer, metric] = imregconfig('multimodal');
     tform = imregtform(warped_plus,fixed_plus,'rigid', optimizer, metric);
     % can visualize at this stage
     registered = imwarp(warped,tform,'OutputView',imref2d(size(fixed)));    
      
     tform.T
    
    % finding diminished ROIs
     z=registered>0;
     S = FindLargestSquares(z);
     [C, H, W, M] = FindLargestRectangles(z,[0 0 1]);
     % icy_imshow(z+M);
     [~, pos] = max(C(:));
     [r, c] = ind2sub(size(S), pos);
     droi_y = c:W(r,c);
     droi_x = r:H(r,c);
     
%      fixed = fixed(droi_x,droi_y);
%      registered = registered(droi_x,droi_y);      
%      iv = zeros(size(fixed,1),size(fixed,2),2,1,1);            
%      iv(:,:,1,1,1) = fixed;
%      iv(:,:,2,1,1) = registered;        
%      icy_imshow(uint16(iv));

% 3-rd channel
     u3 = double(I3);
        %      f3 = 1/1.8; 
        %      a3 = 270;
        a3 = handles.Optosplit_registration_a3;
        f3 = handles.Optosplit_registration_f3; % essentially, hardcoded
     u3 = imresize(u3,f3);
     u3 = imrotate(u3,a3);
%     icy_imshow(fixed);
%     icy_imshow(u3);
               
     warped = u3;
     
     warped = warped+1;
     fixed = fixed+1;
     
     fixed_plus = fixed; % + gsderiv(fixed,2,0);
     warped_plus = warped; % + gsderiv(warped,2,0);
             
     [optimizer, metric] = imregconfig('multimodal');     
%         optimizer.InitialRadius = 0.009;
         optimizer.Epsilon = 1.5e-4;
         %optimizer.GrowthFactor = 1.01;
         optimizer.MaximumIterations = 600;     
          
     tform3 = imregtform(warped_plus,fixed_plus,'affine', optimizer, metric);
     % can visualize at this stage
     registered = imwarp(warped,tform3,'OutputView',imref2d(size(fixed)));    
     
     z=registered>0;
     S = FindLargestSquares(z);
     [C, H, W, M] = FindLargestRectangles(z,[0 0 1]);
     % icy_imshow(z+M);
     [~, pos] = max(C(:));
     [r, c] = ind2sub(size(S), pos);
     roi3y = c:W(r,c);
     roi3x = r:H(r,c);
     
     fixed = fixed(roi3x,roi3y);
     registered = registered(roi3x,roi3y);          
     iv = zeros(size(fixed,1),size(fixed,2),2,1,1);            
     iv(:,:,1,1,1) = fixed;
     iv(:,:,2,1,1) = registered;        
     %icy_imshow(uint16(iv));    
     figure;imagesc(uint8(map(fixed + 2*registered,0,255)));daspect(gca,[1 1 1]);


%--------------------------------------------------------------
function v = load_Optosplit_image_3(handles,full_path_to_image)

% v = [];
% 
roix = handles.Optosplit_registration_roix;
roiy1 = handles.Optosplit_registration_roiy1;
roiy2 = handles.Optosplit_registration_roiy2;
tform = handles.Optosplit_registration_tform;
droi_x = handles.Optosplit_registration_droi_x;
droi_y = handles.Optosplit_registration_droi_y;
tform3 = handles.Optosplit_registration_tform3; 
roi3x = handles.Optosplit_registration_roi3x; 
roi3y = handles.Optosplit_registration_roi3y;
a3 = handles.Optosplit_registration_a3; 
f3 = handles.Optosplit_registration_f3;

if isempty(roix) || isempty(tform3)
    disp('cannot load Optosplit image - registration not set up');
    return;
end

[~,~,Ifull] = bfopen_v(full_path_to_image);
[SX,SY,sc,sz,st3] = size(Ifull);

if sc~=1
    I3 = Ifull(:,:,1,:,:);
    I =  Ifull(:,:,2,:,:);
elseif sz~=1
    I3 = Ifull(:,:,:,1,:);
    I =  Ifull(:,:,:,2,:);    
else
    I3 = Ifull(:,:,:,:,1:2:st3);
    I = Ifull(:,:,:,:,1+(1:2:st3));
end
    [SX,SY,~,~,st] = size(I);
%icy_imshow(I3);
%icy_imshow(I);

        s = get(handles.dst_channels,'String');
            nch1 = int64(str2num(s(1)));
            nch2 = int64(str2num(s(2)));
            nch3 = int64(str2num(s(3)));            
    %
    if length(s)~=3
        s = s(1:3);
        set(handles.dst_channels,'String',s);
    end
    %
    v = [];
    %
    for f=1:st
        u = single(squeeze(I(:,:,1,1,f)));
        u1 = u(roix,roiy1);
        u2 = u(roix,roiy2);
        
        u2_reg = imwarp(u2,tform,'OutputView',imref2d(size(u1)));
        
        u1 = u1(droi_x,droi_y);
        u2 = u2_reg(droi_x,droi_y);
        
        u3 = single(squeeze((I3(:,:,1,1,f))));
        u3 = imresize(u3,f3);
        u3 = imrotate(u3,a3);        
        u3_reg = imwarp(u3,tform3,'OutputView',imref2d(size(u1)));
        u1 = u1(roi3x,roi3y);
        u2 = u2(roi3x,roi3y);
        u3 = u3_reg(roi3x,roi3y);
        
        if isempty(v)
            v = zeros(size(u1,1),size(u1,2),3,1,st);
        end
        v(:,:,nch1,1,f) = u1;
        v(:,:,nch2,1,f) = u2;
        v(:,:,nch3,1,f) = u3;
    end
    %
    
    
    