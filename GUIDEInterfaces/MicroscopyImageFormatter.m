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

% Last Modified by GUIDE v2.5 16-Mar-2021 19:39:00

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

set(handles.image_type,'String',{'Nikon','Optosplit 2 channels','Optosplit 3 channels'});

handles.Optosplit_registration_roix = [];
handles.Optosplit_registration_roiy1 = [];
handles.Optosplit_registration_roiy2 = [];
handles.Optosplit_registration_tform = [];
handles.Optosplit_registration_droi_x = [];
handles.Optosplit_registration_droi_y = [];

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
    handles.polynom_order = 12;

set(handles.umppix_edit,'String',num2str(handles.umppix));
set(handles.offset_edit,'String',num2str(handles.offset));
set(handles.downsample_edit,'String',num2str(handles.downsample));
set(handles.min_per_frame_edit,'String',num2str(handles.min_per_frame));
set(handles.src_channels,'String',src_channels);
set(handles.dst_channels,'String',dst_channels);
set(handles.t_dep_fitting_poly_order,'String',num2str(handles.polynom_order));
set(handles.clean_reference,'Value',1);

set(handles.show_channel,'String',{'1','2','3','4','5','All'});

handles.ref_img = [];
handles.raw_img = [];
handles.corrected_img = [];

handles.raw_image_filename = [];

handles.p_xy = cell(0);
handles.f_t = cell(0);
handles.Eb = cell(0);

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
% hObject    handle to model_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns model_1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from model_1


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
% hObject    handle to model_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns model_2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from model_2


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
% hObject    handle to model_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns model_3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from model_3


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
% hObject    handle to model_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns model_4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from model_4


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
% hObject    handle to model_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns model_5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from model_5


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
    hw = waitbar(0,['introducing corrections to channel ' num2str(c)]);
    for f = 1:st
        I = squeeze(handles.raw_img(:,:,c,1,f)) - handles.offset;
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
    %    
end
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
            hw = waitbar(0,['introducing corrections to channel ' num2str(c)]);
            for f = 1:st
                I = squeeze(handles.raw_img(:,:,c,1,f)) - handles.offset;
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
            v = [];  % to do 
        case 'Nikon'
            v = load_Nikon_image(handles,full_path_to_file);
    end


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
        
        img = single(image(:,:,current_channel_to_show,1,frame));
        %
        t = mean(img(:)) + 1.65*std(img(:));
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
if isempty(directoryname), return, end
if isfolder(directoryname) 
    set(handles.src_dir,'String',directoryname);
    guidata(hObject,handles);
end

% --- Executes on button press in set_dst_dir.
function set_dst_dir_Callback(hObject, eventdata, handles)
directoryname = uigetdir(pwd,'Pick DST Directory');
if isempty(directoryname), return, end
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
    if isempty(handles.raw_img),return, end
    try
        icy_imshow(uint16(handles.raw_img));
    catch
        disp('cannot send image to Icy');
    end


% --- Executes on button press in send_corrected_to_Icy.
function send_corrected_to_Icy_Callback(hObject, eventdata, handles)
    if isempty(handles.corrected_img),return, end
    try
        icy_imshow(uint16(handles.corrected_img));
    catch
        disp('cannot send image to Icy');
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

guidata(hObject,handles);

hold(AXES,'off');
    xlabel(AXES,'time [h]','fontsize',8);
    ylabel(AXES,'offset-subtracted mean ref. intensity','fontsize',8);
    grid(AXES,'on');
legend(AXES,LEGEND);


% --- Executes on selection change in image_type.
function image_type_Callback(hObject, eventdata, handles)
    if 2==get(hObject,'Value') % Optosplit
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
function v = load_Nikon_image(handles,full_path_to_file)
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

[SX,SY,~,~,~] = size(I);

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
%     icy_imshow(uint16(v),['dx = ' num2str(dx) ', dy = ' num2str(dy)]);


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
end

img_data = [];
for k=1:numel(filenames)
    img_data = cat(5,img_data,img_acc{k});
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
            ds = 10;
            rx = (xmax(c)-ds):(xmax(c)+ds);
            ry = (ymax(c)-ds):(ymax(c)+ds);    
            rx(rx<1)=[];    
            ry(ry<1)=[];           
            %            
            u = img_acc{k};
            parfor f=1:st
                u_f = u(:,:,c,1,f);
                sample = u_f(rx,ry,1);
                f_data(k,c,f) = median(sample(:)) - offset;
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
            S = round(handles.umppix*8);
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
            
            save(filespec,'saved_handles');
            
% --------------------------------------------------------------------
function load_settings_mat_Callback(hObject, eventdata, handles)
            
           [filename,pathname] = uigetfile({'*.mat','mat files'}, ...
                'Select settings file',get(handles.src_dir,'String'));
            if filename == 0, return, end                        
%try            
            load([pathname filesep filename]);

            set(handles.dst_channels,'String',saved_handles.dst_channels);
            set(handles.src_channels,'String',saved_handles.src_channels);
            
            set(handles.image_type,'Value',saved_handles.image_type);

            set(handles.model_1,'Value',saved_handles.model_1);
            set(handles.model_2,'Value',saved_handles.model_2);
            set(handles.model_3,'Value',saved_handles.model_3);            
            set(handles.model_4,'Value',saved_handles.model_4);
            set(handles.model_5,'Value',saved_handles.model_5);            
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
            
            if 2 == get(handles.image_type,'Value')
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
            
            guidata(hObject,handles);
                                  
% catch
%      disp('error when loading mat setups');
% end

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
