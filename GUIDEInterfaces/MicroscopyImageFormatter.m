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

% Last Modified by GUIDE v2.5 05-Mar-2021 09:45:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
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

ALYtools_dir = 'C:\Users\alexany\ALYtools';
addpath(ALYtools_dir);
addpath_ALYtools;

this_dir = 'X:\working\fogim\users\Yuriy Alexandrov\For_Jon_Jan_27_2021';

% HARDCODED
set(handles.src_dir,'String',[this_dir filesep 'Alm236_P4_1']);
set(handles.dst_dir,'String',[this_dir filesep 'Alm236_P4_1_formatted']);
set(handles.ref_image_file,'String',[this_dir filesep 'Alm236_P4_1' filesep 'ALM236_P4_1_MMStack_C-9 - added as no.00029.ome.tif']);

set(handles.show_channel,'String',{'1','2','3','4','5','All'});

set(handles.umppix_edit,'String','0.68');
set(handles.offset_edit,'String','100');
set(handles.downsample_edit,'String','1');
set(handles.min_per_frame_edit,'String','5');
% values
handles.umppix = 0.68';
handles.offset = 100';
handles.downsample = 1;
handles.min_per_frame = 5;

handles.ref_img = [];   
handles.raw_img = [];   
handles.corrected_img = [];  

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
% hObject    handle to umppix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of umppix_edit as text
%        str2double(get(hObject,'String')) returns contents of umppix_edit as a double


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
% hObject    handle to min_per_frame_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_per_frame_edit as text
%        str2double(get(hObject,'String')) returns contents of min_per_frame_edit as a double


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
% hObject    handle to downsample_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of downsample_edit as text
%        str2double(get(hObject,'String')) returns contents of downsample_edit as a double


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
% hObject    handle to src_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of src_channels as text
%        str2double(get(hObject,'String')) returns contents of src_channels as a double


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
% hObject    handle to dst_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dst_channels as text
%        str2double(get(hObject,'String')) returns contents of dst_channels as a double


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
% hObject    handle to offset_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offset_edit as text
%        str2double(get(hObject,'String')) returns contents of offset_edit as a double


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

            v = load_Optosplit_image(full_path_to_file);
            if isempty(v), return, end
            %
            handles.raw_img = v;
            guidata(hObject,handles);          
            %
            [~,FNAME,FEXT] = fileparts(full_path_to_file);
            show_image(handles,'raw_img','image_raw',[FNAME FEXT]);
            show_image(handles,'raw_img','image_corrected',[]);                
        
% --- Executes on button press in recalculate_corrections.
function recalculate_corrections_Callback(hObject, eventdata, handles)
% hObject    handle to recalculate_corrections (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function frame_to_show_Callback(hObject, eventdata, handles)
% hObject    handle to frame_to_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_to_show as text
%        str2double(get(hObject,'String')) returns contents of frame_to_show as a double


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
% hObject    handle to show_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns show_channel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from show_channel


% --- Executes during object creation, after setting all properties.
function show_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to show_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t_dep_fitting_poly_order_Callback(hObject, eventdata, handles)
% hObject    handle to t_dep_fitting_poly_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_dep_fitting_poly_order as text
%        str2double(get(hObject,'String')) returns contents of t_dep_fitting_poly_order as a double


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


% --- Executes on button press in do_them_all.
function do_them_all_Callback(hObject, eventdata, handles)
% hObject    handle to do_them_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in load_ref.
function load_ref_Callback(hObject, eventdata, handles)

    full_path_to_file = get(handles.ref_image_file,'String');

        v = load_Optosplit_image(full_path_to_file);
        if isempty(v), return, end
        handles.ref_img = v;
        handles.raw_img = v;
        guidata(hObject,handles);
        %
        [~,FNAME,FEXT] = fileparts(full_path_to_file);
        show_image(handles,'raw_img','image_raw',[FNAME FEXT]);
        show_image(handles,'raw_img','image_corrected',[]);

% ---         
function show_image(handles,what_to_show,where_to_show,TITLE)    

    image = eval(['handles.' what_to_show]);
    ax = eval(['handles.' where_to_show]);

        frame_to_show = str2num(get(handles.frame_to_show,'String'));

        s = get(handles.show_channel,'String');
        ind = get(handles.show_channel,'Value');
        current_channel_to_show = str2num(s{ind});

        img = image(:,:,current_channel_to_show,1,frame_to_show);
        imshow(uint8(map(img,0,255)), 'Parent', ax);
        
        if ~isempty(TITLE)
            title(ax,TITLE);
        end
    
% --- Example function used for debugging with Optosplit data    
function v = load_Optosplit_image(full_path_to_image)

v = [];

    Y_split = 342;
    dx = 8; 
    dy = 8; 
    Lx = 492+4;
    Ly = 326+4;
    rxD = dx:dx+Lx;
    ryD = dy:dy+Ly;
    ryA = Y_split+dy:Y_split+dy+Ly;
    rxA = dx:dx+Lx;

load('tform_D_to_A_P4_1','tform');    

try
    [~,~,I] = bfopen_v(full_path_to_image);
    [sx,sy,sc,sz,st] = size(I);

    %I(:,:,:,:,120:122) = I(:,:,:,:,117:119);
    %I(:,:,:,:,123:126) = I(:,:,:,:,127:130);
    % icy_imshow(I); % it is OK

    v = [];

    for f=1:st
        uD = single(I(rxD,ryD,1,1,f));
        uA = single(I(rxA,ryA,1,1,f));
        [SX,SY]=size(uD);
        uA_reg = imwarp(uA,tform,'OutputView',imref2d(size(uD)));
        %
        uD = uD(2:SX-1,5:SY-1);
        uA_reg = uA_reg(2:SX-1,5:SY-1);

        if isempty(v)
            v = zeros(size(uD,1),size(uD,2),2,1,st);
        end
        v(:,:,1,1,f) = uD;
        v(:,:,2,1,f) = uA_reg;    
    end

catch
    disp(['error when trying to load image ' full_path_to_image]);
end


% --- Executes on button press in set_src_dir.
function set_src_dir_Callback(hObject, eventdata, handles)
% hObject    handle to set_src_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in set_dst_dir.
function set_dst_dir_Callback(hObject, eventdata, handles)
% hObject    handle to set_dst_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in set_ref_image_file.
function set_ref_image_file_Callback(hObject, eventdata, handles)
% hObject    handle to set_ref_image_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in send_raw_to_Icy.
function send_raw_to_Icy_Callback(hObject, eventdata, handles)
    try
        icy_imshow(handles.raw_img);
    catch
        disp('cannot send image to Icy');
    end


% --- Executes on button press in send_corrected_to_Icy.
function send_corrected_to_Icy_Callback(hObject, eventdata, handles)
    try
        icy_imshow(handles.corrected_img);
    catch
        disp('cannot send image to Icy');
    end
