function varargout = t_dependent_Nuclei_ratio_FRET_visualizer(varargin)
% T_DEPENDENT_NUCLEI_RATIO_FRET_VISUALIZER MATLAB code for t_dependent_Nuclei_ratio_FRET_visualizer.fig
%      T_DEPENDENT_NUCLEI_RATIO_FRET_VISUALIZER, by itself, creates a new T_DEPENDENT_NUCLEI_RATIO_FRET_VISUALIZER or raises the existing
%      singleton*.
%
%      H = T_DEPENDENT_NUCLEI_RATIO_FRET_VISUALIZER returns the handle to a new T_DEPENDENT_NUCLEI_RATIO_FRET_VISUALIZER or the handle to
%      the existing singleton*.
%
%      T_DEPENDENT_NUCLEI_RATIO_FRET_VISUALIZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in T_DEPENDENT_NUCLEI_RATIO_FRET_VISUALIZER.M with the given input arguments.
%
%      T_DEPENDENT_NUCLEI_RATIO_FRET_VISUALIZER('Property','Value',...) creates a new T_DEPENDENT_NUCLEI_RATIO_FRET_VISUALIZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before t_dependent_Nuclei_ratio_FRET_visualizer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to t_dependent_Nuclei_ratio_FRET_visualizer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help t_dependent_Nuclei_ratio_FRET_visualizer

% Last Modified by GUIDE v2.5 23-Sep-2018 18:15:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;

gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @t_dependent_Nuclei_ratio_FRET_visualizer_OpeningFcn, ...
                   'gui_OutputFcn',  @t_dependent_Nuclei_ratio_FRET_visualizer_OutputFcn, ...
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


% --- Executes just before t_dependent_Nuclei_ratio_FRET_visualizer is made visible.
function t_dependent_Nuclei_ratio_FRET_visualizer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to t_dependent_Nuclei_ratio_FRET_visualizer (see VARARGIN)

% Choose default command line output for t_dependent_Nuclei_ratio_FRET_visualizer
handles.output = hObject;

handles.TrackPlotter_handles = varargin{1};

[~,NAME,EXT]=fileparts(handles.TrackPlotter_handles.fullfilename);
set(handles.figure1, 'Name', ['vidi : ' NAME EXT]);
    
minmaxlimits = get(handles.TrackPlotter_handles.filter_table,'Data');
set(handles.filter_table,'Data', minmaxlimits);
set(handles.filter_table, 'RowName', [{'time'} handles.TrackPlotter_handles.features]);
set(handles.filter_table, 'ColumnName', {'min','max'});
set(handles.filter_table,'ColumnEditable',false);

set(handles.visual_mode,'String',{'FRET ratio','selection'});

set(handles.filter_table, 'ColumnFormat', {'numeric','numeric'} );

%
% pos_7_post.OME.tif_analysis_output.OME.tiff
% pos_7_post_FRET_ratio_featured_TRACKMATE_OUTPUT.mat
%
tiffname  = strrep(handles.TrackPlotter_handles.fullfilename, ...
            '_FRET_ratio_featured_TRACKMATE_OUTPUT.mat', ...
            '.OME.tif_analysis_output.OME.tiff');
        
try
    hw = waitbar(0,'Loading image file..');
    if ~isempty(hw), waitbar(0.1,hw); drawnow, end; 
    [~,~,I] = bfopen_v(tiffname);
    if ~isempty(hw), delete(hw), drawnow; end;
    
    handles.image = squeeze(I(:,:,1,3,:));

    set(handles.theslider,'Min',1);
    set(handles.theslider,'Max',size(handles.image,3));
    set(handles.theslider,'Value',1);
    
    A = I(:,:,1,3,:);
    lowestValue = min(A(A(:)>0));  
    highestValue = quantile(A(A(:)~=0),0.9);
    cmap = jet(256);
    caxis(handles.image_axes,[lowestValue-2/256, highestValue]);
    cmap(1,:)=[0,0,0];
    colormap(handles.image_axes,cmap);    
    %
    set(handles.image_axes,'XTick',[]);
    set(handles.image_axes,'YTick',[]);
    %
    c = colorbar(handles.image_axes,'TickLabels',{linspace(double(lowestValue),double(highestValue),11)/100});    
    c.Label.String = 'FRET ratio';  
    set(c,'fontsize',24);
    handles.lowestValue = lowestValue;
    handles.highestValue = highestValue;
    %
    handles.image_mask = calculate_image_mask(hObject, eventdata, handles);
    
    % update_button handles structure
    guidata(hObject, handles);    
    %            
catch
     if exist('hw','var') && ~isempty(hw), delete(hw), drawnow; end;
     fh = ancestor(hObject,'figure');     
     delete(fh);
end

%show_frame(hObject, eventdata, handles);

% UIWAIT makes t_dependent_Nuclei_ratio_FRET_visualizer wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = t_dependent_Nuclei_ratio_FRET_visualizer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function theslider_Callback(hObject, eventdata, handles)
% hObject    handle to theslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
show_frame(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function theslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%--------------------------------------------------------------------------------------%
function show_frame(hObject, eventdata, handles) 
    index = round(get(handles.theslider,'Value'));
    AXES = handles.image_axes;
    
    if 1==get(handles.visual_mode,'Value')
        imagesc(double(handles.image(:,:,index)),'Parent',AXES);
    else        
        u = double(handles.image(:,:,index)>0) + double(handles.image_mask(:,:,index));
        v_min = double(handles.lowestValue);
        v_max = 0.95*double(handles.highestValue);
        u = v_min + (v_max-v_min)*u/3;
        imagesc(u,'Parent',AXES);
    end
    
    caxis(AXES,[handles.lowestValue-2/256, handles.highestValue]);
    daspect(AXES,[1 1 1]);
    set(AXES,'XTick',[]);
    set(AXES,'YTick',[]);
    
%--------------------------------------------------------------------------------------%    
function image_mask = calculate_image_mask(hObject, eventdata, handles)
    %
    DR = handles.TrackPlotter_handles.raw_data;
    D = handles.TrackPlotter_handles.track_data;
    mask = handles.TrackPlotter_handles.mask;
    %
    hw = waitbar(0,'Labeling..');    
    %
    L = zeros(size(handles.image));
    for f=1:size(L,3) % over frames
        L(:,:,f) = bwlabel(squeeze(handles.image(:,:,f))>0);
        if ~isempty(hw), waitbar(f/size(L,3),hw); drawnow, end;
    end
    %
    N = size(D,1);    
    image_mask = zeros(size(handles.image),'uint8');
    FRAMES_notOK = [];
    LABELS_notOK = [];
    FRAMES_OK = [];
    LABELS_OK = [];    
    %
    waitbar(0,hw,'Matching..');
    %
    for k=1:N
        if ~isempty(hw), waitbar(k/N,hw); drawnow, end;
            track = DR{k};
            frames = round(track(:,1));                
            x = round(track(:,2));
            y = round(track(:,3));
                x(0==x)=1;
                y(0==y)=1;
            labels = zeros(size(frames));
            for m=1:numel(frames)
                labels(m) = L(x(m),y(m),frames(m));
            end 
        if 1==mask(k) % accepted            
            FRAMES_OK = cat(1,FRAMES_OK,frames);
            LABELS_OK = cat(1,LABELS_OK,labels);
        else % rejected
            FRAMES_notOK = cat(1,FRAMES_notOK,frames);
            LABELS_notOK = cat(1,LABELS_notOK,labels);            
        end
    end
    %
    waitbar(0,hw,'Categorizing..')
    %
    for f=1:size(L,3)
        if ~isempty(hw), waitbar(f/size(L,3),hw); drawnow, end;
        cur_labels = LABELS_OK(FRAMES_OK==f);
        cur_labels = cur_labels(cur_labels~=0);
        cur_mask_OK = ismember(squeeze(L(:,:,f)),cur_labels);
        %
        cur_labels = LABELS_notOK(FRAMES_notOK==f);
        cur_labels = cur_labels(cur_labels~=0);
        cur_mask_notOK = ismember(squeeze(L(:,:,f)),cur_labels);        
        image_mask(:,:,f) = 2*cur_mask_OK + cur_mask_notOK;
    end
    if ~isempty(hw), delete(hw), drawnow; end;
    guidata(hObject, handles);
    show_frame(hObject, eventdata, handles);
          
% --- Executes on selection change in visual_mode.
function visual_mode_Callback(hObject, eventdata, handles)
% hObject    handle to visual_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns visual_mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from visual_mode
show_frame(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function visual_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to visual_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in update_button.
function update_button_Callback(hObject, eventdata, handles)
% hObject    handle to update_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fullfilename = handles.TrackPlotter_handles.fullfilename;

FIGS = findobj(0, 'type', 'figure');
for k=1:size(FIGS,1)
    h = FIGS(k,1);
    cur_handles = guidata(h);
    try
    if strcmp(cur_handles.fullfilename,fullfilename) % mine
        handles.TrackPlotter_handles = cur_handles;
        guidata(hObject, handles); 
        minmaxlimits = get(handles.TrackPlotter_handles.filter_table,'Data');
        set(handles.filter_table,'Data', minmaxlimits);
        handles.image_mask = calculate_image_mask(hObject, eventdata, handles);
        guidata(hObject, handles);
        show_frame(hObject, eventdata, handles);
        break;
    end
    catch
    end
end


% --------------------------------------------------------------------
function image_zoom_OffCallback(hObject, eventdata, handles)
% hObject    handle to image_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hObject
