function varargout = t_dependent_Nuclei_ratio_FRET_ratio_heatmapper(varargin)
% T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER MATLAB code for t_dependent_Nuclei_ratio_FRET_ratio_heatmapper.fig
%      T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER, by itself, creates a new T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER or raises the existing
%      singleton*.
%
%      H = T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER returns the handle to a new T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER or the handle to
%      the existing singleton*.
%
%      T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER.M with the given input arguments.
%
%      T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER('Property','Value',...) creates a new T_DEPENDENT_NUCLEI_RATIO_FRET_RATIO_HEATMAPPER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before t_dependent_Nuclei_ratio_FRET_ratio_heatmapper_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to t_dependent_Nuclei_ratio_FRET_ratio_heatmapper_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help t_dependent_Nuclei_ratio_FRET_ratio_heatmapper

% Last Modified by GUIDE v2.5 01-Nov-2018 16:00:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @t_dependent_Nuclei_ratio_FRET_ratio_heatmapper_OpeningFcn, ...
                   'gui_OutputFcn',  @t_dependent_Nuclei_ratio_FRET_ratio_heatmapper_OutputFcn, ...
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


% --- Executes just before t_dependent_Nuclei_ratio_FRET_ratio_heatmapper is made visible.
function t_dependent_Nuclei_ratio_FRET_ratio_heatmapper_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to t_dependent_Nuclei_ratio_FRET_ratio_heatmapper (see VARARGIN)

% Choose default command line output for t_dependent_Nuclei_ratio_FRET_ratio_heatmapper
handles.output = hObject;

handles.TrackPlotter_handles = varargin{1};

[~,NAME,EXT]=fileparts(handles.TrackPlotter_handles.fullfilename);
set(handles.figure1, 'Name', ['FRET ratio heatmap : ' NAME EXT]);

set(handles.features_chooser,'String',handles.TrackPlotter_handles.features);

set(handles.features_chooser,'Value',6);

% Update handles structure
guidata(hObject, handles);

visualize(handles);

% UIWAIT makes t_dependent_Nuclei_ratio_FRET_ratio_heatmapper wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = t_dependent_Nuclei_ratio_FRET_ratio_heatmapper_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in features_chooser.
function features_chooser_Callback(hObject, eventdata, handles)
% hObject    handle to features_chooser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns features_chooser contents as cell array
%        contents{get(hObject,'Value')} returns selected item from features_chooser
visualize(handles);

% --- Executes during object creation, after setting all properties.
function features_chooser_CreateFcn(hObject, eventdata, handles)
% hObject    handle to features_chooser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ------------------------------------------------------
function visualize(handles)
    tracks = handles.TrackPlotter_handles.MI_tracks;
    display_tracks = handles.TrackPlotter_handles.MI_norm_FRET_ratio;
        
    par_ind = get(handles.features_chooser,'Value');        
              index = 0;
                 switch par_ind                    
                    case 4 % #nghbrs                    
                        if ismember(size(tracks,2),[10 11 15])
                            index = 9;
                        end
                    case 5 % cell density
                        if ismember(size(tracks,2),[10 11 15])
                            index = 10;
                        end
                    case 6 % FRET ratio
                            index = 4;
                    case 8 % donor intensity                     
                        index = 5;
                    case 9 % acceptor intensity
                        index = 6;
                    case 10 % nucleus size
                        index = 7;
                    case 11 % Pearson                      
                        index = 8;
                    case 14 % FRET molar fraction
                        index = 11;
                    case 15 % cell size
                        if 15== size(tracks,2)
                        index = 12;
                        end
                    case 16 % 
                       if 15== size(tracks,2)
                        index = 13;
                        end
                    case 17 % 
                        if 15== size(tracks,2)
                        index = 14;
                        end
                    case 18 % 
                         if 15== size(tracks,2)
                        index = 15;
                        end
                 end                  
   
                 if 0==index, return, end
                 
                 sortarr = [];
                 for k=1:numel(tracks)
                    track = tracks{k};
                    vals = squeeze(track(:,index));
                    sortarr = [sortarr; median(vals)];
                 end
                   
     [~,I]=sort(sortarr,'descend');
     
     M = display_tracks(I,:);
        
%     import bioma.data.*
%     M=map(M,0,1);
%     DMobj = DataMatrix(M);
%     hmo = HeatMap(DMobj,'Colormap',redgreencmap);
imagesc(handles.heatmap_axes,M);

    cmap = jet(256);       
    cmap(1,:)=[0,0,0];
    colormap(handles.heatmap_axes,cmap);

xlabel(handles.heatmap_axes,'time [h]');
str = get(handles.features_chooser,'String'); 
ylabel(handles.heatmap_axes,str{par_ind});
tax =  handles.TrackPlotter_handles.dt*(1:size(M,2));
tax = 0.1*fix(tax*10);
sortarr = 0.1*fix(sortarr*10);
set(handles.heatmap_axes, 'xticklabel', {linspace(min(tax),max(tax),11)}, ... 
    'yticklabel', {flip(linspace(min(sortarr),max(sortarr),11))});
grid on

% --- Executes on button press in selection_only.
function selection_only_Callback(hObject, eventdata, handles)
% hObject    handle to selection_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of selection_only
visualize(handles);