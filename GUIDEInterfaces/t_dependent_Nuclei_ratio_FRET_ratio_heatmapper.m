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

% Last Modified by GUIDE v2.5 21-Apr-2020 14:11:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
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
set(handles.curve_type,'String',{'mean','median'});

% Update handles structure
guidata(hObject, handles);

visualize(hObject,handles);

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
visualize(hObject,handles);

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
function visualize(hObject,handles)

    % to update mask
    fullfilename = handles.TrackPlotter_handles.fullfilename;
    FIGS = findobj(0, 'type', 'figure');
    for k=1:size(FIGS,1)
        h = FIGS(k,1);
        cur_handles = guidata(h);
        try
        if strcmp(cur_handles.fullfilename,fullfilename) % mine
            handles.TrackPlotter_handles = cur_handles;
            guidata(hObject, handles); 
            break;
        end
        catch
        end
    end
    % to update mask

    tracks = handles.TrackPlotter_handles.MI_tracks;
    display_tracks = handles.TrackPlotter_handles.MI_norm_FRET_ratio;
    
    if 1==get(handles.selection_only,'Value') % if masking option is on
        mask = handles.TrackPlotter_handles.mask;
        tracks = tracks(mask);
        display_tracks = display_tracks(mask,:);
    end
    
    track = tracks{1};
    
    par_ind = get(handles.features_chooser,'Value');        
              index = 0;
                 switch par_ind                    
                    case 4 % #nghbrs                    
                        if ismember(size(track,2),[10 11 15])
                            index = 9;
                        end
                    case 5 % cell density
                        if ismember(size(track,2),[10 11 15])
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
                        if 15== size(track,2)
                        index = 12;
                        end
                    case 16 % 
                       if 15== size(track,2)
                        index = 13;
                        end
                    case 17 % 
                        if 15== size(track,2)
                        index = 14;
                        end
                    case 18 % 
                         if 15== size(track,2)
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
                   
     [~,I] = sort(sortarr,'descend');
     
     M = display_tracks(I,:);
        
%     import bioma.data.*
%     M=map(M,0,1);
%     DMobj = DataMatrix(M);
%     hmo = HeatMap(DMobj,'Colormap',redgreencmap);
imagesc(handles.heatmap_axes,M);

    cmap = jet(256);       
    cmap(1,:)=[0,0,0];
    colormap(handles.heatmap_axes,cmap);

%%%% 2 default curves
   cla(handles.curves_axes,'reset');
   N = size(M,1);
if true && N>=2
    
   m1 = M(1:fix(N/2),:);
   m2 = M(fix(N/2):N,:);
   
   mode_ind = get(handles.curve_type,'Value');
   mode_str = get(handles.curve_type,'String');
   if 1 == mode_ind % 'mean'
       curve1 = mean(m1,1);
       curve2 = mean(m2,1);
   else
       curve1 = median(m1,1);
       curve2 = median(m2,1);       
   end
   
t = round((0:length(curve1)-1)*handles.TrackPlotter_handles.dt*60);
plot(handles.curves_axes,t,curve1,'k--',t,curve2,'k:','linewidth',2);
legend(handles.curves_axes,{['upper half: ' mode_str{mode_ind}],['lower half: ' mode_str{mode_ind}]});
xlabel(handles.curves_axes,'time [min]');
set(handles.curves_axes,'yticklabel', []);
axis(handles.curves_axes,[min(t) max(t) min(min(curve1(:)),min(curve2(:))) max(max(curve1(:)),max(curve2(:)))]);
grid(handles.curves_axes,'on'); 

str = get(handles.features_chooser,'String'); 
ylabel(handles.heatmap_axes,str{par_ind});
%
sarr = 0.1*round(sortarr*10);
b1=min(sarr(:));
b2=max(sarr(:));
step = min(length(sarr),fix(length(sarr)/10));
yticks(handles.heatmap_axes,0:step:(length(sarr)));
step_val = (b2 - b1)/10;
yticklabels(handles.heatmap_axes,flip(b1+step_val*(0:10)));
%
b1=min(t);
b2=max(t);
step = fix(length(t)/4);
xticks(handles.heatmap_axes,0:step:(length(t)));
step_val = round((b2 - b1)/4);
xticklabels(handles.heatmap_axes,b1+step_val*(0:4));
%
xlabel(handles.heatmap_axes,'time [min]');
grid(handles.heatmap_axes,'on'); 
    
end






% --- Executes on button press in selection_only.
function selection_only_Callback(hObject, eventdata, handles)
% hObject    handle to selection_only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of selection_only
visualize(hObject,handles);


% --- Executes on button press in add_curve.
function add_curve_Callback(hObject, eventdata, handles)
% hObject    handle to add_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function sorting_parameter_range_Callback(hObject, eventdata, handles)
% hObject    handle to sorting_parameter_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sorting_parameter_range as text
%        str2double(get(hObject,'String')) returns contents of sorting_parameter_range as a double


% --- Executes during object creation, after setting all properties.
function sorting_parameter_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sorting_parameter_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in curve_type.
function curve_type_Callback(hObject, eventdata, handles)
visualize(hObject,handles);

% --- Executes during object creation, after setting all properties.
function curve_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to curve_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clear_curves_axes.
function clear_curves_axes_Callback(hObject, eventdata, handles)
% hObject    handle to clear_curves_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
