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

% Last Modified by GUIDE v2.5 31-Oct-2018 12:27:11

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


function visualize(handles)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tracks = handles.TrackPlotter_handles.raw_data;
    locmaxpeaks = zeros(1,numel(tracks));
    %peakstrengthindicator = zeros(1,numel(tracks));
    maxlength = -Inf;
    for k=1:numel(tracks)
        track = tracks{k};
        FRET_ratio = squeeze(track(:,4));
        % set side parts to zero
            L=length(FRET_ratio);dL = fix(L/6);
            FRET_ratio(1:dL)=0;
            FRET_ratio(L-dL:L)=0;
        % set side parts to zero
        nucleus_size = squeeze(track(:,7));
        [pks,locs]= findpeaks(FRET_ratio,'MINPEAKDISTANCE',2);
        maxpks = max(pks);
        if ~isempty(maxpks)
            locmaxpeaks(k) = locs(pks==maxpks(1)); 
        else
            maxpks=max(FRET_ratio);
            locmaxpeaks(k) = find(FRET_ratio==maxpks(1));
        end
        l=numel(FRET_ratio);
        if l>maxlength
            maxlength=l;
        end
        %peakstrengthindicator(k) = max(pks)/median(pks);
    end
    %[~,indices]=sort(peakstrengthindicator);
    
    par_ind = get(handles.features_chooser,'Value');
    
     D = handles.TrackPlotter_handles.track_data;
%      par_ind = 1; % duration
%      par_ind = 2; % velocity 
%      par_ind = 3; % directionality
%      par_ind = 4; % nnghb          
%      par_ind = 5; % cell density
%      par_ind = 6; % FRET ratio
%      par_ind = 7; % FRET ratio variability     
%      par_ind = 8; % donor 
%      par_ind = 9; % acceptor
%      par_ind = 10; % nuclear size      
%      par_ind = 11; % Pearson corr
%      par_ind = 12; % start time
%      par_ind = 13; % autocorr
%      par_ind = 14; % beta FRET
        
     sortarr = D(:,par_ind+2);
     [~,indices]=sort(sortarr);
     indices=flip(indices); 
    
    Ldiv2 = floor(maxlength/2);
    L = 2*Ldiv2+1;
    M = zeros(numel(tracks),L);
       for k=1:numel(tracks)
            track = tracks{indices(k)};
            FRET_ratio = squeeze(track(:,4));
            l=length(FRET_ratio);
            shift = locmaxpeaks(indices(k));
            c = Ldiv2+1;
            range = c-shift:c+l-shift-1;
            try
                M(k,range)=FRET_ratio';
            catch
                disp([shift l L]);
            end
       end    
       s=M(M~=0);M(M==0)=min(s(:)); 
        
%     import bioma.data.*
%     M=map(M,0,1);
%     DMobj = DataMatrix(M);
%     hmo = HeatMap(DMobj,'Colormap',redgreencmap);
imagesc(handles.heatmap_axes,M);
xlabel(handles.heatmap_axes,'time [h]');
str = get(handles.features_chooser,'String'); 
ylabel(handles.heatmap_axes,str{par_ind});
tax =  handles.TrackPlotter_handles.dt*(1:L);
set(handles.heatmap_axes, 'xticklabel', {linspace(min(tax),max(tax),11)}, ... 
    'yticklabel', {flip(linspace(min(sortarr),max(sortarr),11))});
