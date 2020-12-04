function varargout = t_dependent_Nuclei_ratio_FRET_HTS_heatmapper(varargin)
% T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_HEATMAPPER MATLAB code for t_dependent_Nuclei_ratio_FRET_HTS_heatmapper.fig
%      T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_HEATMAPPER, by itself, creates a new T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_HEATMAPPER or raises the existing
%      singleton*.
%
%      H = T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_HEATMAPPER returns the handle to a new T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_HEATMAPPER or the handle to
%      the existing singleton*.
%
%      T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_HEATMAPPER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_HEATMAPPER.M with the given input arguments.
%
%      T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_HEATMAPPER('Property','Value',...) creates a new T_DEPENDENT_NUCLEI_RATIO_FRET_HTS_HEATMAPPER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before t_dependent_Nuclei_ratio_FRET_HTS_heatmapper_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to t_dependent_Nuclei_ratio_FRET_HTS_heatmapper_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help t_dependent_Nuclei_ratio_FRET_HTS_heatmapper

% Last Modified by GUIDE v2.5 04-Dec-2020 20:01:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @t_dependent_Nuclei_ratio_FRET_HTS_heatmapper_OpeningFcn, ...
                   'gui_OutputFcn',  @t_dependent_Nuclei_ratio_FRET_HTS_heatmapper_OutputFcn, ...
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


% --- Executes just before t_dependent_Nuclei_ratio_FRET_HTS_heatmapper is made visible.
function t_dependent_Nuclei_ratio_FRET_HTS_heatmapper_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to t_dependent_Nuclei_ratio_FRET_HTS_heatmapper (see VARARGIN)

% Choose default command line output for t_dependent_Nuclei_ratio_FRET_HTS_heatmapper
handles.output = hObject;

handles.TrackPlotter_handles = varargin{1};
[~,NAME,EXT]=fileparts(handles.TrackPlotter_handles.fullfilename);
set(handles.figure1, 'Name', ['HTS mapper : ' NAME EXT]);
set(handles.parameter,'String',handles.TrackPlotter_handles.features);

handles.letters = {'A','B','C','D','E','F','G','H'};
handles.numbers = {'1','2','3','4','5','6','7','8','9','10','11','12'};

set(handles.evaluation_time,'String','1');
set(handles.statistic,'String',{'p-value: KS','p-value: t-test','p-value: Wilcoxon','Cohen"s d','|median diff|'});
set(handles.letter_start,'String',handles.letters);
set(handles.letter_end,'String',handles.letters);
set(handles.number_start,'String',handles.numbers);
set(handles.number_end,'String',handles.numbers);
set(handles.remove_outliers,'Value',0);
set(handles.vizualization_mode,'String',{'platemap','chart'});

set(handles.parameter,'Value',6); % FRET ratio
set(handles.number_start,'Value',1);
set(handles.number_end,'Value',12);
set(handles.letter_start,'Value',1);
set(handles.letter_end,'Value',8);
set(handles.vizualization_mode,'Value',2);

handles.max_frame = floor(handles.TrackPlotter_handles.filter_table.Data(1,2)/handles.TrackPlotter_handles.dt);
handles.cur_frame = floor(str2double(get(handles.evaluation_time,'String'))/handles.TrackPlotter_handles.dt);

% define well tokens and setup masked data
    handles.raw_data = handles.TrackPlotter_handles.ST_raw_data(handles.TrackPlotter_handles.mask);
    raw_filenames = handles.TrackPlotter_handles.ST_raw_data_filenames(handles.TrackPlotter_handles.mask);
    %
    all_tokens_10_12 = [];
    all_tokens = [];    
    for k=1:8
        for m=1:12
            if m>=10
                all_tokens_10_12 = [all_tokens_10_12 {[handles.letters{k} '-' handles.numbers{m}]} ];
            end
            all_tokens = [all_tokens {[handles.letters{k} '-' handles.numbers{m}]} ];
        end
    end
    raw_data_tokens = num2cell(zeros(size(raw_filenames)));
    % first pass
    for k=1:numel(raw_data_tokens)
        for m=1:numel(all_tokens)
        if contains(raw_filenames{k},all_tokens{m})
            raw_data_tokens(k) = all_tokens(m);
        end
        end
    end
    % second pass
    for k=1:numel(raw_data_tokens)
        for m=1:numel(all_tokens_10_12)
        if contains(raw_filenames{k},all_tokens_10_12{m})
            raw_data_tokens(k) = all_tokens_10_12(m);
        end
        end
    end
            
    handles.raw_data_tokens = raw_data_tokens;
    
%%%%%%%%%%%%%%%%%% continued 30_11_2020    
    handles.diagram_panel = uipanel(handles.figure1,'visible','off');
    handles.diagram_panel.Position = [0.15 0.05 .7 .8];
    
    handles.platemap_panel = uipanel(handles.figure1,'visible','off');
    handles.platemap_panel.Position = [0.1 0.20 .8 .5];    
%%%%%%%%%%%%%%%%%% continued 30_11_2020    
        
    update_diagram_Callback(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes t_dependent_Nuclei_ratio_FRET_HTS_heatmapper wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = t_dependent_Nuclei_ratio_FRET_HTS_heatmapper_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% % --- Executes on button press in update_diagram.
% function update_diagram_Callback(hObject, eventdata, handles)
% 
%     p = uipanel(handles.figure1);
%     p.Position = [0.15 0.05 .7 .8];
%         
%     [selected_wells,D] = create_heatmap_data(handles);
%     
%     handles.h = heatmap(p,selected_wells,selected_wells,D);
%     handles.h.Colormap = jet;
%     
%     s = get(handles.parameter,'String');
%     param_name = s{get(handles.parameter,'Value')};
%     s = get(handles.statistic,'String');
%     statistic_name = s{get(handles.statistic,'Value')};
%     
%     handles.h.Title = ['t = ' get(handles.evaluation_time,'String'), ...    
%                        'h , ' param_name,' , ',statistic_name];
%     
%     bckg_color = get(handles.figure1,'Color');
%     handles.h.MissingDataColor = bckg_color;    
%     handles.h.GridVisible = 'off';    
% % Update handles structure
% handles.D = D;
% guidata(hObject, handles);

% --- Executes on button press in update_diagram.
function update_diagram_Callback(hObject, eventdata, handles)

    s = get(handles.parameter,'String');
    param_name = s{get(handles.parameter,'Value')};
    
    s = get(handles.statistic,'String');
    statistic_name = s{get(handles.statistic,'Value')};

    mode  = get(handles.vizualization_mode,'value');
    
    if 2==mode % chart    
        set(handles.diagram_panel,'visible','on');  
        set(handles.platemap_panel,'visible','off');
        set(handles.statistic,'String',{'p-value: KS','p-value: t-test','p-value: Wilcoxon','Cohen"s d','|median diff|'});
        set(handles.statistic,'visible','on'); 
        set(handles.generate_t_dependence,'Enable','off');
        set(handles.save_current_time_slice,'Enable','off');
        guidata(hObject, handles);

        [selected_wells,D] = create_heatmap_data(handles);

        handles.h = heatmap(handles.diagram_panel,selected_wells,selected_wells,D);
        handles.h.Colormap = jet;

        handles.h.Title = ['t = ' get(handles.evaluation_time,'String'), ...    
                           'h , ' param_name,' , ',statistic_name];

        bckg_color = get(handles.figure1,'Color');
        handles.h.MissingDataColor = bckg_color;    
        handles.h.GridVisible = 'off';    
        % Update handles structure
        handles.D = D; % look slike, not needed
        guidata(hObject, handles);
    else % platemap
        set(handles.diagram_panel,'visible','off');  
        set(handles.platemap_panel,'visible','on');
        set(handles.statistic,'String',{'mean','std','median','range','skewness','kurtosis'});
        set(handles.statistic,'visible','on');
        set(handles.generate_t_dependence,'Enable','on');
        set(handles.save_current_time_slice,'Enable','on');
        guidata(hObject, handles);

        [~,D] = create_platemap_data(handles);    
        handles.h = heatmap(handles.platemap_panel,handles.numbers,handles.letters,D);
        handles.h.Title = ['t = ' get(handles.evaluation_time,'String'), ...    
                           'h , ' param_name,' , ',statistic_name];    
        bckg_color = get(handles.figure1,'Color');
        handles.h.MissingDataColor = bckg_color;        
        handles.h.Colormap = jet;
    end
    
 guidata(hObject, handles);



% --- Executes on selection change in statistic.
function statistic_Callback(hObject, eventdata, handles)
% hObject    handle to statistic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns statistic contents as cell array
%        contents{get(hObject,'Value')} returns selected item from statistic


% --- Executes during object creation, after setting all properties.
function statistic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to statistic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in letter_start.
function letter_start_Callback(hObject, eventdata, handles)
end_letter = get(handles.letter_end,'Value');
my_letter = get(hObject,'Value');
if my_letter>end_letter
    set(hObject,'Value',end_letter);
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function letter_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to letter_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in letter_end.
function letter_end_Callback(hObject, eventdata, handles)
start_letter = get(handles.letter_start,'Value');
my_letter = get(hObject,'Value');
if my_letter<start_letter
    set(hObject,'Value',start_letter);
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function letter_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to letter_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in number_start.
function number_start_Callback(hObject, eventdata, handles)
end_number = get(handles.number_end,'Value');
my_number = get(hObject,'Value');
if my_number>end_number
    set(hObject,'Value',end_number);
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function number_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in number_end.
function number_end_Callback(hObject, eventdata, handles)
start_number = get(handles.number_start,'Value');
my_number = get(hObject,'Value');
if my_number<start_number
    set(hObject,'Value',start_number);
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function number_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in parameter.
function parameter_Callback(hObject, eventdata, handles)
% hObject    handle to parameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns parameter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from parameter


% --- Executes during object creation, after setting all properties.
function parameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function evaluation_time_Callback(hObject, eventdata, handles)
    suggested_frame = floor(str2double(get(handles.evaluation_time,'String'))/handles.TrackPlotter_handles.dt);
    if suggested_frame >= 1 && suggested_frame <= handles.max_frame
        handles.cur_frame = suggested_frame;
    end
    % Update handles structure
    guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function evaluation_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to evaluation_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function selected_wells = select_wells(handles)
selected_wells = cell(0);
for k=get(handles.letter_start,'Value'):get(handles.letter_end,'Value')
    for m=get(handles.number_start,'Value'):get(handles.number_end,'Value')
        selected_wells = [selected_wells [handles.letters{k} '-' handles.numbers{m}] ];
    end
end

function [selected_wells,D] = create_heatmap_data(handles) 
    %    
    feature_index = get(handles.parameter,'Value');
    
                switch feature_index                    
                    case 4 % #nghbrs                    
                        param_index = 9;
                    case 5 % cell density
                        param_index = 10;
                    case 6 % FRET ratio
                        param_index = 4;
                    case 8 % donor intensity                     
                        param_index = 5;
                    case 9 % acceptor intensity
                        param_index = 6;
                    case 10 % nucleus size
                        param_index = 7;
                    case 11 % Pearson                      
                        param_index = 8;
                    case 14 % FRET molar fraction
                        param_index = 11;
                    case 15 % cell size
                        param_index = 12;
                    case 16 % 
                        param_index = 13;
                    case 17 % 
                        param_index = 14;
                    case 18 % 
                        param_index = 15;
                    case 2 % speed                        
                        %Y = squeeze(handles.velocity_t{k});
                        param_index = 4;
                    case 19 % nuc_cell_are_ratio ??
                        % Y = squeeze(handles.nuc_cell_area_ratio_t{k});                                                
                        param_index = 4;                        
                end                                        
    %
    selected_wells = select_wells(handles);
    N = numel(selected_wells);
    %
    sample = cell(N,1); % statistical samples

    % run over all tracks at frame f 
    ndata = numel(handles.raw_data);   
    for k = 1:ndata        
            track_k = handles.raw_data{k};
            frame_index = find(track_k(:,1)==handles.cur_frame);
            if ~isempty(frame_index)
                value = track_k(frame_index,param_index);
                index = find(ismember(selected_wells,handles.raw_data_tokens{k})); 
                if ~isempty(index)
                    sample{index} = [sample{index} value];
                end
            end
    end
               
    D = nan(N);
    
    hw = waitbar(0,'calculating statistics..','WindowStyle','modal');
    for k=1:numel(sample)
        for m=1:numel(sample)
                        x1 = sample{k};
                        x2 = sample{m};
                        if get(handles.remove_outliers,'Value')
                            x1=rmoutliers(x1,'median');
                            x2=rmoutliers(x2,'median');                            
                        end
            if ~isempty(x1) && ~isempty(x2) && k<m
                switch(get(handles.statistic,'Value'))
                    case 1
                        [~,P] = kstest2(x1,x2);
                    case 2
                        [~,P] = ttest2(x1,x2,'vartype','unequal'); % not sure                        
                    case 3  
                        [P,~] = ranksum(x1,x2);
                    case 4 % Cohen's d
                        N1 = numel(x1);
                        N2 = numel(x2);
                        s = sqrt( 1/(N1+N2)*( (N1-1)*var(x1) + (N2-1)*var(x2) ) );
                        P = abs( mean(x1) - mean(x2) )/s;
                    case 5
                        P = abs(median(x1)-median(x2));
                    case 6 
                        P = nan;
                end
                D(k,m) = P; 
            end            
        end
        if ~isempty(hw), waitbar(k/numel(sample),hw); drawnow, end   
    end
    if ~isempty(hw), delete(hw), drawnow; end
 
% not sure about this, but having "sample" object, one can try more options

% https://stackoverflow.com/questions/7187945/adjust-p-values-for-multiple-comparisons-in-matlab

% load fisheriris
% [pVal tbl stats] = kruskalwallis(meas(:,1), species)   %# Kruskal-Wallis or ANOVA
% title('Sepal Length'), xlabel('Groups'), ylabel('Value')
% 
% [c,m] = multcompare(stats, 'ctype','bonferroni', 'display','on');


% --- Executes on button press in remove_outliers.
function remove_outliers_Callback(hObject, eventdata, handles)
% hObject    handle to remove_outliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of remove_outliers


% continued designing...

function [selected_wells,D] = create_platemap_data(handles) 

    tracks = handles.TrackPlotter_handles.raw_data;
    wells = handles.raw_data_tokens;    
    %
    feature_index = get(handles.parameter,'Value');
    
                switch feature_index                    
                    case 4 % #nghbrs                    
                        param_index = 9;
                    case 5 % cell density
                        param_index = 10;
                    case 6 % FRET ratio
                        param_index = 4;
                    case 8 % donor intensity                     
                        param_index = 5;
                    case 9 % acceptor intensity
                        param_index = 6;
                    case 10 % nucleus size
                        param_index = 7;
                    case 11 % Pearson                      
                        param_index = 8;
                    case 14 % FRET molar fraction
                        param_index = 11;
                    case 15 % cell size
                        param_index = 12;
                    case 16 % 
                        param_index = 13;
                    case 17 % 
                        param_index = 14;
                    case 18 % 
                        param_index = 15;
                    case 2 % speed                        
                        %Y = squeeze(handles.velocity_t{k});
                        param_index = 4;
                    case {1,3,19,7,13} % nuc_cell_are_ratio ??
                        % Y = squeeze(handles.nuc_cell_area_ratio_t{k});                                                
                        param_index = 4;                        
                end                                        

    selected_wells = select_wells(handles);

    C = cell(8,12);

    % run over all tracks at frame f 
    ndata = numel(tracks);   
    for dat = 1:ndata        
            track = tracks{dat};
            frame_index = find(track(:,1)==handles.cur_frame);
            index = find(ismember(selected_wells,wells{dat}));
            if ~isempty(frame_index) && ~isempty(index)
                value = track(frame_index,param_index);
                [k,m] = get_letter_number_indices(handles,wells{dat});
                C{k,m} = [C{k,m}; value];
            end
    end
    %    
    D = NaN(8,12);
    for k=1:8
        for m=1:12
            x = C{k,m};
            if ~isempty(x) && numel(x)>1
                        if get(handles.remove_outliers,'Value')
                            x=rmoutliers(x,'median');
                        end
                    %  set(handles.statistic,'String',{'mean','std','median','range','skewness','kurtosis'});
                    switch(get(handles.statistic,'Value'))
                        case 1
                            D(k,m) = mean(x);
                        case 2
                            D(k,m) = std(x);
                        case 3
                            D(k,m) = median(x);
                        case 4
                            D(k,m) = max(x) - min(x);
                        case 5
                            D(k,m) = skewness(x);
                        case 6
                            D(k,m) = kurtosis(x);                            
                    end
            end                    
        end
    end

function [letter_ind,number_ind] = get_letter_number_indices(handles,well_token)
    s = strsplit(well_token,'-');
    letter_ind = find(ismember(handles.letters,s{1}));
    number_ind = find(ismember(handles.numbers,s{2}));

% --- Executes on selection change in vizualization_mode.
function vizualization_mode_Callback(hObject, eventdata, handles)
update_diagram_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function vizualization_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vizualization_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in generate_t_dependence.
function generate_t_dependence_Callback(hObject, eventdata, handles)
    
t = (0:(handles.max_frame-1))*handles.TrackPlotter_handles.dt;

f_sample = cell(numel(t),1);

    tracks = handles.TrackPlotter_handles.raw_data;
    wells = handles.raw_data_tokens;    
    %
    feature_index = get(handles.parameter,'Value');
    
                switch feature_index                    
                    case 4 % #nghbrs                    
                        param_index = 9;
                    case 5 % cell density
                        param_index = 10;
                    case 6 % FRET ratio
                        param_index = 4;
                    case 8 % donor intensity                     
                        param_index = 5;
                    case 9 % acceptor intensity
                        param_index = 6;
                    case 10 % nucleus size
                        param_index = 7;
                    case 11 % Pearson                      
                        param_index = 8;
                    case 14 % FRET molar fraction
                        param_index = 11;
                    case 15 % cell size
                        param_index = 12;
                    case 16 % 
                        param_index = 13;
                    case 17 % 
                        param_index = 14;
                    case 18 % 
                        param_index = 15;
                    case 2 % speed                        
                        %Y = squeeze(handles.velocity_t{k});
                        param_index = 4;
                    case {1,3,19,7,13} % nuc_cell_are_ratio ??
                        % Y = squeeze(handles.nuc_cell_area_ratio_t{k});                                                
                        param_index = 4;                        
                end                                        

    selected_wells = select_wells(handles);

    hw = waitbar(0,'gathering statistics..','WindowStyle','modal');
    ndata = numel(tracks);   
    for dat = 1:ndata
        %[ dat ndata ]
            track = tracks{dat};
            index = find(ismember(selected_wells,wells{dat}));
            if ~isempty(index)
                for f = 1:numel(f_sample)
                    minmax_f = minmax(track(:,1)');
                    if f>=minmax_f(1) && f<=minmax_f(2)
                        f_ind = find(f==track(:,1));
                        value = track(f_ind,param_index);
                        f_sample{f} = [f_sample{f}; value];
                    end                
                end
            end
        if ~isempty(hw), waitbar(dat/ndata,hw); drawnow, end
    end
    if ~isempty(hw), delete(hw), drawnow; end            
    %
    mode = get(handles.statistic,'Value');
    %
    D = NaN(numel(f_sample),1);
    D_std = D;
    for f=1:numel(f_sample)
            x = f_sample{f};
            if ~isempty(x) && numel(x)>1
                        if get(handles.remove_outliers,'Value')
                            x=rmoutliers(x,'median');
                        end
                    %  set(handles.statistic,'String',{'mean','std','median','range','skewness','kurtosis'});
                    if intersect(mode,[1,2,3])
                        D_std(f) = std(x);
                    end
                    %
                    switch mode
                        case 1
                            D(f) = mean(x);
                        case 2
                            D(f) = D_std(f);
                        case 3
                            D(f) = median(x);
                        case 4
                            D(f) = max(x) - min(x);
                        case 5
                            D(f) = skewness(x);
                        case 6
                            D(f) = kurtosis(x);                            
                    end
            end                    
    end
    
    token1 = selected_wells{1};
    token2 = selected_wells{end};
    h = figure;
    set(h, 'Name', ['HTS mapper : ' token1 ' ... ' token2]);
    
    if intersect(mode,[1,3])
        plot(gca,t,D,'k.-','linewidth',3);
        hold(gca,'on');
        shadedErrorBar(t,D,D_std,'lineProps',{'k.-','linewidth',2,'markersize',24});
        hold(gca,'off');        
    else
        plot(gca,t,D,'k.-','linewidth',3);
    end
    grid(gca,'on');    
    xlabel(gca,'time [h]');
    index = get(handles.statistic,'Value');
    s = get(handles.statistic,'String');
    statistic_name = s{index};
    s = get(handles.parameter,'String');
    parameter_name = s{feature_index};
    ylabel([statistic_name ' ( ' parameter_name ' ) ']);
            
    uiresume(handles.figure1);
    

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function save_current_time_slice_Callback(hObject, eventdata, handles)

    tracks = handles.TrackPlotter_handles.raw_data;
    wells = handles.raw_data_tokens;    
        
    records = [];
    
    nparams = size(tracks{1},2);
    
    selected_wells = select_wells(handles);

   % run over all tracks at frame f 
    hw = waitbar(0,'gathering statistics..','WindowStyle','modal');
    ndata = numel(tracks);   
    for dat = 1:ndata        
            track = tracks{dat};
            frame_index = find(track(:,1)==handles.cur_frame);
            index = find(ismember(selected_wells,wells{dat}));
            if ~isempty(frame_index) && ~isempty(index)
                record = track(frame_index,:);
                records = [records; [wells{dat} num2cell(record)] ];
            end
        if ~isempty(hw), waitbar(dat/ndata,hw); drawnow, end
    end
    if ~isempty(hw), delete(hw), drawnow; end            

[filename, pathname] = uiputfile( ...
       {'*.xls';'*.mat';'*.*'}, ...
        'Save as');
    
    if contains(filename,'xls')
        if 11==nparams
            caption = {'well','frame','xc','yc','FRET ratio','ID','IA','nuc.size','Pearson','#nghb','cell density','beta FRET'};
        elseif 15 == nparams
            caption = {'well','frame','xc','yc','FRET ratio','ID','IA','nuc.size','Pearson','#nghb','cell density','beta FRET','cell area','ref I cell','ref I nuc','ref I nuc/cell'};
        end
        xlswrite([pathname filesep filename],[caption; records]);
    else
        save([pathname filesep filename],'records');
    end
