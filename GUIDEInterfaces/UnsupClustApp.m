function varargout = UnsupClustApp(varargin)
% UNSUPCLUSTAPP MATLAB code for UnsupClustApp.fig
%      UNSUPCLUSTAPP, by itself, creates a new UNSUPCLUSTAPP or raises the existing
%      singleton*.
%
%      H = UNSUPCLUSTAPP returns the handle to a new UNSUPCLUSTAPP or the handle to
%      the existing singleton*.
%
%      UNSUPCLUSTAPP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNSUPCLUSTAPP.M with the given input arguments.
%
%      UNSUPCLUSTAPP('Property','Value',...) creates a new UNSUPCLUSTAPP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UnsupClustApp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UnsupClustApp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UnsupClustApp

% Last Modified by GUIDE v2.5 09-Aug-2021 18:33:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UnsupClustApp_OpeningFcn, ...
                   'gui_OutputFcn',  @UnsupClustApp_OutputFcn, ...
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


% --- Executes just before UnsupClustApp is made visible.
function UnsupClustApp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to UnsupClustApp (see VARARGIN)

% Choose default command line output for UnsupClustApp
handles.output = hObject;

handles.maiden_name = get(handles.figure1,'Name');

set(handles.number_of_clusters,'String',{'1','2','3','4','5','6','7','auto'});
set(handles.embedding_method,'String',{'PCA cov','PCA corr','t-SNE','none'});
set(handles.clustering_method,'String',{'k-means','k-medoids','hierarchical','GMM'});
set(handles.clusters_view,'String',{'C1,C2','C1,C2,C3'});
set(handles.manova1_coords,'Value',1);
%
set(handles.parameters_table,'CellEditCallback',@parameters_table_callback);
set(handles.conditions_table,'CellEditCallback',@conditions_table_callback);

handles.acting_conditions = [];
handles.acting_data = [];
handles.conditions_choice = [];
handles.parameters_choice = [];
handles.conditions_names = [];

        handles.colors = zeros(7,3);
        handles.colors(1,:) = [0 0.4470 0.7410];
        handles.colors(2,:) = [0.8500 0.3250 0.0980];
        handles.colors(3,:) = [0.9290 0.6940 0.1250];
        handles.colors(4,:) = [0.4940 0.1840 0.5560];
        handles.colors(5,:) = [0.4660 0.6740 0.1880];
        handles.colors(6,:) = [0.3010 0.7450 0.9330];
        handles.colors(7,:) = [0.6350 0.0780 0.1840];
        handles.markers = {'o','s','^','d','p','*','+'};
        handles.styles = {'-',':','--','-.','-',':','--'};
        %
        handles.IDX = []; 
        handles.transformed_data = [];       
        %
        handles.n_clusters = 3;        
        %
        set(handles.clusters_conditions_distribution_panel,'Visible','off');
        set(handles.clusters_conditions_distribution_table,'Visible','off');               
        %
        set(handles.clusters_names,'Data',{'type 1','type 2','type 3'}');
        %
        handles.data = [];

if 4==length(varargin) % try to load data
    data = varargin{1};
    parameters = varargin{2};
    conditions = varargin{3};
    datafname = varargin{4};
        if isempty(data), return, end
    %
    handles.data = data;
    handles.parameters = parameters;
    %remove underscores if present
    for k=1:numel(handles.parameters)
        handles.parameters{k} = strrep(handles.parameters{k},'_','-');
    end
    %
    handles.conditions = conditions;
    clear('data');
    clear('parameters');
    clear('conditions');
    %
    cla(handles.clusters_axis,'reset');
    cla(handles.histograms,'reset');
    %
    % set up for neew data
    set(handles.histo_param_1,'String',handles.parameters);
    %
    set(handles.clusters_view,'Value',1);
    rotate3d(handles.clusters_axis,'off');
    rotate3d(handles.histograms,'off');
    %
    set(handles.figure1, 'Name', [handles.maiden_name ' : ' datafname]);
    %
    L = length(handles.parameters);
    % calculate data minmax
    handles.data_minmax = [min(handles.data,[],1)' max(handles.data,[],1)'];
    handles.data_minmax_vanilla = handles.data_minmax; % to roll back to default if needed
    %        
    data = [num2cell(true(L,1)) handles.parameters' ... 
        num2cell(handles.data_minmax(:,1)) num2cell(handles.data_minmax(:,2))];
    set(handles.parameters_table, 'Data', data, ... 
        'ColumnEditable',[true,false]);
    set(handles.parameters_table,'ColumnEditable', [true, false, true, true]);
    %    
    handles.conditions_names = unique(handles.conditions);
    L = length(handles.conditions_names);
    data = [num2cell(true(L,1)) handles.conditions_names];
    set(handles.conditions_table, 'Data', data, ... 
        'ColumnEditable',[true,false]);    
    %
    set(handles.number_of_clusters,'Value',3);
    data = {'type 1','type 2','type 3'}';
    set(handles.clusters_names,'Data', data);
    handles.n_clusters = 3;
    
    handles.conditions_choice = ones(numel(handles.conditions_names),1);
    handles.parameters_choice = ones(numel(handles.parameters),1);
    %
    set(handles.clusters_conditions_distribution_panel,'Visible','off');
    set(handles.clusters_conditions_distribution_table,'Visible','off');    
    
end
        
% Update handles structure

guidata(hObject,handles);   
uiresume(handles.figure1); 

if ~isempty(handles.data)
    calculate_Callback(hObject, eventdata, handles);
end

% UIWAIT makes UnsupClustApp wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = UnsupClustApp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in embedding_method.
function embedding_method_Callback(hObject, eventdata, handles)
% hObject    handle to embedding_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns embedding_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from embedding_method


% --- Executes during object creation, after setting all properties.
function embedding_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to embedding_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in clustering_method.
function clustering_method_Callback(hObject, eventdata, handles)
% hObject    handle to clustering_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns clustering_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from clustering_method


% --- Executes during object creation, after setting all properties.
function clustering_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clustering_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in number_of_clusters.
function number_of_clusters_Callback(hObject, eventdata, handles)
    str = get(hObject,'String');
    token = str{get(hObject,'Value')};
    token_num = str2num(token);
    if isnumeric(token_num)
        handles.n_clusters = token_num;
    else
        %
        % maybe..
        %
    end
    namestr = cell(handles.n_clusters,1);
    for k=1:handles.n_clusters
        namestr{k} = ['type ' num2str(k)];
    end
    %
    set(handles.clusters_names,'Data',namestr);
    set(handles.clusters_names,'ColumnEditable', true);
    guidata(hObject, handles);
    uiresume(handles.figure1);     

% --- Executes during object creation, after setting all properties.
function number_of_clusters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_of_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in clusters_view.
function clusters_view_Callback(hObject, eventdata, handles)
if 1==get(hObject,'Value')
    rotate3d(handles.clusters_axis,'off');
else
    rotate3d(handles.clusters_axis,'on');   
end
visualize_feature_space(hObject,handles);

% --- Executes during object creation, after setting all properties.
function clusters_view_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clusters_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in histo_param_1.
function histo_param_1_Callback(hObject, eventdata, handles)
visualize_histograms(hObject,handles);

% --- Executes during object creation, after setting all properties.
function histo_param_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to histo_param_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in histo_param_2.
function histo_param_2_Callback(hObject, eventdata, handles)
% hObject    handle to histo_param_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns histo_param_2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from histo_param_2


% --- Executes during object creation, after setting all properties.
function histo_param_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to histo_param_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in histo_logscale.
function histo_logscale_Callback(hObject, eventdata, handles)
visualize_histograms(hObject,handles);

% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function load_data_Callback(hObject, eventdata, handles)
    %
    [data,parameters,conditions,datafname] = load_datafile();
    %
    if isempty(data), return, end
    %
    handles.data = data;
    handles.parameters = parameters;
    %remove underscores if present
    for k=1:numel(handles.parameters)
        handles.parameters{k} = strrep(handles.parameters{k},'_','-');
    end
    %
    handles.conditions = conditions;
    clear('data');
    clear('parameters');
    clear('conditions');
    %
    cla(handles.clusters_axis,'reset');
    cla(handles.histograms,'reset');
    %
    guidata(hObject,handles);   
    uiresume(handles.figure1);      
    %    
    setup_for_new_data(hObject,handles,datafname);
    
    % ? doesn't work..
    %calculate_Callback(hObject,eventdata,handles);
        
function setup_for_new_data(hObject,handles,datafname) 
    if isempty(handles.data), return, end
    %
    set(handles.histo_param_1,'String',handles.parameters);
    %
    set(handles.histo_param_1,'Value',1);
    %
    set(handles.clusters_view,'Value',1);
    rotate3d(handles.clusters_axis,'off');
    rotate3d(handles.histograms,'off');
    %
    set(handles.figure1, 'Name', [handles.maiden_name ' : ' datafname]);    
    %
    L = length(handles.parameters);
    % calculate data minmax
    handles.data_minmax = [min(handles.data,[],1)' max(handles.data,[],1)'];
    handles.data_minmax_vanilla = handles.data_minmax; % to roll back to default if needed
    %        
    data = [num2cell(true(L,1)) handles.parameters' ... 
        num2cell(handles.data_minmax(:,1)) num2cell(handles.data_minmax(:,2))];
    set(handles.parameters_table, 'Data', data, ... 
        'ColumnEditable',[true,false]);
    set(handles.parameters_table,'ColumnEditable', [true, false, true, true]);
    %    
    handles.conditions_names = unique(handles.conditions);
    L = length(handles.conditions_names);
    data = [num2cell(true(L,1)) handles.conditions_names];
    set(handles.conditions_table, 'Data', data, ... 
        'ColumnEditable',[true,false]);    
    %
    set(handles.number_of_clusters,'Value',3);
    data = {'type 1','type 2','type 3'}';
    set(handles.clusters_names,'Data', data);
    handles.n_clusters = 3;
    
    handles.conditions_choice = ones(numel(handles.conditions_names),1);
    handles.parameters_choice = ones(numel(handles.parameters),1);
    %
    set(handles.clusters_conditions_distribution_panel,'Visible','off');
    set(handles.clusters_conditions_distribution_table,'Visible','off');         
    %
    guidata(hObject,handles);
    uiresume(handles.figure1);  
    
function [data,parameters,conditions,fname] = load_datafile() 
    data = [];
    parameters = [];
    conditions = [];    
             [fname, fpath] = uigetfile({'*.mat;*.csv;*.xls;*.xlsx'},'Select data file',pwd);
             if fpath == 0, return, end
             try             
                 if contains(fname,'.mat')
                         load(fullfile(fpath,fname),'data','parameters','conditions');
                 else
                    [data,txt,raw] = xlsread(fullfile(fpath,fname));
                    parameters = txt(1,:);
                    parameters = parameters(2:length(parameters));
                    conditions = txt(2:size(txt,1),1);                 
                 end
             catch
                 disp('error while loading data file');
             end
             
%---------------------------
function parameters_table_callback(hObject,callbackdata)
try
    handles = guidata(hObject);
    handles.parameters_choice = cell2mat(handles.parameters_table.Data(:,1));
    % numbers consistency
    minval_nans = isnan(cell2mat(handles.parameters_table.Data(:,3)));
    maxval_nans = isnan(cell2mat(handles.parameters_table.Data(:,4)));    
    handles.parameters_table.Data(find(1==minval_nans),3)=num2cell(handles.data_minmax_vanilla(find(1==minval_nans),1));
    handles.parameters_table.Data(find(1==maxval_nans),4)=num2cell(handles.data_minmax_vanilla(find(1==maxval_nans),2));
    % minmax consistency
    minval_bads = cell2mat(handles.parameters_table.Data(:,3))>=cell2mat(handles.parameters_table.Data(:,4)); % min more than max
    maxval_bads = cell2mat(handles.parameters_table.Data(:,4))<=cell2mat(handles.parameters_table.Data(:,3)); % max less than min    
    handles.parameters_table.Data(find(1==minval_bads),3)=num2cell(handles.data_minmax_vanilla(find(1==minval_bads),1));
    handles.parameters_table.Data(find(1==maxval_bads),4)=num2cell(handles.data_minmax_vanilla(find(1==maxval_bads),2));
    % update
    handles.data_minmax(:,1) = cell2mat(handles.parameters_table.Data(:,3));
    handles.data_minmax(:,2) = cell2mat(handles.parameters_table.Data(:,4));
    %
    acting_parameters = handles.parameters(1==handles.parameters_choice);
    set(handles.histo_param_1,'Value',1);
    set(handles.histo_param_1,'String',acting_parameters);
    %
    % yepp!     
    guidata(hObject, handles);
    uiresume(handles.figure1);    
catch
end

%---------------------------
function conditions_table_callback(hObject,callbackdata)
try
    handles = guidata(hObject);
    handles.conditions_choice = cell2mat(handles.conditions_table.Data(:,1));
    guidata(hObject, handles);
    uiresume(handles.figure1);    
catch
end
        
% --- Executes on button press in set_vanilla_minmax.
function set_vanilla_minmax_Callback(hObject, eventdata, handles)
    handles.data_minmax = handles.data_minmax_vanilla;
    
    choice = handles.parameters_table.Data(:,1);
    
    data = [choice handles.parameters' ... 
        num2cell(handles.data_minmax(:,1)) num2cell(handles.data_minmax(:,2))];
    set(handles.parameters_table, 'Data', data, ... 
        'ColumnEditable',[true,false]);
    set(handles.parameters_table,'ColumnEditable', [true, false, true, true]);    
    guidata(hObject, handles);
    uiresume(handles.figure1); 

% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
    %
    if isempty(handles.data), return, end 
    %        
    acting_conditions_indices = ismember(handles.conditions,handles.conditions_names(1==handles.conditions_choice));
    %
    acting_data = handles.data(:,1==handles.parameters_choice);

    data_minmax = handles.data_minmax(1==handles.parameters_choice,:);

    % min max compliance
    min_compliant_indices = acting_data>=data_minmax(:,1)';
    max_compliant_indices = acting_data<=data_minmax(:,2)';

    L = size(min_compliant_indices,2);
    min_compliant_indices = L==sum(min_compliant_indices,2);
    max_compliant_indices = L==sum(max_compliant_indices,2);

    acting_indices = acting_conditions_indices & min_compliant_indices & max_compliant_indices;
    handles.acting_conditions = handles.conditions(acting_indices);    
    handles.acting_data = handles.data(acting_indices,1==handles.parameters_choice);
    %
    IDX = [];
    data = handles.acting_data;
    if isempty(data), return, end

    embedding_method = get(handles.embedding_method,'Value');
    %
    switch embedding_method 
        case 1 % PCA
            % no truncation
            truncation = size(data,2);
            V = cov(data);
            COEFF = pcacov(V);
            U = data*COEFF;
            data = U(:,1:truncation);
        case 2
            % no truncation
            truncation = size(data,2);            
            V = cov(data);
            SD = sqrt(diag(V));
            R = V./(SD*SD');
            COEFF = pcacov(R);                                    
            U = data*COEFF;
            data = U(:,1:truncation);
        case 3
            % no truncation
            truncation = size(data,2);
            V = cov(data);
            COEFF = pcacov(V);
            U = data*COEFF;
            data = U(:,1:truncation);
            % only after PCA...
            data = tsne(data);
        case 4 % none            
    end
        
    n_clusters = handles.n_clusters;

    switch get(handles.clustering_method,'Value')
        case 1
            % KMEANS
            IDX = kmeans(data,n_clusters);
            % IDX = kmeans(data,n_clusters,'OnlinePhase','on','Distance','correlation');
            % IDX = kmeans(data,n_clusters,'Replicates',10);
            % IDX = kmeans(data,n_clusters,'Start','cluster');
            % IDX = kmeans(data,n_clusters,'OnlinePhase','on','Distance','correlation','Start','cluster','Replicates',10);        
        case 2
            % KMEDOIDS
            IDX = kmedoids(data,n_clusters);
            %IDX = kmedoids(data,n_clusters,'Algorithm','clara');        
        case 3
            % HIERARCHICAL CLUSTERING
            Z = linkage(data,'ward');
            IDX = cluster(Z,'Maxclust',n_clusters);                       
        case 4
            options = statset('MaxIter',1000); 
            gmfit = fitgmdist(data,n_clusters,'Options',options);
            IDX = cluster(gmfit,data);            
    end
    %
    % reassign categories
    numdata = zeros(1,n_clusters);
    for k=1:n_clusters
        numdata(k)=sum(IDX==k);
    end
    sorted_numdata = sort(numdata,'descend');

    lut = zeros(1,n_clusters);
    for k=1:n_clusters
        for m=1:n_clusters
            if numdata(k)==sorted_numdata(m)
                break;
            end
        end
        lut(k)=m;
    end
    %
    for k=1:numel(IDX)
         IDX(k) = lut(IDX(k));
    end
    %
    handles.transformed_data = data;
    handles.IDX = IDX;
    
    guidata(hObject, handles);    
    %
    visualize_feature_space(hObject,handles);
    %
    visualize_histograms(hObject,handles);
    %
    visualize_cluster_condition_distribution(hObject,handles);
    %
    guidata(hObject, handles);
    uiresume(handles.figure1);     

% ------------------------------------------------------------------------    
function visualize_feature_space(hObject,handles)
    %
    data = handles.transformed_data;
    if isempty(data), return, end

    AX = handles.clusters_axis;
    
    cla(AX,'reset');

    clusters_names = get(handles.clusters_names,'Data');
    LEGEND = [];
        
    switch get(handles.manova1_coords,'Value')
        case 1
            % CANONICAL COORDINATES
            [~,~,stats] = manova1(data,handles.IDX);
            %
            for k=1:handles.n_clusters
                C_1 = stats.canon(handles.IDX==k,1);
                C_2 = stats.canon(handles.IDX==k,2);
                if 0 ~= numel(C_1)
                    if 1== get(handles.clusters_view,'Value')
                            plot(AX,C_1,C_2,'color',handles.colors(k,:),'marker',handles.markers{k},'linestyle','none','markersize',4,'linewidth',2);            
                    elseif 2== get(handles.clusters_view,'Value') && 3~=get(handles.embedding_method,'Value')
                            C_3 = stats.canon(handles.IDX==k,3);
                            plot3(AX,C_1,C_2,C_3,'color',handles.colors(k,:),'marker',handles.markers{k},'linestyle','none','markersize',4,'linewidth',2);
                    end
                hold(AX,'on');
                LEGEND = [LEGEND; {[num2str(k)  ' : '  clusters_names{k} ' : ' num2str(numel(C_1))]}];
                end
            end

        case 0
            for k=1:handles.n_clusters
                C_1 = data(handles.IDX==k,1);
                C_2 = data(handles.IDX==k,2);
                if 0 ~= numel(C_1)
                    if 1== get(handles.clusters_view,'Value')
                            plot(AX,C_1,C_2,'color',handles.colors(k,:),'marker',handles.markers{k},'linestyle','none','markersize',4,'linewidth',2);            
                    elseif 2== get(handles.clusters_view,'Value') && 3~=get(handles.embedding_method,'Value')
                            C_3 = data(handles.IDX==k,3);
                            plot3(AX,C_1,C_2,C_3,'color',handles.colors(k,:),'marker',handles.markers{k},'linestyle','none','markersize',4,'linewidth',2);
                    end
                hold(AX,'on');
                LEGEND = [LEGEND; {[num2str(k)  ' : '  clusters_names{k} ' : ' num2str(numel(C_1))]}];
                end
            end    
    end

    hold(AX,'off')
%         xticks(AX,[]);yticks(AX,[]);
%         set(AX,'XColor','none');
%         set(AX,'YColor','none');    
    legend(AX,LEGEND);
    xlabel(AX,'C_1','fontsize',8);
    ylabel(AX,'C_2','fontsize',8);
    zlabel(AX,'C_3','fontsize',8);    
    grid(AX,'on');
    %
    guidata(hObject, handles);    
    uiresume(handles.figure1);

% ------------------------------------------------------------------------    
function visualize_histograms(hObject,handles)
    if isempty(handles.acting_conditions) || isempty(handles.acting_data) || isempty(handles.IDX), return, end
    %
    histo_param_1_index = get(handles.histo_param_1,'Value');
    
    acting_parameters = handles.parameters(1==handles.parameters_choice);
    
    AX = handles.histograms;
    
            for k=1:handles.n_clusters
                X = handles.acting_data(k==handles.IDX,histo_param_1_index);
                histogram(AX,X,'facecolor',handles.colors(k,:),'facealpha',.5,'edgecolor','none','normalization','probability');                
                %[k length(sample)]
                hold(AX,'on');
            end
            hold(AX,'off');
            ylabel(AX,'probability');
            xlabel(AX,acting_parameters{histo_param_1_index});
            grid(AX,'on');
            legend(AX,get(handles.clusters_names,'Data'));
                if get(handles.histo_logscale,'Value')
                    set(AX,'YScale','log');
                end                 

% ------------------------------------------------------------------------    
function visualize_cluster_condition_distribution(hObject, handles)
if isempty(handles.acting_conditions) || isempty(handles.acting_data) || isempty(handles.IDX), return, end
%
acting_conditions_names = handles.conditions_names(1==handles.conditions_choice);

n_cond = numel(acting_conditions_names);

[n_row,n_col] = MCED(n_cond);

cat_data = zeros(n_cond,handles.n_clusters);
for cond=1:n_cond
    cur_cond_name = acting_conditions_names{cond};
    cond_mask = strcmp(handles.acting_conditions,cur_cond_name);
    for clus = 1:handles.n_clusters
        clus_mask = clus==handles.IDX;
        cat_data(cond,clus) = sum(clus_mask & cond_mask);
    end
end

set(handles.clusters_conditions_distribution_table,'Data',cat_data');
set(handles.clusters_conditions_distribution_table,'RowName',get(handles.clusters_names,'Data'));
set(handles.clusters_conditions_distribution_table,'ColumnName',acting_conditions_names);
set(handles.clusters_conditions_distribution_table,'Visible','on');

set(handles.clusters_conditions_distribution_panel,'Visible','on');
for k=1:n_cond
    ax = subplot(n_row,n_col,k,'Parent', handles.clusters_conditions_distribution_panel);
    cla(ax,'reset');
    %
    X = cat_data(k,:);        
    b = bar(diag(X),'stacked');
    % looks correct.. (?)
    for m=1:length(X)
        set(b(m),'FaceColor',handles.colors(m,:));
    end
    %set(ax,'YScale','log');
    title(ax,acting_conditions_names{k});
    grid(ax,'on');   
end
 
% --- Executes on button press in select_all_parameters.
function select_all_parameters_Callback(hObject, eventdata, handles)
    handles.data_minmax = handles.data_minmax_vanilla;

    L = numel(handles.parameters);    
    
    data = [num2cell(true(L,1)) handles.parameters' ... 
        num2cell(handles.data_minmax(:,1)) num2cell(handles.data_minmax(:,2))];
    set(handles.parameters_table, 'Data', data, ... 
        'ColumnEditable',[true,false]);
    set(handles.parameters_table,'ColumnEditable', [true, false, true, true]);   
    set(handles.histo_param_1,'String',handles.parameters);
    %
    handles.parameters_choice = cell2mat(handles.parameters_table.Data(:,1));
    % numbers consistency
    minval_nans = isnan(cell2mat(handles.parameters_table.Data(:,3)));
    maxval_nans = isnan(cell2mat(handles.parameters_table.Data(:,4)));    
    handles.parameters_table.Data(find(1==minval_nans),3)=num2cell(handles.data_minmax_vanilla(find(1==minval_nans),1));
    handles.parameters_table.Data(find(1==maxval_nans),4)=num2cell(handles.data_minmax_vanilla(find(1==maxval_nans),2));
    % minmax consistency
    minval_bads = cell2mat(handles.parameters_table.Data(:,3))>=cell2mat(handles.parameters_table.Data(:,4)); % min more than max
    maxval_bads = cell2mat(handles.parameters_table.Data(:,4))<=cell2mat(handles.parameters_table.Data(:,3)); % max less than min    
    handles.parameters_table.Data(find(1==minval_bads),3)=num2cell(handles.data_minmax_vanilla(find(1==minval_bads),1));
    handles.parameters_table.Data(find(1==maxval_bads),4)=num2cell(handles.data_minmax_vanilla(find(1==maxval_bads),2));
    % update
    handles.data_minmax(:,1) = cell2mat(handles.parameters_table.Data(:,3));
    handles.data_minmax(:,2) = cell2mat(handles.parameters_table.Data(:,4));
    %
    acting_parameters = handles.parameters(1==handles.parameters_choice);
    set(handles.histo_param_1,'Value',1);
    set(handles.histo_param_1,'String',acting_parameters);  
    %
    guidata(hObject, handles);
    uiresume(handles.figure1); 

% --- Executes on button press in unselect_all_parameters.
function unselect_all_parameters_Callback(hObject, eventdata, handles)
    handles.data_minmax = handles.data_minmax_vanilla;
        
    L = numel(handles.parameters);
    
    data = [num2cell(false(L,1)) handles.parameters' ... 
        num2cell(handles.data_minmax(:,1)) num2cell(handles.data_minmax(:,2))];
    set(handles.parameters_table, 'Data', data, ... 
        'ColumnEditable',[true,false]);
    set(handles.parameters_table,'ColumnEditable', [true, false, true, true]);    
    set(handles.histo_param_1,'String',handles.parameters);
    %
    handles.parameters_choice = cell2mat(handles.parameters_table.Data(:,1));
    % numbers consistency
    minval_nans = isnan(cell2mat(handles.parameters_table.Data(:,3)));
    maxval_nans = isnan(cell2mat(handles.parameters_table.Data(:,4)));    
    handles.parameters_table.Data(find(1==minval_nans),3)=num2cell(handles.data_minmax_vanilla(find(1==minval_nans),1));
    handles.parameters_table.Data(find(1==maxval_nans),4)=num2cell(handles.data_minmax_vanilla(find(1==maxval_nans),2));
    % minmax consistency
    minval_bads = cell2mat(handles.parameters_table.Data(:,3))>=cell2mat(handles.parameters_table.Data(:,4)); % min more than max
    maxval_bads = cell2mat(handles.parameters_table.Data(:,4))<=cell2mat(handles.parameters_table.Data(:,3)); % max less than min    
    handles.parameters_table.Data(find(1==minval_bads),3)=num2cell(handles.data_minmax_vanilla(find(1==minval_bads),1));
    handles.parameters_table.Data(find(1==maxval_bads),4)=num2cell(handles.data_minmax_vanilla(find(1==maxval_bads),2));
    % update
    handles.data_minmax(:,1) = cell2mat(handles.parameters_table.Data(:,3));
    handles.data_minmax(:,2) = cell2mat(handles.parameters_table.Data(:,4));
    %
    acting_parameters = handles.parameters(1==handles.parameters_choice);
    set(handles.histo_param_1,'Value',1);
    if ~isempty(acting_parameters)
        set(handles.histo_param_1,'String',acting_parameters);
    end
    %
    guidata(hObject, handles);
    uiresume(handles.figure1); 
    
% -------------------------------------------------------------
function [med1,med2] = MCED(N_)
%MCED maximum close to each other divisors
%   
N = abs(round(N_));
med1 = [];
med2 = [];

H = floor(sqrt(N));

for k=H:-1:1
    if 0==rem(N,k)
       med1 = k;
       med2 = N/k;
       break;
    end
end


% --- Executes on button press in manova1_coords.
function manova1_coords_Callback(hObject, eventdata, handles)
% hObject    handle to manova1_coords (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of manova1_coords
