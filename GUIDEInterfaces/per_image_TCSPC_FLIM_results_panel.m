function varargout = per_image_TCSPC_FLIM_results_panel(varargin)
% PER_IMAGE_TCSPC_FLIM_RESULTS_PANEL MATLAB code for per_image_TCSPC_FLIM_results_panel.fig
%      PER_IMAGE_TCSPC_FLIM_RESULTS_PANEL, by itself, creates a new PER_IMAGE_TCSPC_FLIM_RESULTS_PANEL or raises the existing
%      singleton*.
%
%      H = PER_IMAGE_TCSPC_FLIM_RESULTS_PANEL returns the handle to a new PER_IMAGE_TCSPC_FLIM_RESULTS_PANEL or the handle to
%      the existing singleton*.
%
%      PER_IMAGE_TCSPC_FLIM_RESULTS_PANEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PER_IMAGE_TCSPC_FLIM_RESULTS_PANEL.M with the given input arguments.
%
%      PER_IMAGE_TCSPC_FLIM_RESULTS_PANEL('Property','Value',...) creates a new PER_IMAGE_TCSPC_FLIM_RESULTS_PANEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before per_image_TCSPC_FLIM_results_panel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to per_image_TCSPC_FLIM_results_panel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help per_image_TCSPC_FLIM_results_panel

% Last Modified by GUIDE v2.5 14-Aug-2017 10:17:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @per_image_TCSPC_FLIM_results_panel_OpeningFcn, ...
                   'gui_OutputFcn',  @per_image_TCSPC_FLIM_results_panel_OutputFcn, ...
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


% --- Executes just before per_image_TCSPC_FLIM_results_panel is made visible.
function per_image_TCSPC_FLIM_results_panel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to per_image_TCSPC_FLIM_results_panel (see VARARGIN)

res = varargin{1};

cla(handles.decays,'reset');
cla(handles.residuals,'reset');

                t = res.t;
                y = res.y;
                y_est = res.y_est;
                fitting_mask = res.fitting_mask;
                data_controller =res.data_controller;
                Ndecays = res.Ndecays;
                Nbins = res.Nbins;
                datas = res.datas;
                captions = res.captions;
                table_names = res.table_names;
                IRF = res.IRF;
                
                linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':',...
                '-.','--','-',':','-.','--','-',':','-.'));
                MarkerEdgeColors=jet(Ndecays);  % n is the number of different items you have
                Markers=['o','x','+','*','s','d','v','^','<','>','p','h','.',...
                '+','*','o','x','^','<','h','.','>','p','s','d','v',...
                'o','x','+','*','s','d','v','^','<','>','p','h','.']; 
                %
                while numel(linestyles) < Ndecays
                    linestyles = [linestyles; linestyles];
                end
                while numel(Markers) < Ndecays
                    Markers = [Markers Markers];
                end                
                %
                handles.linestyles = linestyles;
                handles.MarkerEdgeColors = MarkerEdgeColors;
                handles.Markers = Markers;
                handles.t = res.t;
                handles.IRF = res.IRF;
                handles.data_controller = data_controller;
            
                DECAYS = zeros(Nbins,Ndecays);
                DECAYS_THEOR = zeros(Nbins,Ndecays);
                RESIDUALS = zeros(Nbins,Ndecays);
    
                r = y-y_est';
                norm_res = r./sqrt(y_est');
                
                for k=1:Ndecays
                    RANGE = ((k-1)*Nbins+1):(k*Nbins);
                    DECAYS(:,k) = y(RANGE).*fitting_mask; 
                    DECAYS_THEOR(:,k) = y_est(RANGE).*fitting_mask';
                    RESIDUALS(:,k) = - norm_res(RANGE).*fitting_mask; % to make it similar to FLIMfit
                end                 
                %                
                min_decays = min(min(y(y>0)),min(y_est(y_est>0)));
                max_decays = max(max(y(y>0)),max(y_est(y_est>0)));
                min_resdls = min(norm_res);
                max_resdls = max(norm_res);                
                
                handles.RAW_DECAYS = res.raw_decays;
                handles.DECAYS = DECAYS;
                handles.DECAYS_THEOR = DECAYS_THEOR; 
                handles.RESIDUALS = RESIDUALS;
                handles.bckg_corrected_decays = res.bckg_corrected_decays;
                
                handles.min_decays = min_decays;
                handles.max_decays = max_decays;
                handles.min_resdls = min_resdls;
                handles.max_resdls = max_resdls;
                handles.datas = res.datas;
                handles.table_names = res.table_names;
                handles.captions = res.captions;
                
                set(handles.results,'CellEditCallback',@decay_check_callback);
                
                L = size(datas,1);
                set(handles.results, 'Data', [num2cell(true(L,1)) datas]);                    
                set(handles.results,'ColumnName',[' ' captions]);
                
% Choose default command line output for per_image_TCSPC_FLIM_results_panel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

visualize(handles);

% UIWAIT makes per_image_TCSPC_FLIM_results_panel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to thecommand line.
function varargout = per_image_TCSPC_FLIM_results_panel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in display_IRF.
function display_IRF_Callback(hObject, eventdata, handles)
% hObject    handle to display_IRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    visualize(handles);

%-------------------------------------------------------------------------%    
function decay_check_callback(hObject,callbackdata)

    handles = guidata(hObject); 
    visualize(handles);
       
%-------------------------------------------------------------------------%      
function visualize(handles)

     cla(handles.decays,'reset');
     cla(handles.residuals,'reset');

     display_flag = cell2mat(handles.results.Data(:,1));
     
     display_IRF = get(handles.display_IRF,'Value');
    
     DECAYS= handles.DECAYS;
     DECAYS_THEOR= handles.DECAYS_THEOR; 
     RESIDUALS= handles.RESIDUALS;
                
     min_decays= handles.min_decays;
     max_decays= handles.max_decays;
     min_resdls= handles.min_resdls;
     max_resdls= handles.max_resdls;
                
     linestyles= handles.linestyles;
     MarkerEdgeColors= handles.MarkerEdgeColors;
     Markers= handles.Markers;
     t = handles.t;
     IRF = handles.IRF;
     data_controller = handles.data_controller;
     
                Ndecays = size(DECAYS,2);
                for k=1:Ndecays
                    if display_flag(k)
                        semilogy(handles.decays,t,DECAYS(:,k),[linestyles{k} Markers(k)],'Color',MarkerEdgeColors(k,:));
                        hold(handles.decays,'on'); 
                        semilogy(handles.decays,t,DECAYS_THEOR(:,k),'k.-');                              
                        %
                        %plot(handles.residuals,t,RESIDUALS(:,k),[linestyles{k} Markers(k)],'Color',MarkerEdgeColors(k,:));
                        plot(handles.residuals,t,RESIDUALS(:,k),'Color',MarkerEdgeColors(k,:));
                        hold(handles.residuals,'on');
                    end
                end
                
                if display_IRF
                        hold(handles.decays,'on');
                        coef = max_decays/max(IRF);
                        semilogy(handles.decays,t,coef*IRF,'m-');                         
                end
                
                ylabel(handles.decays,'fluorescence intensity'); 
                grid(handles.decays,'on');
                axis(handles.decays,[min(t) max(t) real(min_decays), real(max_decays)]);
                %
                ylabel(handles.residuals,'normalized residuals');
                xlabel(handles.residuals,'time [ps]');
                grid(handles.residuals,'on');
                axis(handles.residuals,[min(t) max(t) real(min_resdls) real(max_resdls)]);
                               
                % HACK.. that I suspect is needed is because of figure-brackets usage...
                LEGEND = [];
                if 1==numel(data_controller.M_filenames)
                    LEGEND = {data_controller.M_filenames{k}, 'theor'};
                else
                    for k=1:numel(data_controller.M_filenames)
                        if display_flag(k)
                            LEGEND = [LEGEND; data_controller.M_filenames{k}; 'theor'];
                        end
                    end;
                end 
                legend(handles.decays,LEGEND);
                title(handles.decays,char(data_controller.per_image_TCSPC_FLIM_fit_model));
        
        
% --- Executes on button press in show_decays_all.
function show_decays_all_Callback(hObject, eventdata, handles)
% hObject    handle to show_decays_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = get(handles.results,'Data');
L = size(data,1);
data(:,1) = num2cell(true(L,1));
set(handles.results, 'Data', data);
%
guidata(hObject, handles);
visualize(handles);

% --- Executes on button press in show_decays_none.
function show_decays_none_Callback(hObject, eventdata, handles)
% hObject    handle to show_decays_none (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = get(handles.results,'Data');
L = size(data,1);
data(:,1) = num2cell(false(L,1));
set(handles.results, 'Data', data);
%
guidata(hObject, handles);
visualize(handles);

function Untitled_1_Callback(hObject, eventdata, handles)
% dummy

% --------------------------------------------------------------------
function save_current_decays_Callback(hObject, eventdata, handles)
% hObject    handle to save_current_decays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;

[fname, fpath] = uiputfile('*.csv','Save decays as..',[dc.DefaultDirectory filesep 'decays']);
            if fpath == 0; return; end
            filespec = fullfile(fpath,fname);
            try
                DECAYS = handles.RAW_DECAYS;
                DECAYS_THEOR = handles.DECAYS_THEOR; % careful - results of fitting of tvb/b corrected data
                datas = handles.datas;
                %
                Nphots = datas(:,size(datas,2)-1);
                %truncate
                Nphots = Nphots(1:numel(Nphots)-1);
                % exporting full theoretical model
                DECAYS_THEOR_PPIX = DECAYS_THEOR./cell2mat(Nphots') + ... 
                    dc.per_image_TCSPC_FLIM_tvb_scaling*dc.per_image_TCSPC_FLIM_tvb + ...
                    dc.per_image_TCSPC_FLIM_background_value;
                %
                t = handles.t';
                caption = {'t'};
                Ndecays = numel(dc.M_filenames);
                for k=1:Ndecays
                    caption = [caption dc.M_filenames{k}];
                end
                for k=1:Ndecays
                    caption = [caption [char(dc.M_filenames{k}) ' theor']];
                end                
                for k=1:Ndecays
                    caption = [caption [char(dc.M_filenames{k}) ' theorppix']];
                end                                
                data = num2cell([t DECAYS DECAYS_THEOR DECAYS_THEOR_PPIX]);
                data = [caption; data];
                xlswrite(filespec,data);                
            catch
                errordlg('Error while trying to save decays, ehm.. write protection?');
            end


% --------------------------------------------------------------------
function save_current_decays_separately_Callback(hObject, eventdata, handles)
save_current_decays_separately(handles,'raw');
                                                                            
% --------------------------------------------------------------------
function save_fitting_results_Callback(hObject, eventdata, handles)
% hObject    handle to save_fitting_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dc = handles.data_controller;
[fname, fpath] = uiputfile('*.xls','Save fitting results as..',[dc.DefaultDirectory filesep 'per_image_TCSPC_FLIM_results']);
            if fpath == 0; return; end
            filespec = fullfile(fpath,fname);
xlswrite( filespec,[handles.captions; handles.datas],char(handles.table_names) );

% --------------------------------------------------------------------
function save_decays_separately_bckg_corrected_Callback(hObject, eventdata, handles)
save_current_decays_separately(handles,'bckg corrected');


function save_current_decays_separately(handles,mode)
% MOCK THIS CAPTIONS TO BE LOADED IN FLIMfit
% Well	490.000000
% WellIndex	0.000000
% Wavelength	490.000000
% Polarisation	2.000000
% AcquisitionStatus	1.000000
% PowerLevel	4608.000000
% AcquisitionTime	13.109596

dc = handles.data_controller;

    directoryname = uigetdir(dc.DefaultDirectory);
    if isnumeric(directoryname) || isempty(directoryname), return, end;

    t = handles.t';
    
    switch mode
        case 'raw'
            DECAYS = handles.RAW_DECAYS;
        case 'bckg corrected'
            DECAYS = handles.bckg_corrected_decays;
    end
        
    IRF = dc.per_image_TCSPC_FLIM_irf; % original
    
                   if isempty(dc.M_imgdata), return, end;
                    try
                        hw = waitbar(0,'..saving segmntations..');
                        for k=1:numel(dc.M_imgdata)
                            fname = char(dc.M_filenames{k});
                            fname = strrep(fname,'.sdt','');
                            fname = strrep(fname,'.ome.tiff','');
                            fname = strrep(fname,'.OME.tiff','');
                            fname = [ directoryname filesep fname '.txt'];
                            %
                                   fid = fopen(fname,'w'); % this train is needed for it to be open in (one specific version of) FLIMfit
                                   fprintf(fid,'%s  %f\r\n','Well',490);
                                   fprintf(fid,'%s  %f\r\n','WellIndex',0);
                                   fprintf(fid,'%s  %f\r\n','Wavelength',490);
                                   fprintf(fid,'%s  %f\r\n','Polarization',2);
                                   fprintf(fid,'%s  %f\r\n','AcquisitionStatus',1);
                                   fprintf(fid,'%s  %f\r\n','PowerLevel',4608);
                                   fprintf(fid,'%s  %f\r\n','AcquisitionTime',13.109596);
                                   %
                                   for m=1:numel(t)
                                    fprintf(fid,'%f  %f\r\n',t(m),DECAYS(m,k));
                                   end
                                   fclose(fid);
                            %
                            if ~isempty(hw), waitbar(k/numel(dc.M_imgdata),hw); drawnow, end;
                        end
                        if ~isempty(hw), delete(hw), drawnow; end;
                        %
                        %irf
                                   fid = fopen([ directoryname filesep 'irf.txt'],'w');
                                   fprintf(fid,'%s  %f\r\n','Well',490);
                                   fprintf(fid,'%s  %f\r\n','WellIndex',0);
                                   fprintf(fid,'%s  %f\r\n','Wavelength',490);
                                   fprintf(fid,'%s  %f\r\n','Polarization',2);
                                   fprintf(fid,'%s  %f\r\n','AcquisitionStatus',1);
                                   fprintf(fid,'%s  %f\r\n','PowerLevel',4608);
                                   fprintf(fid,'%s  %f\r\n','AcquisitionTime',13.109596);
                                   %
                                   for m=1:numel(t)
                                    fprintf(fid,'%f  %f\r\n',t(m),IRF(m));
                                   end
                                   fclose(fid);                        
                        %irf
                        %                        
                    catch
                    end
