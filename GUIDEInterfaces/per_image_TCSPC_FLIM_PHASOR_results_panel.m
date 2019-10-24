function varargout = per_image_TCSPC_FLIM_PHASOR_results_panel(varargin)
% PER_IMAGE_TCSPC_FLIM_PHASOR_RESULTS_PANEL MATLAB code for per_image_TCSPC_FLIM_PHASOR_results_panel.fig
%      PER_IMAGE_TCSPC_FLIM_PHASOR_RESULTS_PANEL, by itself, creates a new PER_IMAGE_TCSPC_FLIM_PHASOR_RESULTS_PANEL or raises the existing
%      singleton*.
%
%      H = PER_IMAGE_TCSPC_FLIM_PHASOR_RESULTS_PANEL returns the handle to a new PER_IMAGE_TCSPC_FLIM_PHASOR_RESULTS_PANEL or the handle to
%      the existing singleton*.
%
%      PER_IMAGE_TCSPC_FLIM_PHASOR_RESULTS_PANEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PER_IMAGE_TCSPC_FLIM_PHASOR_RESULTS_PANEL.M with the given input arguments.
%
%      PER_IMAGE_TCSPC_FLIM_PHASOR_RESULTS_PANEL('Property','Value',...) creates a new PER_IMAGE_TCSPC_FLIM_PHASOR_RESULTS_PANEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before per_image_TCSPC_FLIM_PHASOR_results_panel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to per_image_TCSPC_FLIM_PHASOR_results_panel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help per_image_TCSPC_FLIM_PHASOR_results_panel

% Last Modified by GUIDE v2.5 31-May-2018 15:49:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @per_image_TCSPC_FLIM_PHASOR_results_panel_OpeningFcn, ...
                   'gui_OutputFcn',  @per_image_TCSPC_FLIM_PHASOR_results_panel_OutputFcn, ...
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


% --- Executes just before per_image_TCSPC_FLIM_PHASOR_results_panel is made visible.
function per_image_TCSPC_FLIM_PHASOR_results_panel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to per_image_TCSPC_FLIM_PHASOR_results_panel (see VARARGIN)

res = varargin{1};

cla(handles.phasor_plot,'reset');

                phasors = res{1};
                Tp = res{2};
                histogram2(handles.phasor_plot,phasors(:,1),phasors(:,2),'Normalization','probability','DisplayStyle','tile','EdgeColor','none','LineStyle','none');
                
                hold on;
                theta = 0:pi/1000:pi;
                plot(handles.phasor_plot,(cos(theta)+1)/2,sin(theta)/2,'r-','linewidth',2);
                
                hold on;
                W = 2*pi/Tp;
                
                %tau = 0:500:20000;
                npoints = 50;
                dtau = round(Tp/npoints);
                tau = 0:dtau:Tp;
                
                G = 1./(1+(W*tau).^2);
                S = W*tau.*G;
                plot(handles.phasor_plot,G,S,'m.','markersize',30);
                
                hold on;
                for k=1:numel(tau)
                text(handles.phasor_plot,G(k),S(k),['\bf\it' num2str(tau(k)/1000)],'fontsize',16,'Color','Black');
                hold on;
                end
                                             
                hold off;
                axis([-.05 1 0 0.7]);
                xlabel('G');
                ylabel('S');
                grid on
                set(handles.phasor_plot,'DataAspectRatio', [1 1 1])
                %   

% Choose default command line output for per_image_TCSPC_FLIM_PHASOR_results_panel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes per_image_TCSPC_FLIM_PHASOR_results_panel wait for user response (see UIRESUME)
% uiwait(handles.per_image_TCSPC_FLIM_PHASOR_results);


% --- Outputs from this function are returned to the command line.
function varargout = per_image_TCSPC_FLIM_PHASOR_results_panel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Copy_to_clipboard_Callback(hObject, eventdata, handles)
% hObject    handle to Copy_to_clipboard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

print('-clipboard','-dmeta','-noui');
