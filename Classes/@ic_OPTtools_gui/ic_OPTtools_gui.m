classdef ic_OPTtools_gui
        
    % Copyright (C) 2013 Imperial College London.
    % All rights reserved.
    %
    % This program is free software; you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation; either version 2 of the License, or
    % (at your option) any later version.
    %
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License along
    % with this program; if not, write to the Free Software Foundation, Inc.,
    % 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
    %
    % This software tool was developed with support from the UK 
    % Engineering and Physical Sciences Council 
    % through  a studentship from the Institute of Chemical Biology 
    % and The Wellcome Trust through a grant entitled 
    % "The Open Microscopy Environment: Image Informatics for Biological Sciences" (Ref: 095931).
   
    properties
        
        window = [];
        
    end
    
    methods
      
        function obj = ic_OPTtools_gui(wait,require_auth)
                                                    
            if nargin < 1
                wait = false;
            end
                
            if nargin < 2
                require_auth = false;
            end
            
            if ~isdeployed
                addpath_ALYtools;
            else
                wait = true;
            end
                                  
            profile = ic_OPTtools_profile_controller();
            profile.load_profile();
            

            % Try and read in version number
            try
                v = textread(['GeneratedFiles' filesep 'version.txt'],'%s');
                v = v{1};
            catch
                v = '[unknown version]';
            end                                              

            obj.window = figure( ...
                'Name', ['ic_OPTtools ' v], ...
                'NumberTitle', 'off', ...
                'MenuBar', 'none', ...
                'Position', [500 500 420 .1], ...                
                'Toolbar', 'none', ...
                'DockControls', 'off', ...                
                'Resize', 'off', ...                
                'HandleVisibility', 'off', ...
                'Visible','off');
            
            handles = guidata(obj.window); 
                                                        
            handles.version = v;
            handles.window = obj.window;
            handles.use_popup = true;
                                                                                   
            handles.data_controller = ic_OPTtools_data_controller(handles);                        
            handles.omero_data_manager = ic_OPTtools_omero_data_manager(handles);            
            
            handles = obj.setup_menu(handles);            
            
            handles.menu_controller = ic_OPTtools_front_end_menu_controller(handles);
                        
            guidata(obj.window,handles);
                        
            if ~isdeployed
            
                loadOmero();

                % find path to OMEuiUtils.jar - approach copied from
                % bfCheckJavaPath

                % first check it isn't already in the dynamic path
                jPath = javaclasspath('-dynamic');
                utilJarInPath = false;
                for i = 1:length(jPath)
                    if strfind(jPath{i},'OMEuiUtils.jar');
                        utilJarInPath = true;
                        break;
                    end
                end

                if ~utilJarInPath
                    path = which('OMEuiUtils.jar');
                    if isempty(path)
                        path = fullfile(fileparts(mfilename('fullpath')), 'OMEuiUtils.jar');
                    end
                    if ~isempty(path) && exist(path, 'file') == 2
                        javaaddpath(path);
                    else 
                         assert('Cannot automatically locate an OMEuiUtils JAR file');
                    end
                end

                % verify that enough memory is allocated
                bfCheckJavaMemory();

                % load both bioformats & OMERO
                autoloadBioFormats = 1;

                % load the Bio-Formats library into the MATLAB environment
                status = bfCheckJavaPath(autoloadBioFormats);
                assert(status, ['Missing Bio-Formats library. Either add loci_tools.jar '...
                    'to the static Java path or add it to the Matlab path.']);
   
            end

            % initialize logging
            loci.common.DebugTools.enableLogging('INFO');
            java.lang.System.setProperty('javax.xml.transform.TransformerFactory', 'com.sun.org.apache.xalan.internal.xsltc.trax.TransformerFactoryImpl');            
                        
            close all;
            
            set(obj.window,'Visible','on');
            set(obj.window,'CloseRequestFcn',@obj.close_request_fcn);
                        
            if wait
                waitfor(obj.window);
            end                            
            
        end
        
        function vx = split_ver(obj,ver)
            % Convert version string into a number
            tk = regexp(ver,'([0-9]+).([0-9]+).([0-9]+)','tokens');
            if ~isempty(tk{1})
                tk = tk{1};
                vx = str2double(tk{1})*1e6 + str2double(tk{2})*1e3 + str2double(tk{3});
            else 
                vx = 0;
            end
        end
        
        function close_request_fcn(obj,~,~)
            
           if isempty(obj.window)
               disp('no window..');
               return;
           end
                        
           handles = guidata(obj.window);
           
           if isempty(handles)
               disp('close_request_fcn: hadles empty');
               delete(obj.window);
               return;
           end
           
           try
           client = handles.omero_data_manager.client;            
            if ~isempty(client)                
                %
                disp('Closing OMERO session');
                client.closeSession();
                %
                 handles.omero_data_manager.session = [];
                 handles.omero_data_manager.client = [];
            end
           catch
               disp('not closing Omero..');
           end
            
           try
            % Make sure we clean up all the left over classes
            names = fieldnames(handles);
                      
            for i=1:length(names)
                % Check the field is actually a handle and isn't the window
                % which we need to close right at the end
                if ~strcmp(names{i},'window') && all(ishandle(handles.(names{i})))
                    delete(handles.(names{i}));
                end
            end
           catch
               disp('not deleting handles..');
           end
            
            % Finally actually close window
            delete(handles.window);
            
        end
        
    end
    
end
