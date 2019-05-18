classdef ALYtools_menu_controller < handle
        
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
                
        menu_OMERO_login;
        menu_OMERO_Working_Data_Info;
        menu_OMERO_Set_Dataset;
        menu_OMERO_Switch_User;
        menu_OMERO_Connect_To_Another_User;
        menu_OMERO_Connect_To_Logon_User;
        menu_OMERO_Reset_Logon; 
            %                   
            menu_OMERO_load;
            
        % holy cows
        omero_data_manager;     
        data_controller;
        
        window;
        version;

    %================================= file                
        menu_file_Data_Info;
        menu_file_new_window;
        menu_file_load;
        menu_file_load_multiple;

    %================================= segmentation    
    menu_segmentation_adjust;
    menu_segmentation_save;
    menu_segmentation_load;
    menu_segmentation_clear;
    
    %================================= analysis        
    menu_analysis_current;
    menu_analysis_batch;
        
    %================================= settings    
    menu_settings_microns_per_pixel;
    menu_settings_problem;
    menu_settings_problem_dependent;
    menu_settings_general;
    menu_settings_save;
    menu_settings_load;
    
    menu_settings_problem_FungusDependentGranuleRelease;
    menu_settings_problem_CIDR;
    menu_settings_problem_TTO; 
    menu_settings_problem_PR;
    menu_settings_problem_HL1;
    menu_settings_problem_NucCyt;
    menu_settings_problem_MPHG;  
    menu_settings_problem_Sparks;
    menu_settings_problem_Experimental; 
    menu_settings_problem_per_image_TCSPC_FLIM;
    menu_settings_problem_per_image_TCSPC_FLIM_PHASOR;
    menu_settings_problem_t_dependent_Nuclei_ratio_FRET;
    menu_settings_problem_Image_Tiling;
    menu_settings_prblm_AI_Powered_2D_SMLM_Reconstruction;        
        
    %================================= visualization    
    menu_visualization_setup_Icy_directory;
    menu_visualization_start_Icy;
    %
    menu_visualization_send_original_to_Icy;
                                
    end
    
    properties(SetObservable = true)

    end
        
    methods
        
        %------------------------------------------------------------------        
        function obj = ALYtools_menu_controller(handles)
            
            assign_handles(obj,handles);
            set_callbacks(obj);
            
            obj.data_controller.menu_controller = obj;    
            
            set(obj.window,'Name',['ALYtools ' obj.version ' : ' obj.data_controller.problem]);
                       
            % automatically start Icy
            try
                icy_figure;
                icy_closeall;
            catch                                 
                if ispc && ~isempty(obj.data_controller.IcyDirectory)
                        dos([obj.data_controller.IcyDirectory filesep 'icy']);
                    elseif ismac && ~isempty(obj.data_controller.IcyDirectory)                    
                        unix(['open ' obj.data_controller.IcyDirectory filesep 'icy']); % ?                    
                    elseif isunix
                        unix('icy');                                                
                else
                    % msgbox('error - Icy directory was not set up');
                end 
            end
            %  automatically start Icy
        end
        %------------------------------------------------------------------
        function set_callbacks(obj)
            
             mc = metaclass(obj);
             obj_prop = mc.Properties;
             obj_method = mc.Methods;
                          
             % Search for properties with corresponding callbacks
             for i=1:length(obj_prop)
                prop = obj_prop{i}.Name;
                if strncmp(prop,'menu_',5)
                    method = [prop '_callback'];
                    matching_methods = findobj([obj_method{:}],'Name',method);
                    if ~isempty(matching_methods)               
                        eval(['set(obj.' prop ',''Callback'',@obj.' method ')' ]);
                    end
                end          
             end
             
        end
        %------------------------------------------------------------------                                       
        function menu_file_new_window_callback(obj,~,~)
            ALYtools();
        end
                        
        %------------------------------------------------------------------
        % OMERO
        %------------------------------------------------------------------
        function menu_OMERO_login_callback(obj,~,~)
            obj.omero_data_manager.Omero_logon();
            
            if ~isempty(obj.omero_data_manager.session)
                props = properties(obj);
                OMERO_props = props( strncmp('menu_OMERO',props,10) );
                for i=1:length(OMERO_props)
                    set(obj.(OMERO_props{i}),'Enable','on');
                end
            end            
        end
        %------------------------------------------------------------------
        function menu_OMERO_Set_Dataset_callback(obj,~,~)            
            infostring = obj.omero_data_manager.Set_Dataset();
            if ~isempty(infostring)
                set(obj.menu_OMERO_Working_Data_Info,'Label',infostring,'ForegroundColor','blue');
                set(obj.menu_Batch_Indicator_Src,'Label',infostring,'ForegroundColor','blue');                
            end;
        end                        
        %------------------------------------------------------------------        
        function menu_OMERO_Reset_Logon_callback(obj,~,~)
            obj.omero_data_manager.Omero_logon();
        end
        %------------------------------------------------------------------        
        function menu_OMERO_Switch_User_callback(obj,~,~)
            %delete([ pwd '\' obj.omero_data_manager.omero_logon_filename ]);
            obj.omero_data_manager.Omero_logon_forced();
        end        
        %------------------------------------------------------------------
        function menu_OMERO_Connect_To_Another_User_callback(obj,~,~)
            obj.omero_data_manager.Select_Another_User();
            obj.omero_data_manager.dataset = [];
            obj.data_controller.proj = [];
            obj.data_controller.volm = [];            
            obj.data_controller.on_proj_and_volm_clear;
            set(obj.menu_OMERO_Working_Data_Info,'Label','...','ForegroundColor','red');
        end                            
        %------------------------------------------------------------------
        function menu_OMERO_Connect_To_Logon_User_callback(obj,~,~)            
            obj.omero_data_manager.userid = obj.omero_data_manager.session.getAdminService().getEventContext().userId;
            obj.omero_data_manager.dataset = [];
            obj.data_controller.proj = [];
            obj.data_controller.volm = [];            
            obj.data_controller.on_proj_and_volm_clear;
            set(obj.menu_OMERO_Working_Data_Info,'Label','...','ForegroundColor','red');
        end  
         %------------------------------------------------------------------
        function menu_OMERO_set_single_callback(obj, ~, ~)                                               
            infostring = obj.data_controller.OMERO_load_single(obj.omero_data_manager,true); % verbose
            if ~isempty(infostring)
                set(obj.menu_OMERO_Working_Data_Info,'Label',infostring,'ForegroundColor','blue','Enable','on');
                set(obj.menu_file_Working_Data_Info,'Label','...','Enable','off');                
                obj.data_controller.current_filename = [];
                set(obj.menu_settings_Zrange,'Label','Z range : full');
                obj.data_controller.Z_range = []; % no selection                                    
            end;            
        end
         %------------------------------------------------------------------
        function menu_OMERO_set_multiple_callback(obj, ~, ~)
            obj.data_controller.OMERO_load_multiple(obj.omero_data_manager);
        end
        
        %------------------------------------------------------------------
        % OMERO
        %------------------------------------------------------------------                                
                                        
         %------------------------------------------------------------------        
        function menu_tools_preferences_callback(obj,~,~)
            profile = ALYtools_profile_controller();
            profile.set_profile();
        end        
         %------------------------------------------------------------------                
        
       function menu_file_load_callback(obj, ~, ~)
           try
            [filename,pathname] = uigetfile({'*.tif;*.tiff;*.gif;*.avi;*.czi;*.png;*.bmp;*.jpg;*.jpeg;*.lsm;*.sdt;*.nd2','Image Files'}, ...
                'Select data file',obj.data_controller.DefaultDirectory);
           catch
            [filename,pathname] = uigetfile({'*.tif;*.tiff;*.gif;*.avi;*.czi;*.png;*.bmp;*.jpg;*.jpeg;*.lsm;*.sdt;*.nd2','Image Files'}, ...
                'Select data file',pwd);               
           end
            if filename == 0, return, end;       
            %
            obj.data_controller.load(filename,pathname,true); % verbose true
            set(obj.menu_settings_microns_per_pixel,'Label',['Microns per pixel ' num2str(obj.data_controller.microns_per_pixel)]);
        end        
         %------------------------------------------------------------------                

        function menu_file_load_multiple_callback(obj, ~, ~)
            try
                [filenames,pathname] = uigetfile({'*.tif;*.tiff;*.gif;*.avi;*.czi;*.png;*.bmp;*.jpg;*.jpeg;*.lsm;*.sdt;*.nd2;*.txt','Image Files'},...
                'Select data files',obj.data_controller.DefaultDirectory,'MultiSelect','on');
            catch
                [filenames,pathname] = uigetfile({'*.tif;*.tiff;*.gif;*.avi;*.czi;*.png;*.bmp;*.jpg;*.jpeg;*.lsm;*.sdt;*.nd2;*.txt','Image Files'},...
                'Select data files',pwd,'MultiSelect','on');                
            end
            if isempty(filenames), return, end;       
            if isnumeric(filenames) && 0==filenames, return, end;
                        
            %obj.data_controller.load_multiple(filenames,pathname,true); % verbose true
            obj.data_controller.load_multiple(filenames,pathname,false); % verbose true
            set(obj.menu_settings_microns_per_pixel,'Label',['Microns per pixel ' num2str(obj.data_controller.microns_per_pixel)]);
        end                 %------------------------------------------------------------------                
         
         
    %================================= segmentation
    
    function menu_segmentation_adjust_callback(obj, ~, ~)
        obj.data_controller.adjust_segmentation;
    end        
    %------------------------------------------------------------------                    
    function menu_segmentation_save_callback(obj, ~, ~)
        obj.data_controller.save_segmentation;
    end        
    %------------------------------------------------------------------                    
    function menu_segmentation_load_callback(obj, ~, ~)
        obj.data_controller.load_segmentation;
    end            
    %------------------------------------------------------------------
    function menu_segmentation_clear_callback(obj, ~, ~)
        obj.data_controller.M_sgm = [];
    end            
    %------------------------------------------------------------------    
    
    %================================= analysis    
    
    function menu_analysis_current_callback(obj, ~, ~)
        obj.data_controller.analyze_current;
    end
    
    function menu_analysis_batch_callback(obj, ~, ~)
        verbose = true;
        % get list of files...
        try
            [filenames, pathname] = uigetfile({'*.tif;*.tiff;*.gif;*.avi;*.czi;*.png;*.bmp;*.sdt;*.jpg;*.jpeg;*.lsm','Image Files'},'Select data files', ...
            obj.data_controller.DefaultDirectory,'MultiSelect','on');
        catch
            [filenames, pathname] = uigetfile({'*.tif;*.tiff;*.gif;*.avi;*.czi;*.png;*.bmp;*.jpg;*.jpeg;*.lsm','Image Files'},'Select data files', ...
            pwd,'MultiSelect','on');            
        end
        if pathname == 0, return, end;
        %
        obj.data_controller.run_batch(filenames, pathname, verbose);
    end
         %------------------------------------------------------------------                    
        
    %================================= settings
        
    function menu_settings_problem_FungusDependentGranuleRelease_callback(obj, ~, ~) 
         set(obj.menu_settings_problem,'Label',['Problem = ' 'Fungus Dependent Granule Release']);
         obj.data_controller.problem = 'Fungus Dependent Granule Release';
         set(obj.window,'Name',['ALYtools ' obj.version ' : ' obj.data_controller.problem]);
    end
    function menu_settings_problem_CIDR_callback(obj, ~, ~) 
         set(obj.menu_settings_problem,'Label',['Problem = ' 'CIDR']);        
         obj.data_controller.problem = 'CIDR';
         set(obj.window,'Name',['ALYtools ' obj.version ' : ' obj.data_controller.problem]);         
    end
    function menu_settings_problem_TTO_callback(obj, ~, ~) 
         set(obj.menu_settings_problem,'Label',['Problem = ' 'TTO']);        
         obj.data_controller.problem = 'TTO';
         set(obj.window,'Name',['ALYtools ' obj.version ' : ' obj.data_controller.problem]);         
    end
    function menu_settings_problem_PR_callback(obj, ~, ~) 
         set(obj.menu_settings_problem,'Label',['Problem = ' 'PR']);        
         obj.data_controller.problem = 'PR';
         set(obj.window,'Name',['ALYtools ' obj.version ' : ' obj.data_controller.problem]);         
    end
    function menu_settings_problem_HL1_callback(obj, ~, ~) 
         set(obj.menu_settings_problem,'Label',['Problem = ' 'HL1']);        
         obj.data_controller.problem = 'HL1';
         set(obj.window,'Name',['ALYtools ' obj.version ' : ' obj.data_controller.problem]);         
    end        
    function menu_settings_problem_NucCyt_callback(obj, ~, ~) 
         set(obj.menu_settings_problem,'Label',['Problem = ' 'NucCyt']);        
         obj.data_controller.problem = 'NucCyt'; 
         set(obj.window,'Name',['ALYtools ' obj.version ' : ' obj.data_controller.problem]);         
    end            
    function menu_settings_problem_MPHG_callback(obj, ~, ~) 
         set(obj.menu_settings_problem,'Label',['Problem = ' 'MPHG']);        
         obj.data_controller.problem = 'MPHG';
         set(obj.window,'Name',['ALYtools ' obj.version ' : ' obj.data_controller.problem]);         
    end                
    function menu_settings_problem_Sparks_callback(obj, ~, ~) 
         set(obj.menu_settings_problem,'Label',['Problem = ' 'Sparks']);        
         obj.data_controller.problem = 'Sparks';
         set(obj.window,'Name',['ALYtools ' obj.version ' : ' obj.data_controller.problem]);         
    end                    
    function menu_settings_problem_Experimental_callback(obj, ~, ~) 
         set(obj.menu_settings_problem,'Label',['Problem = ' 'Experimental']);        
         obj.data_controller.problem = 'Experimental';
         set(obj.window,'Name',['ALYtools ' obj.version ' : ' obj.data_controller.problem]);         
    end            
    function menu_settings_problem_per_image_TCSPC_FLIM_callback(obj, ~, ~) 
         set(obj.menu_settings_problem,'Label',['Problem = ' 'per_image_TCSPC_FLIM']);        
         obj.data_controller.problem = 'per_image_TCSPC_FLIM';
         set(obj.window,'Name',['ALYtools ' obj.version ' : ' obj.data_controller.problem]);         
    end 
    function menu_settings_problem_per_image_TCSPC_FLIM_PHASOR_callback(obj, ~, ~) 
         set(obj.menu_settings_problem,'Label',['Problem = ' 'per_image_TCSPC_FLIM_PHASOR']);        
         obj.data_controller.problem = 'per_image_TCSPC_FLIM_PHASOR';
         set(obj.window,'Name',['ALYtools ' obj.version ' : ' obj.data_controller.problem]);         
    end         
     function menu_settings_problem_t_dependent_Nuclei_ratio_FRET_callback(obj, ~, ~) 
         if isempty(strfind(version('-java'),'1.8')), errordlg('please switch to Java 1.8 to enable TrackMate - can not continue'); end %#ok<*CPROPLC>
         set(obj.menu_settings_problem,'Label',['Problem = ' 't_dependent_Nuclei_ratio_FRET']);        
         obj.data_controller.problem = 't_dependent_Nuclei_ratio_FRET';
         set(obj.window,'Name',['ALYtools ' obj.version ' : ' obj.data_controller.problem]);         
     end                        
    function menu_settings_problem_Image_Tiling_callback(obj, ~, ~) 
         set(obj.menu_settings_problem,'Label',['Problem = ' 'Image_Tiling']);        
         obj.data_controller.problem = 'Image_Tiling';
         set(obj.window,'Name',['ALYtools ' obj.version ' : ' obj.data_controller.problem]);         
    end             
    function menu_settings_prblm_AI_Powered_2D_SMLM_Reconstruction_callback(obj, ~, ~) 
         set(obj.menu_settings_problem,'Label',['Problem = ' 'AI_Powered_2D_SMLM_Reconstruction']);        
         obj.data_controller.problem = 'AI_Powered_2D_SMLM_Reconstruction';
         set(obj.window,'Name',['ALYtools ' obj.version ' : ' obj.data_controller.problem]);         
    end
           
    %------------------------------------------------------------------                        
    function menu_settings_problem_dependent_callback(obj, ~, ~)
        
        if strcmp(obj.data_controller.problem,'Fungus Dependent Granule Release')
             Granule_Fungus_Cells_Problem_Specific_settings(obj.data_controller);
        elseif strcmp(obj.data_controller.problem,'TTO')
             TTO_Problem_Specific_settings(obj.data_controller);
        elseif strcmp(obj.data_controller.problem,'CIDR')
             CIDR_Problem_Specific_settings(obj.data_controller);
        elseif strcmp(obj.data_controller.problem,'PR')
             PR_Problem_Specific_settings(obj.data_controller); 
        elseif strcmp(obj.data_controller.problem,'HL1')
             HL1_Problem_Specific_settings(obj.data_controller);                 
        elseif strcmp(obj.data_controller.problem,'NucCyt')
             NucCyt_Problem_Specific_settings(obj.data_controller);                              
        elseif strcmp(obj.data_controller.problem,'MPHG')
             MPHG_Problem_Specific_settings(obj.data_controller);                                           
        elseif strcmp(obj.data_controller.problem,'per_image_TCSPC_FLIM') || strcmp(obj.data_controller.problem,'per_image_TCSPC_FLIM_PHASOR') 
             per_image_TCSPC_FLIM_Problem_Specific_settings(obj.data_controller);                                                        
        elseif strcmp(obj.data_controller.problem,'Image_Tiling')
             ImageTiling_Problem_Specific_settings(obj.data_controller);                                                                
        elseif strcmp(obj.data_controller.problem,'AI_Powered_2D_SMLM_Reconstruction')
             AI_Powered_2d_SMLM_reconstruction_settings(obj.data_controller);                                                        
        end
                                
    end        
         %------------------------------------------------------------------                             
         function menu_settings_general_callback(obj, ~, ~)
             ALYtools_General_Settings(obj.data_controller);
         end
        
    %================================= visualization
    
         %------------------------------------------------------------------                
        function menu_visualization_setup_Icy_directory_callback(obj, ~,~)
            [path] = uigetdir(obj.data_controller.DefaultDirectory,'Guide to Icy directory');
            if path ~= 0
                obj.data_controller.IcyDirectory = path;                
            end                       
        end
         %------------------------------------------------------------------                
        function menu_visualization_start_Icy_callback(obj, ~,~)
            
            if ~isempty(obj.data_controller.IcyDirectory)                
                if ispc
                    dos([obj.data_controller.IcyDirectory filesep 'icy']);
                elseif ismac                    
                    unix(['open ' obj.data_controller.IcyDirectory filesep 'icy']); % ?                    
                elseif isunix
                    % ?
                end                                
            else
                msgbox('error - Icy directory was not set up');
            end
            
        end
        
        %------------------------------------------------------------------          
    function menu_visualization_send_original_to_Icy_callback(obj, ~, ~)
        obj.data_controller.send_original_to_Icy;
    end
                
         %------------------------------------------------------------------                
            
        function menu_settings_microns_per_pixel_callback(obj, ~, ~)       
             value = enter_value();
             if isempty(value) || ~isnumeric(value) || value==obj.data_controller.microns_per_pixel || value <=0, return, end;
             obj.data_controller.microns_per_pixel = value;            
             set(obj.menu_settings_microns_per_pixel,'Label',['Microns per pixel ' num2str(obj.data_controller.microns_per_pixel)]);
        end        
        
         %------------------------------------------------------------------                                                                           
         function menu_settings_save_callback(obj, ~, ~)                    
            [fname, fpath] = uiputfile('*.xml','Save Settings as..',[obj.data_controller.DefaultDirectory filesep 'ALYtools_settings']);
            if fpath == 0; return; end
            filespec = fullfile(fpath,fname);
            try
                obj.data_controller.save_settings(filespec);
            catch
                errordlg('Error while trying to save settings file, ehm.. write protection?');
            end
         end
         %------------------------------------------------------------------                         
         function menu_settings_load_callback(obj, ~, ~)                    
             [fname, fpath] = uigetfile({'*.xml'},'Select settings file',obj.data_controller.DefaultDirectory);
             if fpath == 0, return, end;
             filespec = fullfile(fpath,fname);
             try
                obj.data_controller.load_settings(filespec);
                set(obj.menu_settings_microns_per_pixel,'Label',['Microns per pixel ' num2str(obj.data_controller.microns_per_pixel)]);
             catch
                errordlg('Error while trying to load settings file');
             end
         end
                  
    %================================= % VANITY       
    
        %------------------------------------------------------------------
        function menu_help_about_callback(obj, ~, ~)
            % to do
        end            
        %------------------------------------------------------------------
        function menu_help_tracker_callback(obj, ~, ~)
            % to do
        end            
        %------------------------------------------------------------------
        function menu_help_bugs_callback(obj, ~, ~)
            % to do
        end
                            
    end
    
end
