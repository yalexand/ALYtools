function handles = setup_menu(obj,handles)

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
           
        
    microns_per_pixel = handles.data_controller.microns_per_pixel;
    problem = handles.data_controller.problem;
        
    %================================= file

    menu_file = uimenu(obj.window,'Label','File');
    
    handles.menu_file_Data_Info = uimenu(menu_file,'Label','...','ForegroundColor','red','Enable','off');    
    handles.menu_file_new_window = uimenu(menu_file,'Label','New Window','Accelerator','N');
    handles.menu_file_load = uimenu(menu_file,'Label','Load Single','Separator','on');
    handles.menu_file_load_multiple = uimenu(menu_file,'Label','Load Multiple','Separator','on');
            
    %================================= OMERO
    
    menu_OMERO = uimenu(obj.window,'Label','OMERO');

    handles.menu_OMERO_Working_Data_Info = uimenu(menu_OMERO,'Label','...','ForegroundColor','red','Enable','off');
    handles.menu_OMERO_login = uimenu(menu_OMERO,'Label','Log in to OMERO');    

    menu_OMERO_Set_Data = uimenu(menu_OMERO,'Label','Set User');
    handles.menu_OMERO_Switch_User = uimenu(menu_OMERO_Set_Data,'Label','Switch User...','Separator','on','Enable','off');    

    handles.menu_OMERO_Connect_To_Another_User = uimenu(menu_OMERO_Set_Data,'Label','Connect to another User...','Enable','off');    
    handles.menu_OMERO_Connect_To_Logon_User = uimenu(menu_OMERO_Set_Data,'Label','Connect to Logon User...','Enable','off');    

    handles.menu_OMERO_Reset_Logon = uimenu(menu_OMERO_Set_Data,'Label','Restore Logon','Separator','on','Enable','off');                
    
    handles.menu_OMERO_load = uimenu(menu_OMERO,'Label','Load Data','Separator','on','Enable','off');
    
    %================================= segmentation
    
    menu_segmentation = uimenu(obj.window,'Label','Segmentation');
    handles.menu_segmentation_adjust = uimenu(menu_segmentation,'Label','Adjust');
    handles.menu_segmentation_save = uimenu(menu_segmentation,'Label','Make and save','Separator','on');
    handles.menu_segmentation_load = uimenu(menu_segmentation,'Label','Load and use');
    handles.menu_segmentation_clear = uimenu(menu_segmentation,'Label','Clear loaded');    
    
    %================================= analysis    
    
    menu_analysis = uimenu(obj.window,'Label','Analysis');   
    handles.menu_analysis_current = uimenu(menu_analysis,'Label','Run Current');       
    handles.menu_analysis_batch = uimenu(menu_analysis,'Label','Set Batch and Run');           
        
    %================================= settings
    
    menu_settings = uimenu(obj.window,'Label','Settings');        
    handles.menu_settings_microns_per_pixel = uimenu(menu_settings,'Label',['Microns per pixel = ' num2str(microns_per_pixel)]);                
    
    handles.menu_settings_problem = uimenu(menu_settings,'Label',['Application = ' problem],'Separator','on');                
    handles.menu_settings_problem_FungusDependentGranuleRelease = uimenu(handles.menu_settings_problem,'Label','Fungus Dependent Granule Release','Separator','on');
    handles.menu_settings_problem_CIDR = uimenu(handles.menu_settings_problem,'Label','CIDR','Separator','on');
    handles.menu_settings_problem_TTO = uimenu(handles.menu_settings_problem,'Label','TTO','Separator','on');
    handles.menu_settings_problem_PR = uimenu(handles.menu_settings_problem,'Label','PR','Separator','on');
    handles.menu_settings_problem_HL1 = uimenu(handles.menu_settings_problem,'Label','HL1','Separator','on');
    handles.menu_settings_problem_NucCyt = uimenu(handles.menu_settings_problem,'Label','NucCyt','Separator','on');
    handles.menu_settings_problem_MPHG = uimenu(handles.menu_settings_problem,'Label','MPHG','Separator','on');
    handles.menu_settings_problem_Sparks = uimenu(handles.menu_settings_problem,'Label','Sparks','Separator','on');    
    handles.menu_settings_problem_Experimental = uimenu(handles.menu_settings_problem,'Label','Experimental','Separator','on');
    handles.menu_settings_problem_per_image_TCSPC_FLIM = uimenu(handles.menu_settings_problem,'Label','per_image_TCSPC_FLIM','Separator','on');
    handles.menu_settings_problem_per_image_TCSPC_FLIM_PHASOR = uimenu(handles.menu_settings_problem,'Label','per_image_TCSPC_FLIM_PHASOR','Separator','on');
    handles.menu_settings_problem_t_dependent_Nuclei_ratio_FRET = uimenu(handles.menu_settings_problem,'Label','t_dependent_Nuclei_ratio_FRET','Separator','on');
    handles.menu_settings_problem_Image_Tiling = uimenu(handles.menu_settings_problem,'Label','Image_Tiling','Separator','on');
    handles.menu_settings_prblm_AI_Powered_2D_SMLM_Reconstruction = uimenu(handles.menu_settings_problem,'Label','AI_Powered_2D_SMLM_Reconstruction','Separator','on');
    handles.menu_settings_problem_OPT_ZFish_Embryo = uimenu(handles.menu_settings_problem,'Label','OPT_ZFish_Embryo','Separator','on');
    handles.menu_settings_problem_SIFNE = uimenu(handles.menu_settings_problem,'Label','SIFNE','Separator','on');
    
    handles.menu_settings_problem_dependent = uimenu(menu_settings,'Label','Application-Specific','Separator','on');                    
    handles.menu_settings_general = uimenu(menu_settings,'Label','General','Separator','on');                        
    handles.menu_settings_save = uimenu(menu_settings,'Label','Save Settings as...','Separator','on');
    handles.menu_settings_load = uimenu(menu_settings,'Label','Load Settings');
            
    %================================= visualization
    
    menu_visualization = uimenu(obj.window,'Label','Icy');
    menu_visualization_Icy_setup = uimenu(menu_visualization,'Label','Setup');        
    handles.menu_visualization_setup_Icy_directory = uimenu(menu_visualization_Icy_setup,'Label','Set Icy Directory');    
    handles.menu_visualization_start_Icy = uimenu(menu_visualization_Icy_setup,'Label','Start Icy');
    %
    handles.menu_visualization_send_original_to_Icy = uimenu(menu_visualization,'Label','Send Original','Separator','on');
    
    %================================= help   
    
%     menu_help = uimenu(obj.window,'Label','Help');
%     handles.menu_help_about = uimenu(menu_help,'Label','About...');
%     handles.menu_help_tracker = uimenu(menu_help,'Label','Open Issue Tracker...');
%     handles.menu_help_bugs = uimenu(menu_help,'Label','File Bug Report...');
        
end


