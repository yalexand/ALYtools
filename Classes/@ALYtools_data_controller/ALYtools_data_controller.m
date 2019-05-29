
classdef ALYtools_data_controller < handle 
    
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
   
    properties(Constant)
        data_settings_filename = 'ALYtools_data_settings.xml';
    end
    
    properties(SetObservable = true)
            
        microns_per_pixel = 1;
        downsampling = 1;
        %
        imgdata = [];
        
        problem = 'Fungus Dependent Granule Release';
        
        scene_axes = [];
        
        OPT_data_controller = [];
        NumWorkers = 2;
                                                                       
    end                    
    
    properties(Transient)
        
        DefaultDirectory = ['C:' filesep];
        RootDirectory = ['C:' filesep];
        
        IcyDirectory = [];
                
        current_filename;
        
        omero_Image_IDs = [];
                
        send_analysis_output_to_Icy = true;
        save_analysis_output_as_OMEtiff = true;
        save_analysis_output_as_xls = true;
        %
        send_original_to_Icy_on_show = false;
        ALYtools_always_on_top = false;

        % multiple image loading - data and filenames
        M_imgdata = [];
        M_filenames = [];
        M_sgm = []; % segmentation masks
        %
        
        
        %%%%%%%%%%%%%%% PROBLEM SPECIFIC - 'Fungus Dependent Granule Release'
        
        % NB - SIZES IN MICRONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% granules
        g_scale = 10;
        g_threshold = 5; 
        g_smoothing = 2;
        g_min_area = 9;
        g_rel_bg_scale = 3; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fungus
        f_scale = 40;
        f_threshold = 0.7; 
        f_smoothing = 5;
        f_min_area = 400;
        f_rel_bg_scale = 3; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells
        c_scale = 80;
        c_threshold = 3; 
        c_smoothing = 10;
        c_min_area = 500;
        c_rel_bg_scale = 4; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells        
        
        GRANULE_TO_CELL_DISTANCE_THRESHOLD = 20;
        %%%%%%%%%%%%%%% PROBLEM SPECIFIC - 'Fungus Dependent Granule Release'
        
        saved_segmentation_setups;
        
        % CIDR - "collagen induced DDR1 redistribution"
        CIDR_ref_channel = 1; % supposed to be overridden
        % segmentation
        CIDR_ripple_scale = 1;
        CIDR_ripple_threshold = 0.05;
        CIDR_rough_scale = 6;
        CIDR_rough_threshold = 0.3;        
        % analysis
        CIDR_smoothing_scale = 16;
        CIDR_Nmax = 500;
                
        % TTO - "T-tubule organization"
        TTO_ref_channel = 1; % supposed to be overridden
        % segmentation
        TTO_threshold = 0.1;
        TTO_smoothing_scale = 120;
        % analysis
        TTO_Lmin = 180;
        TTO_Lmax = 300;
        TTO_Nmax = 40000;
        TTO_fL = 100;
        TTO_fH = 6;
        TTO_df = 2;

        % Peak/Ridge
        PR_ref_channel = 1; % supposed to be overridden
        PR_coloc_channel = 2; % supposed to be overridden
        % segmentation       
        PR_K = 2.5;
        PR_S1 = 7;
        PR_S2 = 20;
        PR_a = 0.3;
        PR_t = 0.3;
        PR_mode = 'Peak';
        PR_min_size = 4;
        
        % HL1 - segmentation
        HL1_sgm_minimal_radius = 4;
        HL1_sgm_stdT_threshold = 5; % gray levels
        HL1_sgm_I_threshold = 0.7; % multiplier for golbal sum.img "std"        
        HL1_sgm_mode = 'Intensity std only'; % default
        HL1_sgm_modes  = {'Intensity std only' 'Time std only' 'Intensity AND Time std' 'Intensity OR Time std'};
        HL1_sgm_fill_holes = true;
        
        % HL1
        HL1_binning_radius = 6;
        HL1_ref_channel = 2;
        HL1_TD_avreraging_size = 40; % high pass
        HL1_TD_smoothing_size = 0; % low pass
        HL1_df = 2;
        HL1_invert_input = false;
        %
        %minpkdst = 11; % min or max? when T = 14
        %minpkdst = 23; % min or max? when T = 30
        %minpkdst = 290; % min or max? when T ~ 300                
        HL1_calculate_phasemap = false;        
        HL1_MINPEAKDISTANCE = 11;
        HL1_MINPEAKHEIGHT_std_factor = 0.5;
        HL1_min_std_T = 5;
        HL1_dynamic_amplitude_threshold = 6;
        %
        HL1_delineate_wave_front = false;
        HL1_create_LF_movie = true;        
        HL1_create_frame_indexed_excited_regions_movie = true;
        H1_isochrones_map_step = 2; %frames        
        %        
        % NC - problem specific
        NC_chNuc = 1;
        NC_chCell = 2;
        NC_bckg_subtraction_proportion = 1.; % 100% of background intensity will be subtracted        
        NC_bckg_dilation_size = Inf;
        %
        %NC - segmentation
        NC_cell_smoothing_radius = 5;
        NC_cell_rg_std_factor_int = 0.5;
        NC_cell_rg_std_factor_out = 0.5;
        NC_cell_overpeak_ratio = 1; % 0<param<1, used: thershold = graypeak + param*(graypeak - image_min);
        %
        NC_nuc_scale = 200;
        NC_nuc_threshold = 0.05; 
        NC_nuc_smoothing = 5;
        NC_nuc_min_area = 100;
        NC_nuc_rel_bg_scale = 2.5; 
        % 
        NC_nuc_breacking_distmap_smoothing_scale = 6;                        
        %        
        % experimental
        expar = [];
              
            MPHG_chNuc = 1;
            MPHG_chCell = 4;
            MPHG_chGran = 2;
                        
            % adjusted for 0.44u/pix             
            % segment nukes            
            MPHG_nuc_scale = 27*0.43945;
            MPHG_nuc_threshold = 1.;
            MPHG_nuc_smoothing_scale = 4*0.43945;
            MPHG_nuc_min_area = 400*0.43945*0.43945;
            MPHG_nuc_rel_bg_scale = 2.5; 
            
            % segment cells and clumps            
            MPHG_cell_smoothing_radius = 3*0.43945;
            MPHG_cell_rg_std_factor_int = 0.5;
            MPHG_cell_rg_std_factor_out = 0.5;
            MPHG_cell_overpeak_ratio = 1.5;                
             
            % segment granules
            MPHG_minimal_gran_size = 2*0.43945*0.43945;
            MPHG_gran_overpeak_ratio = 1;        
            
            per_image_TCSPC_FLIM_t = [];
            per_image_TCSPC_FLIM_irf = [];
            per_image_TCSPC_FLIM_irf_shift = 0; %[ps]
            per_image_TCSPC_FLIM_rep_rate = 20; % [MHz]
            per_image_TCSPC_FLIM_irf_filename = [];
            global_fit_models = {'1exp', '2exp', 'sk2', 'sk2 only', 'double sk2'};
            %
            per_image_TCSPC_FLIM_fit_model = '2exp';
            per_image_TCSPC_FLIM_Tmin = 0; % [ps]
            per_image_TCSPC_FLIM_Tmax = 0; % [ps]
            per_image_TCSPC_FLIM_background_value = 0; % intensity
            per_image_TCSPC_FLIM_saturation_value = Inf;            
            
            per_image_TCSPC_FLIM_tvb = [];
            per_image_TCSPC_FLIM_tvb_filename = [];
            
            per_image_TCSPC_FLIM_fixed_tauD = 0;
            
            per_image_TCSPC_FLIM_nonimaging = false; % not saved
            per_image_TCSPC_FLIM_conv_irf_pp_69_70 = true; % saved
            per_image_TCSPC_FLIM_irf_background = 0;
            per_image_TCSPC_FLIM_tvb_scaling = 0; 
            per_image_TCSPC_FLIM_weights_resampling_factor = 2;             
            
            per_image_TCSPC_FLIM_Ref_lifetime = 0;
            per_image_TCSPC_FLIM_averaging_sigma = 0;
            %
            external_segmentations = [];
            
            %
            FijiScriptsDirectory = 'C:\Fiji.app\scripts';            
            %
            t_dependent_Nuclei_ratio_FRET_TIMESTEP = 7/60;
            TrackMate_RADIUS = 5;
            TrackMate_TARGET_CHANNEL = 1;
            TrackMate_THRESHOLD = 0;
            TrackMate_DO_MEDIAN_FILTERING = true;
            TrackMate_QUALITY = 1.0;
            TrackMate_ALLOW_TRACK_SPLITTING = false;
            TrackMate_ALLOW_TRACK_MERGING = false;
            TrackMate_TRACK_DISPLACEMENT = 2.0;    
            TrackMate_LINKING_MAX_DISTANCE = 10; % pixels
            TrackMate_GAP_CLOSING_MAX_DISTANCE = 10; % pixels
            TrackMate_MAX_FRAME_GAP = 1; % frames
            %
            % PARAMETERS USED TO CALCULATE FRET MOLAR FRACTION
            % excitation exctinction coefficient
            t_dependent_Nuclei_ratio_FRET_eps_A = 3000;       %   %% 83400;       % EYFP - this is on peak
            t_dependent_Nuclei_ratio_FRET_eps_D = 32500;       % ECFP [M-1 cm-1]
            % optical system transmission
            t_dependent_Nuclei_ratio_FRET_t_A  = 1;
            t_dependent_Nuclei_ratio_FRET_t_D  = 1;
            % quantum yield at room temperature
            t_dependent_Nuclei_ratio_FRET_Q_A = 0.61; % [ns] EYFP  
            t_dependent_Nuclei_ratio_FRET_Q_D = 0.40; % [ns] ECFP  
            % lifetimes
            t_dependent_Nuclei_ratio_FRET_tau_A = 3.100; % [ns] EYFP  
            t_dependent_Nuclei_ratio_FRET_tau_D = 3.000; % [ns] ECFP 
            % FRETting donor lifetime
            t_dependent_Nuclei_ratio_FRET_tau_FRET = 0.700; % [ns] - HIGH FRET
            % spectral leakage coefficients
            t_dependent_Nuclei_ratio_FRET_K_DA = 0;
            t_dependent_Nuclei_ratio_FRET_K_AD = 0;
            % ratio of excitation intensities in the donor and acceptor absorption bands
            % phi = IexD/IexA
            t_dependent_Nuclei_ratio_FRET_phi = 1;
            % fraction of functional donor and acceptor
            t_dependent_Nuclei_ratio_FRET_b_A = 1; 
            t_dependent_Nuclei_ratio_FRET_b_D = 1;
            %
            t_dependent_Nuclei_ratio_FRET_donor_channel = 1;
            %
            ImageTiling_Ncols = 0;
            ImageTiling_Nrows = 0;
            ImageTiling_Ovlp_X = 0.5;
            ImageTiling_Ovlp_Y = 0.5;            
            ImageTiling_QT = 0;
            ImageTiling_mode = 'bleached_fluor'; %'brightfield'
            %
            AI_Powered_2D_SMLM_Reconstruction_network = []; % path to file
            AI_Powered_2D_SMLM_Reconstruction_upscale_factor = 5;
            AI_Powered_2D_SMLM_Reconstruction_vicinity_half_width = 5;
            AI_Powered_2D_SMLM_Reconstruction_pixel_size = 107; % nm/pixel
            AI_Powered_2D_SMLM_Reconstruction_extraction_scale = 2; %non-super res pixels            
            AI_Powered_2D_SMLM_Reconstruction_extraction_scale_ratio = 3.5; % ditto
            AI_Powered_2D_SMLM_Reconstruction_extraction_threshold = 6; % std factor?
            AI_Powered_2D_SMLM_Reconstruction_extraction_method = 'Linear_TopHat';
            AI_Powered_2D_SMLM_Reconstruction_max_distance_to_spurious_pixl = 50; % nm
            AI_Powered_2D_SMLM_Reconstruction_time_dependent_block_size = 200; % frames            
            AI_Powered_2D_SMLM_Reconstruction_image_formation_method = 'ASH';
            AI_Powered_2D_SMLM_Reconstruction_image_formation_scale = 3;
            AI_Powered_2D_SMLM_Reconstruction_NA = 1;
            AI_Powered_2D_SMLM_Reconstruction_extraction_methods = {'Linear_TopHat','DOG','Primitive_Linear_Tophat'};
            AI_Powered_2D_SMLM_Reconstruction_image_formation_methods = {'ASH','Smoothed Histogram'};
            AI_Powered_2D_SMLM_Reconstruction_wavelength = 509; % nm
            AI_Powered_2D_SMLM_Reconstruction_min_sigma = 0; % pixel
            AI_Powered_2D_SMLM_Reconstruction_max_sigma = 7; % pixel

    end    
        
    properties(Transient,Hidden)
        % Properties that won't be saved to a data_settings_file etc.
        
        menu_controller;
        
        h_Icy_segmentation_adjustment = [];
        
        run_headless = false;
        
    end
    
    events
        
    end
            
    methods
        
        function obj = ALYtools_data_controller(RUN_HEADLESS,varargin)            
            %   
            obj.run_headless = RUN_HEADLESS;
            %
            if ~isempty(varargin{1})
                handles = args2struct(varargin);                
                assign_handles(obj,handles);                     
                %                    
                scene_panel = handles.scene_panel;              
                obj.scene_axes = axes( 'Parent', scene_panel,'ActivePositionProperty', 'OuterPosition' );
                axis image; set(obj.scene_axes,'XTick',[],'YTick',[]);            
            end
                        
            try 
                obj.load_settings([pwd filesep obj.data_settings_filename]);
            catch
            end
            
            %obj.OPT_data_controller = ic_OPTtools_data_controller(true,[]); %headless
                        
            obj.RootDirectory = pwd;
            
            if (ispc || ismac) && isempty(obj.IcyDirectory) && ~isdeployed && ~obj.run_headless
                hw = waitbar(0,'looking for Icy directory..');
                waitbar(0.1,hw);                
                if ispc
                       prevdir = pwd;
                       cd('c:\');
                       [~,b] = dos('dir /s /b icy.exe');
                       if ~strcmp(b,'File Not Found')
                            filenames = textscan(b,'%s','delimiter',char(10));
                            s = char(filenames{1});
                            s = s(1,:);
                            s = strsplit(s,'icy.exe');
                            obj.IcyDirectory = s{1};
                       end
                       cd(prevdir);
                elseif ismac
                    % to do
                else
                    % to do
                end                
                delete(hw); drawnow;                               
            end
            %
            % if there are more than 4 cores, one can try using parfor
            % par_pool = gcp;
            % obj.NumWorkers = par_pool.NumWorkers;
                                                                                                                                                   
        end
                
%-------------------------------------------------------------------------%
        function delete(obj)
            
            if ~obj.run_headless
                ButtonName = questdlg('Do you want to save current settings?', ...
                             'Settings saving on exit', ...
                             'Yes', 'No - quit without settings saving', 'Yes');
                if strcmp(ButtonName,'Yes')                                                 
                    try
                        obj.save_settings([obj.RootDirectory filesep obj.data_settings_filename]);
                    catch
                        disp('Error while trying to save settings: not saved');
                    end
                    close all; % radical
                end
            end
            
        end       
%-------------------------------------------------------------------------%                
         function clear_all(obj,~)
             obj.imgdata = [];
             %
             % to do
             %
         end
%-------------------------------------------------------------------------%        
        function load_multiple(obj,filenames,pathname,verbose,~)
                                    
            obj.M_imgdata = [];
            obj.M_filenames = [];
            
            if ischar(filenames) %single file
                if ~obj.open_image([pathname filenames]), return, end;
                obj.M_imgdata{1} = obj.imgdata;
                obj.show_data(obj.send_original_to_Icy_on_show);
                obj.M_filenames{1} = filenames;
                obj.current_filename = filenames; % 24.04.2019 : likely low risk
            else
                obj.imgdata = [];
                hw = [];
                if verbose, hw = waitbar(0,'Loading files, please wait'); end;
                for k=1:numel(filenames)       
                    if ~obj.open_image([pathname filesep char(filenames(k))]), break, end;
                    obj.M_imgdata{k} = obj.imgdata;
                    obj.M_filenames{k} = filenames(k);
                    if ~obj.run_headless, obj.show_data(obj.send_original_to_Icy_on_show); end;
                    if ~isempty(hw), waitbar(k/numel(filenames),hw); drawnow, end; 
                end
                if ~isempty(hw), delete(hw), drawnow; end;
                % setting obj.current_filename - presuming current folder's name is relevant
                strs = strsplit(pathname,filesep);
                obj.current_filename = char(strs(numel(strs)));
                if isempty(obj.current_filename)
                    obj.current_filename = char(strs(numel(strs)-1));
                end
            end
            
            %obj.imgdata = [];
            obj.DefaultDirectory = pathname;            
            
        end      
%-------------------------------------------------------------------------%        
        function load(obj,filename,pathname,verbose,~)
            obj.M_imgdata = [];
            obj.open_image([pathname filesep filename]);                 
            obj.current_filename = filename; 
            if ~isempty(obj.scene_axes)
                obj.show_data(obj.send_original_to_Icy_on_show);
            end
            obj.DefaultDirectory = pathname;            
        end      
        
 %-------------------------------------------------------------------------%           
        function ret = open_image(obj,full_filename,~)
            
            ret = true;
            
            full_filename = strrep(full_filename,[filesep filesep],filesep);
            
            obj.imgdata = [];
            obj.OPT_data_controller.delays = [];
            obj.OPT_data_controller.proj = [];
            obj.OPT_data_controller.volm = [];
            
            % check if that is maybe OPT projections data
            try
                infostring = obj.OPT_data_controller.Set_Src_Single(full_filename,false);
                if ~isempty(infostring)
                    obj.imgdata = obj.OPT_data_controller.proj(:,:,1); % to visualize
                    return;                    
                end
            catch
            end
                        
            extension = full_filename(length(full_filename)-2:length(full_filename));            
            if strcmp(extension,'txt')
                try
                    data = importdata(full_filename);
                    if isnumeric(data)
                        data_v = data(:,2);
                    else
                        nrows = size(data.data,1);
                        if nrows - 256 < 20 % hack
                            L = 256;
                        end
                         offset = nrows - L + 1 ;
                         data_v = data.data(offset:nrows,2);
                    end
                    I = zeros(1,1,1,1,numel(data_v));
                    I(1,1,1,1,:) = data_v;
                    obj.per_image_TCSPC_FLIM_nonimaging = true;
                catch
                    ret = false;
                    errordlg('*.txt file seems to be corrupt - can not load');
                    return;
                end               
            else                          
                try
                    [uppixX,uppixZ,I] = bfopen_v(full_filename);
                    if ~isempty(uppixX)
                        %obj.microns_per_pixel = uppixX;
                    else
                        disp('No spatial resolution information found, set obj.microns_per_pixel to settings value');
                    end;
                    obj.per_image_TCSPC_FLIM_nonimaging = false;
                catch
                    errordlg('error on loading - can not continue');
                    ret = false;
                    return; % mmm
                end
            end
            
            [sX,sY,sZ,sC,sT] = size(I);
            
            switch obj.problem 
                 case 'Fungus Dependent Granule Release'
                     if sC~=3
                         errordlg('3 channel image is expected, can not load'),
                         return,
                     end
                     obj.imgdata = I;
                     
                 case 'PR' 
                     if sZ~=1 % collapse Z - don't do 3d analysis
                         obj.imgdata = reshape(squeeze(sum(I,3)),[sX,sY,1,sC,1]);
                     else
                         obj.imgdata = I;
                         if obj.PR_ref_channel > sC
                             obj.PR_ref_channel = 1;
                         end
                        if obj.PR_coloc_channel > sC
                             obj.PR_coloc_channel = 1;
                         end                         
                     end
                     
                case 'CIDR' 
                     obj.imgdata = I; 
                         if obj.CIDR_ref_channel > sC
                             obj.CIDR_ref_channel = 1;
                         end 
                         
                case 'TTO'                    
                     obj.imgdata = I;
                         if obj.TTO_ref_channel > sC
                             obj.TTO_ref_channel = 1;
                         end 
                         
                case 'HL1'                    
                     obj.imgdata = I;                        
                         if obj.HL1_ref_channel > sC
                             obj.HL1_ref_channel = 1;
                         end                                               
                         
                case { 'Experimental', 'per_image_TCSPC_FLIM','per_image_TCSPC_FLIM_PHASOR' 't_dependent_Nuclei_ratio_FRET', 'Image_Tiling', 'AI_Powered_2D_SMLM_Reconstruction'}
                    
                     % if sT<6, don't believe it is T, reinterpret as channels
                     if sT < 6 && sC ==1
                        I = reshape(I,[sX,sY,sZ,sT,sC]);
                     end
                     
                     % if sC>12, don't believe it is C, reinterpret as Z
                     if sC > 12 && sZ ==1
                        I = reshape(I,[sX,sY,sC,sZ,sT]);
                     end
                     
                     %                   
                     obj.imgdata = I;                        
                                          
                case 'MPHG'
                    
                     % if sT<6, don't believe it is T, reinterpret as channels
                     if sT < 6 && sC ==1
                        I = reshape(I,[sX,sY,sZ,sT,sC]);
                     end
                     %                   
                     obj.imgdata = I;                        

                case 'Sparks'   
                    
                     % if sT<6, don't believe it is T, reinterpret as channels
                     if sT < 6 && sC ==1
                        I = reshape(I,[sX,sY,sZ,sT,sC]);
                     end
                     %                   
                     obj.imgdata = I;                        
                                          
                case 'NucCyt'   
                    
                     % if sT<6, don't believe it is T, reinterpret as channels
                     if sT < 6 && sC ==1
                        I = reshape(I,[sX,sY,sZ,sT,sC]);
                     end
                     %                   
                     obj.imgdata = I;                                                                
            end  
            
            if ~obj.ALYtools_always_on_top && ~isempty(obj.menu_controller)
                WinOnTop(obj.menu_controller.window,false); % undo "always on top" after 1st load - too annoying
            end
        end       
%-------------------------------------------------------------------------%        

        function pixsize = u2pix(obj,usize,~)
            pixsize = max(fix(usize/obj.microns_per_pixel/obj.downsampling),1);
        end
        function pix2size = u22pix2(obj,u2size,~)
            pix2size = max(fix(u2size/obj.microns_per_pixel/obj.downsampling/obj.microns_per_pixel/obj.downsampling),1);
        end
        
%-------------------------------------------------------------------------%        
        function show_data(obj,send_original_to_Icy,~)

            [~,~,~,sC,sT] = size(obj.imgdata);
            
            switch obj.problem
                
                case 'Fungus Dependent Granule Release'
                    I = double(obj.imgdata);
                    if 3~=sC, obj.imgdata = []; cla(obj.scene_axes); errordlg('3 channel image is expected'); return; end
                    image(cat(3,uint8(map(I(:,:,1,1,1),0,255)),uint8(map(I(:,:,1,2,1),0,255)),uint8(map(I(:,:,1,3,1),0,255))),'Parent',obj.scene_axes); 
                    
                case {'CIDR' 'TTO' 'PR' 'HL1' 'Experimental' 'NucCyt' 'MPHG' 'Sparks','Image_Tiling','AI_Powered_2D_SMLM_Reconstruction'}
                    if sT>10
                        I = double(obj.imgdata(:,:,:,:,1:10));
                    else
                        I = double(obj.imgdata);                        
                    end
                    if 3 == sC || 4 == sC
                        image(cat(3,uint8(map(I(:,:,1,1,1),0,255)),uint8(map(I(:,:,1,2,1),0,255)),uint8(map(I(:,:,1,3,1),0,255))),'Parent',obj.scene_axes);
                    elseif (1 == sC || 2 == sC) && 1==sT
                        image(uint8(map(squeeze(I(:,:,1,1,1)),0,255)),'Parent',obj.scene_axes);
                    elseif (1 == sC) && sT>1
                        u = squeeze(sum(I,5));
                        image(uint8(map(u,0,255)),'Parent',obj.scene_axes);
                    end  
                case {'per_image_TCSPC_FLIM','per_image_TCSPC_FLIM_PHASOR'}
                    image(uint8(map(squeeze(sum(I,5)),0,255)),'Parent',obj.scene_axes);
                case 't_dependent_Nuclei_ratio_FRET'                    
                    image(uint8(map(double(squeeze(obj.imgdata(:,:,1,1,1))),0,255)),'Parent',obj.scene_axes);                                        
            end
            %
            daspect(obj.scene_axes,[1 1 1]);
            set( obj.scene_axes, 'xticklabel', [], 'yticklabel', [] );
            xlabel(obj.scene_axes,obj.current_filename);
            
            if send_original_to_Icy && ~obj.per_image_TCSPC_FLIM_nonimaging 
                obj.send_original_to_Icy;
            end
        end
%-------------------------------------------------------------------------%                                        
        function Set_Fungus_Dependent_Granule_Release_segmentation_default(obj,~,~)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% granules
            obj.g_scale = 10;
            obj.g_threshold = 5; 
            obj.g_smoothing = 2;
            obj.g_min_area = 9;
            obj.g_rel_bg_scale = 3; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fungus
            obj.f_scale = 40;
            obj.f_threshold = 0.7; 
            obj.f_smoothing = 5;
            obj.f_min_area = 400;
            obj.f_rel_bg_scale = 3; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells
            obj.c_scale = 80;
            obj.c_threshold = 2; 
            obj.c_smoothing = 10;
            obj.c_min_area = 100;
            obj.c_rel_bg_scale = 3; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells                    
        end
%-------------------------------------------------------------------------% 
        function save_Fungus_Dependent_Granule_Release_segmentation_setups(obj,~,~)
            obj.saved_segmentation_setups = [];
            obj.saved_segmentation_setups.g_scale = obj.g_scale;
            obj.saved_segmentation_setups.g_threshold = obj.g_threshold;
            obj.saved_segmentation_setups.g_smoothing = obj.g_smoothing;
            obj.saved_segmentation_setups.g_min_area = obj.g_min_area;
            obj.saved_segmentation_setups.g_rel_bg_scale = obj.g_rel_bg_scale;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fungus
            obj.saved_segmentation_setups.f_scale = obj.f_scale;
            obj.saved_segmentation_setups.f_threshold = obj.f_threshold;
            obj.saved_segmentation_setups.f_smoothing = obj.f_smoothing;
            obj.saved_segmentation_setups.f_min_area = obj.f_min_area;
            obj.saved_segmentation_setups.f_rel_bg_scale = obj.f_rel_bg_scale;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells
            obj.saved_segmentation_setups.c_scale = obj.c_scale;
            obj.saved_segmentation_setups.c_threshold = obj.c_threshold;
            obj.saved_segmentation_setups.c_smoothing = obj.c_smoothing;
            obj.saved_segmentation_setups.c_min_area = obj.c_min_area;
            obj.saved_segmentation_setups.c_rel_bg_scale = obj.c_rel_bg_scale;
            
            % TTO
            obj.saved_segmentation_setups.TTO_threshold = obj.TTO_threshold;
            obj.saved_segmentation_setups.TTO_smoothing_scale = obj.TTO_smoothing_scale;
        end
%-------------------------------------------------------------------------% 
        function restore_Fungus_Dependent_Granule_Release_segmentation_setups(obj,~,~)
            obj.g_scale = obj.saved_segmentation_setups.g_scale;
            obj.g_threshold = obj.saved_segmentation_setups.g_threshold;
            obj.g_smoothing = obj.saved_segmentation_setups.g_smoothing;
            obj.g_min_area = obj.saved_segmentation_setups.g_min_area;
            obj.g_rel_bg_scale = obj.saved_segmentation_setups.g_rel_bg_scale;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fungus
            obj.f_scale = obj.saved_segmentation_setups.f_scale;
            obj.f_threshold = obj.saved_segmentation_setups.f_threshold;
            obj.f_smoothing = obj.saved_segmentation_setups.f_smoothing;
            obj.f_min_area = obj.saved_segmentation_setups.f_min_area;
            obj.f_rel_bg_scale = obj.saved_segmentation_setups.f_rel_bg_scale;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells
            obj.c_scale = obj.saved_segmentation_setups.c_scale;
            obj.c_threshold = obj.saved_segmentation_setups.c_threshold;
            obj.c_smoothing = obj.saved_segmentation_setups.c_smoothing;
            obj.c_min_area = obj.saved_segmentation_setups.c_min_area;
            obj.c_rel_bg_scale = obj.saved_segmentation_setups.c_rel_bg_scale;            
        end
%-------------------------------------------------------------------------% 
        function restore_TTO_segmentation_setups(obj,~,~)
            obj.TTO_threshold = obj.saved_segmentation_setups.TTO_threshold;
            obj.TTO_smoothing_scale = obj.saved_segmentation_setups.TTO_smoothing_scale;
        end
%-------------------------------------------------------------------------%         
        function Set_TTO_segmentation_default(obj,~,~)
                    obj.TTO_threshold = 0.1;
                    obj.TTO_smoothing_scale = 120;
        end
%-------------------------------------------------------------------------%        
        function Set_CIDR_segmentation_default(obj,~,~)
            % segmentation
            obj.CIDR_ripple_scale = 1;
            obj.CIDR_ripple_threshold = 0.05;
            obj.CIDR_rough_scale = 6;
            obj.CIDR_rough_threshold = 0.3;                    
        end        
%-------------------------------------------------------------------------%                
        function Set_PR_segmentation_default(obj,~,~)
            obj.PR_K = 2.5;
            obj.PR_S1 = 7;
            obj.PR_S2 = 20;
            obj.PR_a = 0.3;
            obj.PR_t = 0.3;
            obj.PR_mode = 'Peak';
            obj.PR_min_size = 4;
        end                
%-------------------------------------------------------------------------%                
        function Set_HL1_segmentation_default(obj,~,~)
            obj.HL1_sgm_minimal_radius = 4;
            obj.HL1_sgm_stdT_threshold = 5; % gray levels
            obj.HL1_sgm_I_threshold = 0.7; % multiplier for golbal sum.img "std"        
            obj.HL1_sgm_mode = 'Intensity std only'; % default
            obj.HL1_sgm_modes  = {'Intensity std only' 'Time std only' 'Intensity AND Time std' 'Intensity OR Time std'};
            obj.HL1_sgm_fill_holes = true;            
        end                        
%-------------------------------------------------------------------------%
        function adjust_segmentation(obj,~,~)

            if isempty(obj.imgdata), errordlg('no data!'), return, end;
            
            switch obj.problem

                case 'Fungus Dependent Granule Release'                    
                     obj.save_Fungus_Dependent_Granule_Release_segmentation_setups;
                     Granule_Fungus_Cells_Segmentation_settings(obj);
                    
                case 'CIDR'                    
                     CIDR_Segmentation_settings(obj);
                    
                case 'TTO'                    
                     TTO_Segmentation_settings(obj);                                       
                     
                case 'PR'                    
                     PR_Segmentation_settings(obj);                                         
                
                case 'HL1'
                     HL1_Segmentation_settings(obj);

                case 'NucCyt'
                     NucCyt_Segmentation_settings(obj);                                   

                case 'MPHG'
                     MPHG_Segmentation_settings(obj);                     

                case 'Sparks'
                     %Sparks_Segmentation_settings(obj);                     
                                          
                case 'Experimental'
                     % Experimental_Segmentation_settings(obj);
                     obj.do_Experimental_Segmentation(true);                                     
                     
                case {'per_image_TCSPC_FLIM','per_image_TCSPC_FLIM_PHASOR'}
                    % parasytise
                    if isempty(obj.M_imgdata), return, end;
                    k = numel(obj.M_imgdata); % last
                    obj.PR_ref_channel = 1;
                    obj.imgdata = squeeze(sum(obj.M_imgdata{k},5));
                    PR_Segmentation_settings(obj);     
                    
                case 't_dependent_Nuclei_ratio_FRET'
                    % t_dependent_Nuclei_ratio_FRET_Segmentation_settings(obj);
                    
                case 'ImageTiling'
                    obj.do_ImageTiling_Segmentation(true);
                    
                case 'AI_Powered_2D_SMLM_Reconstruction'
                    obj.do_AI_Powered_2D_SMLM_Reconstruction_Segmentation(true);                                                                                
            end

        end    
%-------------------------------------------------------------------------%                                        
        function ret = do_Fungus_Dependent_Granule_Release_segmentation(obj,send_to_Icy,~)
                         
            granules    = squeeze(double(obj.imgdata(:,:,1,1,1)));
            fungus      = squeeze(double(obj.imgdata(:,:,1,2,1)));
            cells       = squeeze(double(obj.imgdata(:,:,1,3,1)));            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% granules
            segmented_granules = nth_segmentation(granules, ...
                obj.u2pix(obj.g_scale), ...
                obj.g_rel_bg_scale, ...
                obj.g_threshold, ...
                obj.u2pix(obj.g_smoothing), ...
                obj.u22pix2(obj.g_min_area));
            segmented_granules(segmented_granules~=0)=1;
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fungus
            segmented_fungus = nth_segmentation(fungus, ...
                obj.u2pix(obj.f_scale), ...
                obj.f_rel_bg_scale, ...
                obj.f_threshold, ...
                obj.u2pix(obj.f_smoothing), ...
                obj.u22pix2(obj.f_min_area));            
            segmented_fungus(segmented_fungus~=0)=1;
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells
            segmented_cells = nth_segmentation(cells, ...
                obj.u2pix(obj.c_scale), ...
                obj.c_rel_bg_scale, ...
                obj.c_threshold, ...
                obj.u2pix(obj.c_smoothing), ...
                obj.u22pix2(obj.c_min_area));                        
            segmented_cells(segmented_cells~=0)=1;
                        
          % (allegedly) an improved version of watershed on cells...
            z2 = segmented_cells;
            D = bwdist(~z2); %distance map
            D = medfilt2(D,[5 5]); % a bit voluntaristic..
                D = -D;
                D(~z2) = -Inf;
            	L = watershed(D);
            % remove background    
            stats = regionprops(L,'Area');    
            bckgind = find([stats.Area]==max([stats.Area]));
            L(L==bckgind) = 0;
            %remove small stuff
            idx = find([stats.Area] > obj.u22pix2(obj.c_min_area));
            L = ismember(L,idx);
            L = bwlabel(L);
            %
            idemp = bwmorph(L,'thicken',Inf);            
            segmented_cells(idemp==0)=0;
            segmented_cells = imclearborder(segmented_cells);                              
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ret = cat(3,segmented_granules,segmented_fungus,segmented_cells);
            
            if send_to_Icy                
                try
                    icyvol(:,:,1,1,1) = segmented_granules;
                    icyvol(:,:,2,1,1) = segmented_fungus;
                    icyvol(:,:,3,1,1) = segmented_cells;
                    if isempty(obj.h_Icy_segmentation_adjustment)
                        obj.h_Icy_segmentation_adjustment = icy_imshow(icyvol,'Segmentation: Fungus Dependent Granule Release');                    
                    else
                        icy_imshow(obj.h_Icy_segmentation_adjustment,icyvol,'Segmentation: Fungus Dependent Granule Release');                    
                    end
                catch
                    errodlg('problem with Icy, - might be not running');
                end                
            end
            
        end
%-------------------------------------------------------------------------%                                    
        function analyze_current(obj,~,~)        

            if isempty(obj.imgdata), errordlg('no data!'), return, end
            
            switch obj.problem
                
                case 'Fungus Dependent Granule Release' 
                    [datas, captions, table_names, fig] = obj.analyze_Fungus_Dependent_Granule_Release;
                case 'CIDR'                                      
                    [datas, captions, table_names, fig] = obj.analyze_CIDR;
                case 'TTO'                                      
                    [datas, captions, table_names, fig] = obj.analyze_TTO;
                case 'PR'                                      
                    [datas, captions, table_names, fig] = obj.analyze_PR;
                case 'HL1'                                      
                    [datas, captions, table_names, fig] = obj.analyze_HL1; 
                case 'Experimental'                    
                    [datas, captions, table_names, fig] = obj.analyze_Experimental; 
                case 'NucCyt'
                    [datas, captions, table_names, fig] = obj.analyze_NC; 
                case 'MPHG'
                    [datas, captions, table_names, fig] = obj.analyze_MPHG;                     
                case 'Sparks'
                    [datas, captions, table_names, fig] = obj.analyze_Sparks;                    
                case 'per_image_TCSPC_FLIM'
                    [datas, captions, table_names, fig] = obj.analyze_per_image_TCSPC_FLIM;
                case 'per_image_TCSPC_FLIM_PHASOR'
                    [datas, captions, table_names, fig] = obj.analyze_per_image_TCSPC_FLIM_PHASOR;
                case 't_dependent_Nuclei_ratio_FRET'
                    [datas, captions, table_names, fig] = obj.analyze_t_dependent_Nuclei_ratio_FRET;
                case 'Image_Tiling'
                    [datas, captions, table_names, fig] = obj.analyze_ImageTiling;
                case 'AI_Powered_2D_SMLM_Reconstruction'
                    [datas, captions, table_names, fig] = obj.analyze_AI_Powered_2D_SMLM_Reconstruction;
                    
            end % switch

            timestamp = datestr(now,'yyyy-mm-dd HH-MM-SS');
            dirname = [obj.RootDirectory filesep ['ALYtools Analysis ' timestamp]];
            %%%%
            if strcmp(obj.problem,'t_dependent_Nuclei_ratio_FRET') % special case where lot of data are saved
                mkdir(obj.RootDirectory,['ALYtools Analysis ' timestamp]);
                fig = obj.t_dependent_Nuclei_ratio_FRET_postprocess(fig,dirname);                
            end
            %%%%
            if (obj.save_analysis_output_as_xls && ~isempty(datas)) || obj.save_analysis_output_as_OMEtiff 
                mkdir(obj.RootDirectory,['ALYtools Analysis ' timestamp]);
            end                   
                    %  
                    if strcmp(obj.problem,'Fungus Dependent Granule Release') % several tables
                        if obj.save_analysis_output_as_xls && ~isempty(datas)
                            xlsname = [dirname filesep obj.current_filename '_analysis_data.xls']; 
                            for k=1:numel(datas)
                                if ispc
                                    xlswrite( xlsname,[captions{k}; datas{k}],char(table_names{k}) );
                                else
                                    xlwrite( xlsname,[captions{k}; datas{k}],char(table_names{k}) );
                                end
                            end                                        
                        end
                    elseif strcmp(obj.problem,'CIDR') || strcmp(obj.problem,'TTO') ... 
                            || strcmp(obj.problem,'PR') || strcmp(obj.problem,'NucCyt') ...
                            || strcmp(obj.problem,'MPHG') || strcmp(obj.problem,'Experimental') ...
                            || strcmp(obj.problem,'Sparks') || strcmp(obj.problem,'per_image_TCSPC_FLIM') ... 
                            || strcmp(obj.problem,'Image_Tiling') ...
                            || strcmp(obj.problem,'AI_Powered_2D_SMLM_Reconstruction')
                        if obj.save_analysis_output_as_xls && ~isempty(datas)
                            if ~isempty(obj.current_filename)
                                xlsname = [dirname filesep '_' obj.problem '_' obj.current_filename '_analysis_data.xls'];
                            else
                                xlsname = [dirname filesep obj.problem '_analysis_data.xls'];
                            end
                            if ispc
                                try
                                    xlswrite( xlsname,[captions; datas],char(table_names) );
                                catch
                                    disp('can not write output as xls, save as mat file instead');                                    
                                    fname = xlsname(1:length(xlsname)-4);
                                    save([fname '.mat'],'captions','datas');
                                    cell2csv([fname '.csv'],[captions; datas]);
                                end
                            else
                                try
                                    xlwrite( xlsname,[captions; datas],char(table_names) );
                                catch
                                    disp('can not write output as xls, save as mat file instead');
                                    fname = xlsname(1:length(xlsname)-4);
                                    save([fname '.mat'],'captions','datas');
                                    cell2csv([fname '.csv'],[captions; datas]);                                    
                                end
                            end
                        end
                    end
                    %
                    if obj.send_analysis_output_to_Icy && ~isempty(fig)
                        try
                            if strcmp(obj.problem,'HL1')
                                icy_imshow(fig.icyvol,[obj.current_filename ' - HL1 results']);                             
                                icy_imshow(fig.rotvol,[obj.current_filename ' - HL1 dynamics, T = ' num2str(fig.T)]);
                                if obj.HL1_create_LF_movie
                                    icy_imshow(fig.LFvol,[obj.current_filename ' - HL1 trend']);     
                                end
                                if ~isempty(fig.isochrones)
                                    icy_imshow(fig.isochrones,[obj.current_filename ' - HL1 isochrones']);     
                                end
                            elseif strcmp(obj.problem,'per_image_TCSPC_FLIM')
                                % fig is many images so throw them all
%                                 for k=1:numel(fig)
%                                     badge = [ char(obj.M_filenames{k}) ' at ' timestamp];
%                                     icy_imshow(fig{k},badge);
%                                 end                                
                            else
                                icy_imshow(fig,[obj.current_filename ' ' obj.problem ' analysis']);
                            end
                        catch
                            errordlg('error occurred when sending data to Icy - some visuals may be missing');
                        end;
                    end
                    %
                    if obj.save_analysis_output_as_OMEtiff  && ~isempty(fig) && isnumeric(fig)                   
                                if strcmp(obj.problem,'HL1')
                                    ometiffsavename = [dirname filesep obj.current_filename '_analysis_output.OME.tiff'];
                                    bfsave(fig.icyvol,ometiffsavename,'Compression','LZW','dimensionOrder','XYCTZ');
                                    %
                                    % another 2 interesting files
                                    ometiffsavename_rotor = [dirname filesep obj.current_filename '_analysis_output_ROTOR.OME.tiff'];
                                    ometiffsavename_trend = [dirname filesep obj.current_filename '_analysis_output_TREND.OME.tiff'];
                                    ometiffsavename_isochrones = [dirname filesep obj.current_filename '_analysis_output_ISOCHRONES.OME.tiff'];
                                    bfsave(fig.rotvol,ometiffsavename_rotor,'Compression','LZW','BigTiff', true,'dimensionOrder','XYCZT');
                                    if obj.HL1_create_LF_movie
                                        bfsave(fig.LFvol,ometiffsavename_trend,'Compression','LZW','BigTiff', true,'dimensionOrder','XYCZT');
                                    end
                                    if ~isempty(fig.isochrones)
                                        bfsave(fig.isochrones,ometiffsavename_isochrones,'Compression','LZW','BigTiff', true,'dimensionOrder','XYCZT');
                                    end
                                    %
                                else
                                    ometiffsavename = [dirname filesep obj.current_filename '_analysis_output.OME.tiff'];
                                    bfsave(fig,ometiffsavename,'Compression','LZW','dimensionOrder','XYCTZ');
                                end
                    end                                                                                        
        end
%-------------------------------------------------------------------------%
        function run_batch(obj,filenames, pathname,verbose,~)
            
            switch obj.problem
                            
                case {'TTO' 'CIDR' 'PR' 'HL1' 'Experimental' 'NucCyt' 'MPHG' 'Sparks' 't_dependent_Nuclei_ratio_FRET','AI_Powered_2D_SMLM_Reconstruction'}
                    
                    dirname = [];           
                    cmnxlsname = [];
                    if obj.save_analysis_output_as_OMEtiff || obj.save_analysis_output_as_xls || strcmp(obj.problem,'t_dependent_Nuclei_ratio_FRET')
                        timestamp = datestr(now,'yyyy-mm-dd HH-MM-SS');                        
                        mkdir(obj.RootDirectory,['ALYtools Analysis ' timestamp]);
                        dirname = [obj.RootDirectory filesep ['ALYtools Analysis ' timestamp]];
                        s = strsplit(pathname,filesep);
                        cmnxlsname = [dirname filesep char(s{length(s)-1}) '_joint_analysis_data.xls'];                        
                    end
                    
                    CMN = [];
                    
                    filenames = sort_nat(cellstr(filenames));

                    hw = [];
                    waitmsg = 'Loading, segmenting and analyzing files...';
                    if verbose
                        hw = waitbar(0,waitmsg);
                    end 

                    for k = 1:numel(filenames)
                        
                        obj.current_filename = char(filenames{k});        
                        obj.open_image([pathname filesep obj.current_filename]);                        
                        %
                        obj.show_data(false); % not sending to Icy
                        %
                        success = true; % we are optimists
                        try
                            if strcmp(obj.problem,'TTO')
                                [data, caption, table_name, fig] = obj.analyze_TTO;
                            elseif strcmp(obj.problem,'CIDR')
                                [data, caption, table_name, fig] = obj.analyze_CIDR;
                            elseif strcmp(obj.problem,'PR')
                                [data, caption, table_name, fig] = obj.analyze_PR;
                            elseif strcmp(obj.problem,'HL1')
                                [data, caption, table_name, fig] = obj.analyze_HL1;
                            elseif strcmp(obj.problem,'Experimental')
                                [data, caption, table_name, fig] = obj.analyze_Experimental;                                
                            elseif strcmp(obj.problem,'NucCyt')
                                [data, caption, table_name, fig] = obj.analyze_NC;
                            elseif strcmp(obj.problem,'MPHG')
                                [data, caption, table_name, fig] = obj.analyze_MPHG;
                            elseif strcmp(obj.problem,'Sparks')
                                [data, caption, table_name, fig] = obj.analyze_Sparks;
                            elseif strcmp(obj.problem,'t_dependent_Nuclei_ratio_FRET')
                                [data, caption, table_name, fig] = obj.analyze_t_dependent_Nuclei_ratio_FRET;
                                fig = obj.t_dependent_Nuclei_ratio_FRET_postprocess(fig,dirname);
                            elseif strcmp(obj.problem,'AI_Powered_2D_SMLM_Reconstruction')
                                [data, caption, table_name, fig] = obj.analyze_AI_Powered_2D_SMLM_Reconstruction;
                            end
                        catch
                            disp(['failed to analyze the file ' obj.current_filename]);
                            success = false;
                        end
                        %                  
                        if success
                            if obj.save_analysis_output_as_xls && ~isempty(data) && ~isempty(caption) && ~isempty(table_name)
                                xlsname = [dirname filesep obj.current_filename '_analysis_data.xls']; 
                                if ispc
                                    try
                                        xlswrite( xlsname,[caption; data],char(table_name) );
                                    catch
                                        disp('can not write output as xls, save as mat file instead');                                    
                                        fname = xlsname(1:length(xlsname)-4);
                                        save([fname '.mat'],'caption','data');
                                        cell2csv([fname '.csv'],[caption; data]);                                        
                                    end
                                else
                                    try
                                        xlwrite( xlsname,[caption; data],char(table_name) );
                                    catch
                                        disp('can not write output as xls, save as mat file instead');
                                        fname = xlsname(1:length(xlsname)-4);
                                        save([fname '.mat'],'caption','data');
                                        cell2csv([fname '.csv'],[caption; data]);                                      
                                    end
                                end
                            end
                            %
                            if obj.save_analysis_output_as_xls && ~isempty(data) && ~isempty(caption) && ~isempty(table_name)
                                if isempty(CMN)                        
                                    CMN = [caption; data];
                                else
                                    CMN = [CMN; data];
                                end                               
                            end
                            %
                            if obj.send_analysis_output_to_Icy && ~isempty(fig)
                                try
                                    if strcmp(obj.problem,'HL1')
                                        icy_imshow(fig.icyvol,[obj.current_filename ' - HL1 results']);                             
                                        icy_imshow(fig.rotvol,[obj.current_filename ' - HL1 dynamics, T = ' num2str(fig.T)]);
                                        if obj.HL1_create_LF_movie
                                            icy_imshow(fig.LFvol,[obj.current_filename ' - HL1 trend']);
                                        end
                                        if ~isempty(fig.isochrones)
                                            icy_imshow(fig.isochrones,[obj.current_filename ' - HL1 iscochrones']);
                                        end
                                    else % only 1 image is sent to Icy
                                        icy_imshow(fig,[obj.current_filename ' ' obj.problem ' analysis']);
                                    end
                                catch
                                end;
                            end
                            %
                            if obj.save_analysis_output_as_OMEtiff  && ~isempty(fig)   
                                if strcmp(obj.problem,'HL1')
                                    ometiffsavename = [dirname filesep obj.current_filename '_analysis_output.OME.tiff'];
                                    bfsave(fig.icyvol,ometiffsavename,'Compression','LZW','dimensionOrder','XYCTZ');
                                    %
                                    % another 2 interesting files
                                    ometiffsavename_rotor = [dirname filesep obj.current_filename '_analysis_output_ROTOR.OME.tiff'];
                                    ometiffsavename_trend = [dirname filesep obj.current_filename '_analysis_output_TREND.OME.tiff'];
                                    ometiffsavename_isochrones = [dirname filesep obj.current_filename '_analysis_output_ISOCHRONES.OME.tiff'];
                                    bfsave(fig.rotvol,ometiffsavename_rotor,'Compression','LZW','BigTiff', true,'dimensionOrder','XYCZT');
                                    if obj.HL1_create_LF_movie
                                        bfsave(fig.LFvol,ometiffsavename_trend,'Compression','LZW','BigTiff', true,'dimensionOrder','XYCZT');
                                    end
                                    if ~isempty(fig.isochrones)
                                        bfsave(fig.isochrones,ometiffsavename_isochrones,'Compression','LZW','BigTiff', true,'dimensionOrder','XYCZT');                                        
                                    end
                                    %
                                else
                                    ometiffsavename = [dirname filesep obj.current_filename '_analysis_output.OME.tiff'];
                                    bfsave(fig,ometiffsavename,'Compression','LZW','dimensionOrder','XYCTZ');
                                end
                            end                        
                        end
                        %                                                                                                
                        if ~isempty(hw), waitbar(k/numel(filenames),hw); drawnow, end 
                        
                    end     
                    
                    if obj.save_analysis_output_as_xls && ~isempty(CMN) && ~isempty(caption) && ~isempty(table_name) && ... % common excel data file
                            ~strcmp(obj.problem,'AI_Powered_2D_SMLM_Reconstruction') && ...
                            ~strcmp(obj.problem,'ImageTiling')                            
                       if ispc
                        xlswrite( cmnxlsname,CMN );
                       else
                        xlwrite( cmnxlsname,CMN );                           
                       end
                    end
                    
                    if ~isempty(hw), delete(hw), drawnow; end
                    
                
                case 'Fungus Dependent Granule Release'
                                      
                    dirname = [];           
                    cmnxlsname = [];
                    if obj.save_analysis_output_as_OMEtiff || obj.save_analysis_output_as_xls
                        timestamp = datestr(now,'yyyy-mm-dd HH-MM-SS');                        
                        mkdir(obj.RootDirectory,['ALYtools Analysis ' timestamp]);
                        dirname = [obj.RootDirectory filesep ['ALYtools Analysis ' timestamp]];
                        s = strsplit(pathname,filesep);
                        cmnxlsname = [dirname filesep char(s{length(s)-1}) '_joint_analysis_data.xls'];                        
                    end
                    
                    filenames = sort_nat(cellstr(filenames));

                    hw = [];
                    waitmsg = 'Loading, segmenting and analyzing files...';
                    if verbose
                        hw = waitbar(0,waitmsg);
                    end            
                                        
                    CMN_GRAN_XLS = [];
                    CMN_CELL_XLS = [];
                    gran_table_name = [];
                    cell_table_name = [];
                                                                                    
                    for k = 1:numel(filenames)
                        
                        obj.current_filename = char(filenames{k});
                        obj.open_image([pathname filesep obj.current_filename]);                    
                        %
                        obj.show_data(false); % not sending to Icy
                        %
                        success = true; % we are optimists
                        try
                            [datas, captions, table_names, fig] = obj.analyze_Fungus_Dependent_Granule_Release;
                        catch
                            disp(['failed to analyze the file ' obj.current_filename]);
                            success = false;
                        end
                        %                  
                        if success
                            if obj.save_analysis_output_as_xls
                                xlsname = [dirname filesep obj.current_filename '_analysis_data.xls']; 
                                for m=1:numel(datas)
                                    if ispc
                                        xlswrite( xlsname,[captions{m}; datas{m}],char(table_names{m}) );
                                    else
                                        xlwrite( xlsname,[captions{m}; datas{m}],char(table_names{m}) );
                                    end
                                end                                        
                            end
                            %
                            if obj.save_analysis_output_as_xls
                                if isempty(CMN_GRAN_XLS) 
                                    gran_table_name = char(table_names{1});
                                    cell_table_name = char(table_names{2});                                
                                    CMN_GRAN_XLS = [captions{1}; datas{1}];
                                    CMN_CELL_XLS = [captions{2}; datas{2}];
                                else
                                    CMN_GRAN_XLS = [CMN_GRAN_XLS; datas{1}];
                                    CMN_CELL_XLS = [CMN_CELL_XLS; datas{2}];
                                end                               
                            end
                            %
                            if obj.send_analysis_output_to_Icy
                                try
                                icy_imshow(fig,[obj.current_filename ' ' obj.problem ' analysis']);
                                catch
                                end
                            end
                            %
                            if obj.save_analysis_output_as_OMEtiff                    
                                ometiffsavename = [dirname filesep obj.current_filename '_analysis_output.OME.tiff'];
                                bfsave(fig,ometiffsavename,'Compression','LZW','dimensionOrder','XYCTZ');
                            end                        
                        end
                        %                                                                                                
                        if ~isempty(hw), waitbar(k/numel(filenames),hw); drawnow, end     
                    end
                    
                    if obj.save_analysis_output_as_xls % common excel data file
                      if ispc
                        xlswrite( cmnxlsname,CMN_GRAN_XLS,gran_table_name );
                        xlswrite( cmnxlsname,CMN_CELL_XLS,cell_table_name );
                      else
                        xlwrite( cmnxlsname,CMN_GRAN_XLS,gran_table_name );
                        xlwrite( cmnxlsname,CMN_CELL_XLS,cell_table_name );                          
                      end                        
                    end
                    
                    if ~isempty(hw), delete(hw), drawnow; end
                                        
            end % switch
                                    
        end        
%-------------------------------------------------------------------------%
        function [datas, captions, table_names, fig] = analyze_Fungus_Dependent_Granule_Release(obj,~,~)             
                        
            datas = [];
            captions = [];
            table_names = 'ALYtools data';
            fig = [];
            
            sgm = obj.do_Fungus_Dependent_Granule_Release_segmentation(false); % not sending to Icy
            segmented_granules = squeeze(sgm(:,:,1));
            segmented_fungus = squeeze(sgm(:,:,2));
            segmented_cells = squeeze(sgm(:,:,3));
            
            [sizeX,sizeY,~] = size(segmented_cells);
                        
            FL = bwlabel(segmented_fungus);
            GL = bwlabel(segmented_granules);
            CL = bwlabel(segmented_cells);            

            cell_stats = regionprops(CL,'Area','Centroid','EquivDiameter');    
            granule_stats = regionprops(GL,'Area','Centroid');    
            % 
            indiv_cells = zeros(sizeX,sizeY,numel(cell_stats));
            for c=1:numel(cell_stats)
                indiv_cells(:,:,c) = (CL==c);
            end

            dist_from_cells = zeros(sizeX,sizeY,numel(cell_stats));
            for c=1:numel(cell_stats)
                dist_from_cells(:,:,c) = bwdist(indiv_cells(:,:,c));
            end

            dist_from_cell_centers = zeros(sizeX,sizeY,numel(cell_stats));
            for c=1:numel(cell_stats)
                dist_from_cell_centers(:,:,c) = bwdist(bwmorph(indiv_cells(:,:,c),'thin',Inf));
            end

            gX_dist_from_cell_centers = zeros(sizeX,sizeY,numel(cell_stats));
            gY_dist_from_cell_centers = zeros(sizeX,sizeY,numel(cell_stats));
            for c=1:numel(cell_stats)
                gX_dist_from_cell_centers(:,:,c) = bwdist(bwmorph(indiv_cells(:,:,c),'thin',Inf));
                [gX_dist_from_cell_centers(:,:,c),gY_dist_from_cell_centers(:,:,c)] = gsderiv(dist_from_cell_centers(:,:,c),1,1);
            end

            dist_to_fungi = bwdist(FL);
            [gXf,gYf] = gsderiv(dist_to_fungi,1,1);

            cell_fungi_cosines = zeros(sizeX,sizeY,numel(cell_stats));
            for c=1:numel(cell_stats)
                gXc = gX_dist_from_cell_centers(:,:,c);
                gYc = gY_dist_from_cell_centers(:,:,c);    
                cell_fungi_cosines(:,:,c) = (gXc.*gXf+gYc.*gYf)./sqrt(gXc.*gXc+gYc.*gYc)./sqrt(gXf.*gXf+gYf.*gYf);
            end

            GCM = zeros(numel(granule_stats),numel(cell_stats));
            for g=1:numel(granule_stats),
                for c=1:numel(cell_stats),
                     x = fix(granule_stats(g).Centroid(2));
                     y = fix(granule_stats(g).Centroid(1));
                     % GCM(g,c) = dist_from_cell_centers(x,y,c); % wrong - fixed.. ;)        
                     GCM(g,c) = dist_from_cells(x,y,c);        
                end    
            end
                                    
            GRANULES = [];
            granimg = zeros(size(GL));                          
            gran_ref = squeeze(double(obj.imgdata(:,:,1,1,1)));                                                                                       
            
            for g=1:numel(granule_stats),                
                dists = GCM(g,:);
                nearest_cell_index = min(find(dists==min(dists(:)))); % external min is the fix for (the very rare) case of identical distances
                x = fix(granule_stats(g).Centroid(2));
                y = fix(granule_stats(g).Centroid(1));
                area = granule_stats(g).Area;
                dist_nearest_cell = dist_from_cell_centers(x,y,nearest_cell_index);
                dist_border_nearest_cell = dist_from_cells(x,y,nearest_cell_index);
                fungi_cosine = cell_fungi_cosines(x,y,nearest_cell_index);    
                %projected_on_cell = indiv_cells(x,y,nearest_cell_index);

                if isnan(fungi_cosine)                                    
                    fungi_cosine = -1;
                end
                fungi_cell_angle = 180 - acos(fungi_cosine)/acos(-1)*180;
                
                intensity_average = mean(gran_ref(GL==g));

                % threshold on border distance
                if dist_border_nearest_cell < obj.u2pix(obj.GRANULE_TO_CELL_DISTANCE_THRESHOLD)
                    granimg(GL==g) = fungi_cosine;                    
                    cur_gran.index = g;
                    cur_gran.x = x;
                    cur_gran.y = y;
                    cur_gran.area = area;
                    cur_gran.nearest_cell_index = nearest_cell_index;
                    cur_gran.dist_nearest_cell = dist_nearest_cell;
                    cur_gran.dist_border_nearest_cell = dist_border_nearest_cell;
                    cur_gran.intensity_average = intensity_average;
                    cur_gran.fungi_cosine = fungi_cosine;
                    cur_gran.fungi_cell_angle = fungi_cell_angle;
                    %
                    tofungdist = dist_to_fungi(GL==g);
                    cur_gran.minimal_distance_to_fungus = min(tofungdist(:));
                    cur_gran.mean_distance_to_fungus = mean(tofungdist(:));                                        
                    GRANULES = [GRANULES; cur_gran ];
                end
            end
            
            CELLS = [];            
            for c=1:numel(cell_stats)
                index = c;
                x = fix(cell_stats(c).Centroid(2));
                y = fix(cell_stats(c).Centroid(1));
                area = cell_stats(c).Area;                
                equiv_diameter = cell_stats(c).EquivDiameter;                
                % number of granules                
                lut = ([GRANULES.nearest_cell_index]==c);
                granules_number  = numel(lut(lut~=0));
                %
                areaS = [GRANULES.area];
                distanceS = [GRANULES.dist_nearest_cell];
                fungi_cosineS = [GRANULES.fungi_cosine];
                fungi_cell_angleS = [GRANULES.fungi_cell_angle];
                granule_minimal_distance_to_funguS = [GRANULES.minimal_distance_to_fungus];
                granule_mean_distance_to_funguS = [GRANULES.mean_distance_to_fungus];
                %
                granules_average_area = mean(areaS(lut~=0));
                granules_average_distance_to_nearest_cell = mean(distanceS(lut~=0));
                granules_average_fungus_cosine = mean(fungi_cosineS(lut~=0));
                granules_average_fungi_cell_angle = mean(fungi_cell_angleS(lut~=0));
                granules_minimal_distance_to_fungus = min(granule_minimal_distance_to_funguS(lut~=0));
                granules_mean_distance_to_fungus = mean(granule_mean_distance_to_funguS(lut~=0));                
                %
                tofungdist = dist_to_fungi(1==squeeze(indiv_cells(:,:,c)));
                minimal_distance_to_fungus = min(tofungdist(:));
                mean_distance_to_fungus = mean(tofungdist(:));
                %                                
                CELLS = [CELLS; {obj.current_filename,...                    
                    index,...
                    y,...
                    x,...
                    area,...
                    equiv_diameter,...
                    granules_number,...
                    granules_average_area,...
                    granules_average_distance_to_nearest_cell,...
                    granules_average_fungus_cosine,...
                    granules_average_fungi_cell_angle,...
                    minimal_distance_to_fungus,...
                    mean_distance_to_fungus,...
                    granules_minimal_distance_to_fungus,...
                    granules_mean_distance_to_fungus,...                   
                    }];                           
            end                        

            % prepare outputs suitable for excel..
            
            GR_EX = [];
            CL_EX = CELLS;
            
            for g=1:numel(GRANULES),
                GR_EX = [GR_EX; {obj.current_filename, ...
                    GRANULES(g).index, ...
                    GRANULES(g).y, ...
                    GRANULES(g).x, ...
                    GRANULES(g).area, ...
                    GRANULES(g).nearest_cell_index, ...
                    GRANULES(g).dist_nearest_cell, ...
                    GRANULES(g).dist_border_nearest_cell, ...
                    GRANULES(g).intensity_average, ...
                    GRANULES(g).fungi_cosine, ...
                    GRANULES(g).fungi_cell_angle, ...
                    GRANULES(g).minimal_distance_to_fungus, ...
                    GRANULES(g).mean_distance_to_fungus, ...                    
                    }];                                        
            end
                                                                                   
           capGR = {'filename', ...
               'index', ...
               'CGX', ...
               'CGY', ...
               'Area', ...
               'Cell Index', ...
               'Distance to Cell Centre', ...
               '(Outer) Distance to Cell Border', ...
               'Mean Intensity', ...
               'Fungi Cosine (Raw)', ...
               'Fungi to Cell Angle', ...
               'Min Distance to Fungus', ...
               'Mean Distance to Fungus'};
                      
           capCL = {'filename', ...
               'index', ...
               'CGX', ...
               'CGY', ...
               'Area', ...
               'EquivDiameter', ...
               '#Granules', ...
               'Granules Mean Area', ...
               'Granules Mean Distrance to cell centre', ...
               'Granules Mean Fungi Cosine', ...
               'Granules Mean Fungi to Cell Angle', ...
               'Min Distance to Fungus', ...
               'Mean Distance to Fungus', ...
               'Min Granule Distance to Fungus', ...
               'Mean Granule Distance to Fungus'};
                                          
           datas = {GR_EX CL_EX};
           captions = {capGR capCL};
           table_names = {'Granules' 'Cells'};
           
           icyvol(:,:,1,1,1) = granimg;
           icyvol(:,:,2,1,1) = segmented_fungus;
           icyvol(:,:,3,1,1) = segmented_cells;           
           fig = icyvol;
            
        end    

%-------------------------------------------------------------------------%
        function [datas, captions, table_names, fig] = analyze_TTO(obj,~,~)
            
                datas = [];
                captions = [];
                table_names = 'ALYtools data';
                fig = [];            
                            
                u = squeeze(double(obj.imgdata(:,:,1,obj.TTO_ref_channel,:)));
                
                mask = obj.do_TTO_Segmentation(false); % not showing in Icy
                                                                                                                                      
                y0 = 1;
                x0 = 1;
                [w,h] = size(u);
                
                Lmin = obj.TTO_Lmin;
                Lmax = obj.TTO_Lmax;
                Nmax = obj.TTO_Nmax;
                fL = obj.TTO_fL;
                fH = obj.TTO_fH;
                df = obj.TTO_df;
                                
                % sanity check
                if Nmax < 10 || Lmin < (fL-2) || Lmin < 32 || Lmax < (Lmin+8) || fL < fH || df > 16
                    return;
                end
                              
                FPEAKS = [];
                CONFIDENCES = [];
                
                    t = 0;
                     hw = waitbar(0,'Conducting T-tubule organization analysis.. please wait');
                    while t < Nmax   
                        x1 = x0+floor(rand*w);
                        x2 = x0+floor(rand*w);
                        y1 = y0+floor(rand*h);
                        y2 = y0+floor(rand*h);
                        l_cur = sqrt((x2-x1)^2+(y2-y1)^2);
                        if  l_cur > Lmin && l_cur < Lmax
                            [x,y] = bresenham(x1,y1,x2,y2);                            
                            line_OK = true; % we are optimists
                            for m = 1 : length(x)
                                if 0==mask(x(m),y(m))
                                    line_OK = false;
                                    break;
                                end
                            end
                            if line_OK                                
%                                 for m = 1 : length(x) % visualize lines
%                                   patch_res(x(m),y(m)) = 1;                                  
%                                 end
                                t = t + 1;
                                %                                 
                                profile = zeros(1,length(x));
                                for m = 1 : length(x)
                                    profile(m) = u(x(m),y(m));                                  
                                end
                                %
                                [~,~,fpeak,confidence] = analyze_profile_periodicity_via_Fourier(profile,fL,fH,df); % frequency in 1/pix                                
                                FPEAKS(1,t) = fpeak;
                                CONFIDENCES(1,t) = confidence;                                
                                %
                                if ~isempty(hw), waitbar(t/Nmax,hw); drawnow, end;
                            end
                        end                    
                    end 
                    if ~isempty(hw), delete(hw), drawnow; end;
                %
                texture_organization = mean(CONFIDENCES(:));

%                 thresh = quantile(CONFIDENCES(:),0.8);
%                 strange_data =[];
%                 for k=1:length(FPEAKS)
%                     if CONFIDENCES(k) > thresh
%                         strange_data = [strange_data; FPEAKS(k)];
%                     end;
%                 end                
%                 figure(); hist(strange_data,200);
                
                icyvol(:,:,1,1,1) = u;
                icyvol(:,:,2,1,1) = mask*texture_organization;
                fig = icyvol;            
                
               datas = {obj.current_filename, texture_organization};
               captions = {'filename','texture organization'};
               table_names = {'texture organization'};                                
                                
        end                
%-------------------------------------------------------------------------%
        function [datas, captions, table_names, fig] = analyze_CIDR(obj,~,~)

                if isempty(obj.CIDR_ref_channel), errordlg('no data'); return; end; % param
            
                datas = [];
                captions = [];
                table_names = 'ALYtools data';
                fig = [];
                                            
                sgm = obj.do_CIDR_Segmentation(false); % not send to Icy
                %
                %get patches
                L_vanilla = bwlabel(sgm);
                L = bwlabel(imfill(sgm,'holes'));
                patch_stats = regionprops(L,'Area','BoundingBox');    
                
                d = obj.CIDR_smoothing_scale; % will be used in 0-crossing analysis     % param
                ref = squeeze(double(obj.imgdata(:,:,1,obj.CIDR_ref_channel,:)));
                HF_ref  = ref - box_average(ref,d);

                Lmin = fix(1.25*d); % should be a bit longer than averaging size
                
                Nmax = obj.CIDR_Nmax; % param
                                
                patch_res = zeros(size(sgm));
                
                % storage for quantifiers
                LD = zeros(1,numel(patch_stats));

                % one can use this transform to quantify patches, as well?
                % entropy = entropyfiltmex(uint8(map(HF_ref,0,255)),size(HF_ref),ones(d,d),256);
                                
                PATCHES = []; % output data
                
                for p = 1 : numel(patch_stats)
                                        
                    y0 = max(1,floor(patch_stats(p).BoundingBox(1)));
                    x0 = max(1,floor(patch_stats(p).BoundingBox(2)));
                    h = floor(patch_stats(p).BoundingBox(3));
                    w = floor(patch_stats(p).BoundingBox(4));
                    
                    t = 0; % number of successful tries of line-probing
                    NZC = 0; % total number of 0-crossings
                    L_lines = 0; % total length of thrown lines
                    
                    hw = waitbar(0,['Conducting CIDR analysis, patch ' num2str(p) ', ...please wait']);
                    while t < Nmax   
                        x1 = x0+floor(rand*w);
                        x2 = x0+floor(rand*w);
                        y1 = y0+floor(rand*h);
                        y2 = y0+floor(rand*h);
                        l_cur = sqrt((x2-x1)^2+(y2-y1)^2);
                        if  l_cur >= Lmin
                            [x,y] = bresenham(x1,y1,x2,y2);                            
                            line_OK = true; % we are optimists
                            for m = 1 : length(x)
                                if p~=L_vanilla(x(m),y(m))
                                    line_OK = false;
                                    break;
                                end
                            end
                            if line_OK                                
%                                 for m = 1 : length(x) % visualize lines
%                                   patch_res(x(m),y(m)) = 1;                                  
%                                 end
                                t = t + 1;
                                %
                                % count number of zero-crossing & update stats
                                profile = zeros(1,length(x));
                                for m = 1 : length(x)
                                    profile(m) = HF_ref(x(m),y(m));                                  
                                end
                                NZC = NZC + zc(profile);
                                L_lines = L_lines + l_cur;  
                                %
                                if ~isempty(hw), waitbar(t/Nmax,hw); drawnow, end;
                            end
                        end                    
                    end
                    if ~isempty(hw), delete(hw), drawnow; end;
                    
                    % linear density for this patch
                    LD(p) = NZC/L_lines;
                
                    % if one uses "entropy"
                    % entropy_p_sample = entropy(L_vanilla == p);
                    % entropy_p_mean = mean(entropy_p_sample(:));
                    
                    PATCHES = [PATCHES; {obj.current_filename, ...
                        p, ...
                        patch_stats(p).Area, ...
                        y0,x0, ... % a bit confusing - to use Icy convention, one needs to switch axes
                        h,w, ...
                        LD(p)}];                        
%                        LD(p)},entropy_p_mean]; % if one uses "entropy"
                                                                
                end % for p = 1 : numel(patch_stats)
                
                for p = 1 : numel(patch_stats)
                    patch_res(L_vanilla==p) = LD(p);
                end                           
                
               icyvol(:,:,1,1,1) = HF_ref;
               icyvol(:,:,2,1,1) = L;                              
               icyvol(:,:,3,1,1) = patch_res;
               %icyvol(:,:,4,1,1) = entropy.*sgm;   % if one uses "entropy"
               fig = icyvol;
               
               datas = PATCHES;
               captions = {'filename','index','Area','BndRec X0','BndRec Y0','BndRec W','BndRec H','ZC Linear Density'};
               %captions = {'filename','index','Area','BndRec X0','BndRec Y0','BndRec W','BndRec H','ZC Linear Density','Entropy'}; % if one uses "entropy"
               table_names = {'PATCHES'};                                
               
        end        
%-------------------------------------------------------------------------%  
        function [datas, captions, table_names, fig] = analyze_PR(obj,~,~)

                if isempty(obj.PR_ref_channel), errordlg('no data'); return; end; 
            
                datas = [];
                captions = [];
                table_names = 'ALYtools data';
                fig = [];
                     
                [~,~,~,sC,~]=size(obj.imgdata);
                
                ref = squeeze(double(obj.imgdata(:,:,1,obj.PR_ref_channel,:)));
                coloc = squeeze(double(obj.imgdata(:,:,1,obj.PR_coloc_channel,:)));
                
                sgm = obj.do_PR_Segmentation(false); % not send to Icy
                %
                L = bwlabel(sgm);
                puncta_stats = regionprops(L,'Area','Centroid','Perimeter','Eccentricity','Orientation');    
                                
                PUNCTA = []; % output data
                
                OC_img = zeros(size(ref));
                r_img = zeros(size(ref));
                pval_img = zeros(size(ref));
                
                for p = 1 : numel(puncta_stats)
                                 
                    intensity_sample_p = ref(L==p);
                    mean_intensity_p = mean(intensity_sample_p(:));
                    
                    % in-home Colocalization etude - starts
                    A = intensity_sample_p; 
                    a = A(:);
                    B = coloc(L==p); 
                    b = B(:);
                    %
                    [r,pval] = corr(a,b);
                    % overlap coef
                    OC = sum(a.*b)/sqrt(sum(a.*a)*sum(b.*b));
                    %
                    OC_img(L==p) = OC; % overlap
                    r_img(L==p) = r; % pearson
                    pval_img(L==p) = pval; % p-value
                    % in-home Colocalization etude - ends
                                        
                    PUNCTA = [PUNCTA; {obj.current_filename, ...
                        p, ...
                        puncta_stats(p).Area, ...
                        puncta_stats(p).Centroid(2), ...
                        puncta_stats(p).Centroid(1), ...
                        puncta_stats(p).Perimeter, ...
                        puncta_stats(p).Eccentricity, ...
                        puncta_stats(p).Orientation, ...
                        mean_intensity_p, ...
                        mean_intensity_p*puncta_stats(p).Area, ...
                        OC,r,pval}];
                                                                
                end % for p = 1 : numel(puncta_stats)
                               
               icyvol(:,:,1,1,1) = ref;
               icyvol(:,:,2,1,1) = OC_img;
               icyvol(:,:,3,1,1) = r_img;
               icyvol(:,:,4,1,1) = pval_img;
               
               fig = icyvol;
               
               datas = PUNCTA;
               captions = {'filename','index','Area','CX','CY', ...
                   'Perimeter','Eccentricity','Orientation','Mean Intensity','Total Intensity','Overlap coeff','Pearson coef','p-value'};
               table_names = {'PUNCTA'};                                               
        end 
        
%-------------------------------------------------------------------------%        
        function [datas, captions, table_names, fig] = analyze_HL1(obj,~,~)
            
            if isempty(obj.HL1_ref_channel), errordlg('no data'); return; end;
            
            datas = [];
            captions = [];
            table_names = 'ALYtools data';
            fig = [];            

            [sX,sY,sZ,sC,sT] = size(obj.imgdata);
            if sT < 16, return, end; % no point
            %
            data = squeeze(obj.imgdata(:,:,1,obj.HL1_ref_channel,:));
            
            if obj.HL1_invert_input
                data = max(data(:)) - data;
            end
            
            avr_size = obj.HL1_TD_avreraging_size;
            df = obj.HL1_df;
                                
            LF = zeros(size(data));
            isochrones = [];
            rotor = zeros(size(data));
            rotor_mask = zeros(sX,sY);
            P = zeros(sX,sY);
            A = zeros(sX,sY);
            R = zeros(sX,sY);
            F0 = zeros(sX,sY);
            phasemap = zeros(sX,sY);
            org = squeeze(sum(data,3));
            %
            med_T = zeros(sX,sY);
            std_T = zeros(sX,sY);
            std_PKS_HEIGHT = zeros(sX,sY);            
            
            try
                rotor_mask = obj.do_HL1_Segmentation(false);
            catch
                errordlg('Segmentation error- can not continue');
            end

            % Fourier analysis on original data
            hw = waitbar(0,'Fourier analysis - 1) dominant frequency, 2) periodicity - please wait');                             
            for x=1:sX,
                if ~isempty(hw), waitbar(x/sX,hw); drawnow, end;   
                for y=1:sY,
                    %
                    if rotor_mask(x,y)                        
                        s = double(squeeze(data(x,y,:)));
                        fs = fft(s);
                        psd = abs(fs).^2;
                        L = length(psd);
                        psd = psd(fix(L/2)+1:L);
                        % f0, periodicity index
                        [F0(x,y),R(x,y)] = estimate_1d_psd_periodicity(psd,df,3*df); % 3*df - cutoff                                                              
                        %
                    end
                end
            end
            if ~isempty(hw), delete(hw), drawnow; end;
            
            % de-trending % smoothing
            hw = waitbar(0,'Time domain detrending, smoothing, calculating Amplitude and Power - please wait');
            detrended_data = zeros(size(data));            
            for x=1:sX,
                if ~isempty(hw), waitbar(x/sX,hw); drawnow, end;   
                for y=1:sY,
                    %
                    if rotor_mask(x,y)                        
                        %
                        s = double(squeeze(data(x,y,:)));                                                                                                                        
                        if obj.HL1_TD_smoothing_size > 0
                            %s = smooth(s,obj.HL1_TD_smoothing_size);                                                
                            [~,lf_xy] = TD_high_pass_filter(s,obj.HL1_TD_smoothing_size);
                            s = lf_xy;                            
                        end                        
                        [hf_xy,lf_xy] = TD_high_pass_filter(s,avr_size);
                            % ehm, this is faster, but there is problem on t-borders, suspect padding & 'valid' is needed
                            %lf_xy = smooth(s,avr_size);
                            %hf_xy = s - lf_xy;                        
                        LF(x,y,:)=lf_xy'/lf_xy(1);                                                                        
                        %
                        % Amplitude & Power
                        min_xy = min(hf_xy);
                        max_xy = max(hf_xy);                                
                        a_xy = max_xy - min_xy; 
                        A(x,y) = a_xy;
                        % power
                        p_xy = 1/length(hf_xy)*sum(hf_xy.*hf_xy);                   
                        P(x,y) = p_xy; 
                        %
                        detrended_data(x,y,:) = hf_xy;
                    end
                end
            end
            if ~isempty(hw), delete(hw), drawnow; end;
                                                
            % spatial binning  
            binned_detrended_data = zeros(size(data));            
            hw = waitbar(0,'Binning - please wait');
            SE = strel('disk',obj.HL1_binning_radius,0); % NB!
            se = double(SE.getnhood);
            pix_num = conv2(double(rotor_mask),se,'same'); % denominator                                   
                for t=1:sT
                    if ~isempty(hw), waitbar(t/sT,hw); drawnow, end;
                    plane_convolved = conv2(double(squeeze(detrended_data(:,:,t))),se,'same');                                                                                
                    pix_num(pix_num==0) = 1; % numerator is also 0 here, so just avoiding 0/0
                    plane_convolved = plane_convolved./pix_num;                                                            
                    %plane_convolved(~rotor_mask) = 0;
                    binned_detrended_data(:,:,t) = plane_convolved;                                        
                end
            if ~isempty(hw), delete(hw), drawnow; end;            
            clear('detrended_data');

            hw = waitbar(0,'Normalizing, constructing dynamic movie, (maybe even conducting Peak Analysis), - please wait');                 
            phasemap_mask = rotor_mask;
            first_peak_arrival_time = zeros(sX,sY);
            for x=1:sX,
                if ~isempty(hw), waitbar(x/sX,hw); drawnow, end;   
                for y=1:sY,
                    %
                    if rotor_mask(x,y)
                        s = squeeze(binned_detrended_data(x,y,:));
                        s = s - mean(s(:));
                        min_xy = min(s(:));
                        max_xy = max(s(:));                                
                        a_xy = max_xy - min_xy;                         
                        %
                        rotor_signal = (s-min_xy)/a_xy;                        
                        rotor(x,y,:) = rotor_signal;       
                        %
                        if obj.HL1_calculate_phasemap
                            z = rotor_signal;
                            [pks,locs]= findpeaks(z,'MINPEAKHEIGHT',(median(z)+obj.HL1_MINPEAKHEIGHT_std_factor*std(z)),'MINPEAKDISTANCE',obj.HL1_MINPEAKDISTANCE);
                            %
                            if ~isempty(pks)
                                peaktopeak = diff(locs);
                                med_T(x,y) = median(peaktopeak);
                                std_T(x,y) = std(peaktopeak);
                                std_PKS_HEIGHT(x,y) = std(pks);
                                first_peak_arrival_time(x,y) = locs(1);
                            else
                                phasemap_mask(x,y) = 0;
                            end
                        end
                    end
                end
            end
            if ~isempty(hw), delete(hw), drawnow; end;
            clear('binned_detrended_data');
            
            T = 0;            
            if obj.HL1_calculate_phasemap
                % measure dominant period
                sample = med_T(med_T~=0);
                N = hist(sample(:),1:sT);
                T = find(N==max(N(:)));
                                
                good_periodicity    = std_T < obj.HL1_min_std_T; % NB - this is anti-threshold, values BELOW are allowed
                good_amplitude      = A > obj.HL1_dynamic_amplitude_threshold;     
                phasemap_mask       = phasemap_mask & good_periodicity & good_amplitude;
                                
                hw = waitbar(0,'Building phase map - please wait');
                for x=1:sX,
                    if ~isempty(hw), waitbar(x/sX,hw); drawnow, end;
                    for y=1:sY,
                        if phasemap_mask(x,y)
                            arrival = first_peak_arrival_time(x,y);
                            if arrival == 0, arrival = 1; end; % mmmmm 
                            if arrival > T, arrival = 1; end;
                            phasemap(x,y) = arrival;                            
                        end
                    end
                end
                if ~isempty(hw), delete(hw), drawnow; end; % phase map             
            end
             
            % on the rotor, delineate biggest object brighter than 1/2
            if obj.HL1_delineate_wave_front
                radius = fix(min(sX,sY)/16);
                SE = strel('disk',radius,0);
                hw = waitbar(0,'Defining front location, - please wait');
                for t=1:sT
                    if ~isempty(hw), waitbar(t/sT,hw); drawnow, end;
                    z1 = squeeze(rotor(:,:,t));
                    z = z1 > 1/2;
                    %retain biggest object
                    [L, num] = bwlabel(z);
                    count_pixels_per_obj = sum(bsxfun(@eq,L(:),1:num));
                    [~,ind] = max(count_pixels_per_obj);
                    z = (L==ind);
                    % close
                    z = imclose(z,SE);
                    % fill holes
                    z = imfill(z,'holes');
                    z = imdilate(bwperim(z),ones(3,3));
                    z1(z~=0) = 0;
                    rotor(:,:,t) = z1;                
                end
                if ~isempty(hw), delete(hw), drawnow; end;
            end
            %
            if obj.HL1_create_frame_indexed_excited_regions_movie
                isochrones = zeros(size(data),'uint16');
                radius = fix(min(sX,sY)/16);
                SE = strel('disk',radius,0);
                hw = waitbar(0,'Creating frame-indexed activity movie, - please wait');
                for t=1:sT
                    if ~isempty(hw), waitbar(t/sT,hw); drawnow, end;
                    z1 = squeeze(rotor(:,:,t));
                    z = z1 > 1/2;
                    %retain biggest object
                    [L, num] = bwlabel(z);
                    count_pixels_per_obj = sum(bsxfun(@eq,L(:),1:num));
                    [~,ind] = max(count_pixels_per_obj);
                    z = (L==ind);
                    % close
                    z = imclose(z,SE);
                    % fill holes
                    z = imfill(z,'holes');
                    isochrones(:,:,t) = sT - z*t;                
                end
                if ~isempty(hw), delete(hw), drawnow; end;                               

                % oops - resample
                sIso = sT;
                if obj.H1_isochrones_map_step > 1
                    sIso = ceil(sT/obj.H1_isochrones_map_step);
                    isochrones2 = zeros([sX sY sIso],'uint16');
                    for k=1:sIso
                        isochrones2(:,:,k) = squeeze(isochrones(:,:,1+(k-1)*obj.H1_isochrones_map_step));
                    end
                    isochrones = isochrones2;
                end
                isochrones = reshape(isochrones,[sX sY 1 1 sIso]);                                 
            end
            
            icyvol(:,:,1,1,1) = org;
            icyvol(:,:,2,1,1) = A;
            icyvol(:,:,3,1,1) = P;
            icyvol(:,:,4,1,1) = F0; 
            icyvol(:,:,5,1,1) = R;
            
            if obj.HL1_calculate_phasemap
                icyvol(:,:,6,1,1) = med_T;
                icyvol(:,:,7,1,1) = std_T;
                icyvol(:,:,8,1,1) = std_PKS_HEIGHT;            
                icyvol(:,:,9,1,1) = phasemap;
            end
                        
            %rotvol = reshape(rotor,[sX sY 1 1 sT]); 
            %LFvol = reshape(LF,[sX sY 1 1 sT]);                                     
            rotvol = reshape(uint16(1000*rotor),[sX sY 1 1 sT]);             
            LFvol = reshape(uint16(1000*LF),[sX sY 1 1 sT]); 
                        
            fig.icyvol = icyvol;
            fig.rotvol = rotvol;
            fig.isochrones = isochrones;            
                if obj.HL1_create_LF_movie
                    fig.LFvol = LFvol;
                else
                    fig.LFvol = [];                
                end
            fig.T = T;        
            
        end
%-------------------------------------------------------------------------%
        function send_original_to_Icy(obj,~,~)
                    
                    [sX,sY,sZ,sC,sT] = size(obj.imgdata);
                    icyvol = reshape(obj.imgdata,[sX,sY,sC,sZ,sT]);
                   
                    try
                        icy_imshow(icyvol,[ 'Original ' obj.current_filename ]);
                    catch
                        errordlg('problem with Icy, - might be not running');
                    end;            
        end        
%-------------------------------------------------------------------------%                    
        function save_settings(obj,fname,~)        
            settings = [];
            settings.DefaultDirectory = obj.DefaultDirectory;
            settings.IcyDirectory = obj.IcyDirectory;
            settings.microns_per_pixel = obj.microns_per_pixel;            
            %            
            settings.send_analysis_output_to_Icy = obj.send_analysis_output_to_Icy;
            settings.save_analysis_output_as_OMEtiff = obj.save_analysis_output_as_OMEtiff;
            settings.save_analysis_output_as_xls = obj.save_analysis_output_as_xls;
            %
            settings.send_original_to_Icy_on_show = obj.send_original_to_Icy_on_show;
            settings.ALYtools_always_on_top = obj.ALYtools_always_on_top;            
            %
            % problem
            settings.problem = obj.problem;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% granules
            settings.g_scale = obj.g_scale;
            settings.g_threshold = obj.g_threshold;
            settings.g_smoothing = obj.g_smoothing;
            settings.g_min_area = obj.g_min_area;
            settings.g_rel_bg_scale = obj.g_rel_bg_scale;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fungus
            settings.f_scale = obj.f_scale;
            settings.f_threshold = obj.f_threshold;
            settings.f_smoothing = obj.f_smoothing;
            settings.f_min_area = obj.f_min_area;
            settings.f_rel_bg_scale = obj.f_rel_bg_scale;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells
            settings.c_scale = obj.c_scale;
            settings.c_threshold = obj.c_threshold;
            settings.c_smoothing = obj.c_smoothing;
            settings.c_min_area = obj.c_min_area;
            settings.c_rel_bg_scale = obj.c_rel_bg_scale; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells                    
            settings.GRANULE_TO_CELL_DISTANCE_THRESHOLD = obj.GRANULE_TO_CELL_DISTANCE_THRESHOLD;            
            %
            % TTO
            settings.TTO_ref_channel = obj.TTO_ref_channel;
            % segmentation
            settings.TTO_threshold = obj.TTO_threshold;
            settings.TTO_smoothing_scale = obj.TTO_smoothing_scale;
            % analysis
            settings.TTO_Lmin = obj.TTO_Lmin;
            settings.TTO_Lmax = obj.TTO_Lmax;
            settings.TTO_Nmax = obj.TTO_Nmax;
            settings.TTO_fL = obj.TTO_fL;
            settings.TTO_fH = obj.TTO_fH;
            settings.TTO_df = obj.TTO_df;
            %
            % CIDR
            settings.CIDR_ref_channel = obj.CIDR_ref_channel;
            % segmentation
            settings.CIDR_ripple_scale = obj.CIDR_ripple_scale;
            settings.CIDR_ripple_threshold = obj.CIDR_ripple_threshold;
            settings.CIDR_rough_scale = obj.CIDR_rough_scale;
            settings.CIDR_rough_threshold = obj.CIDR_rough_threshold;   
            % analysis
            settings.CIDR_smoothing_scale = obj.CIDR_smoothing_scale;
            settings.CIDR_Nmax = obj.CIDR_Nmax; 
            %
            % PR
            settings.PR_ref_channel = obj.PR_ref_channel;
            settings.PR_coloc_channel = obj.PR_coloc_channel;
            % segmentation       
            settings.PR_K = obj.PR_K;
            settings.PR_S1 = obj.PR_S1;
            settings.PR_S2 = obj.PR_S2;
            settings.PR_a = obj.PR_a;
            settings.PR_t = obj.PR_t;
            settings.PR_mode = obj.PR_mode;
            settings.PR_min_size = obj.PR_min_size; 
            %
            % HL1
            % segmentation
            settings.HL1_sgm_minimal_radius = obj.HL1_sgm_minimal_radius;
            settings.HL1_sgm_stdT_threshold = obj.HL1_sgm_stdT_threshold;
            settings.HL1_sgm_I_threshold = obj.HL1_sgm_I_threshold;
            settings.HL1_sgm_mode = obj.HL1_sgm_mode;                
            %
            settings.HL1_binning_radius = obj.HL1_binning_radius;            
            settings.HL1_ref_channel = obj.HL1_ref_channel;
            settings.HL1_TD_avreraging_size = obj.HL1_TD_avreraging_size;
            settings.HL1_df = obj.HL1_df;
            settings.HL1_min_std_T = obj.HL1_min_std_T; 
            settings.HL1_dynamic_amplitude_threshold = obj.HL1_dynamic_amplitude_threshold; 
            settings.HL1_calculate_phasemap = obj.HL1_calculate_phasemap;
            settings.HL1_TD_smoothing_size = obj.HL1_TD_smoothing_size;
            settings.HL1_invert_input = obj.HL1_invert_input;
            settings.HL1_MINPEAKDISTANCE = obj.HL1_MINPEAKDISTANCE;
            settings.HL1_MINPEAKHEIGHT_std_factor = obj.HL1_MINPEAKHEIGHT_std_factor;
            %
            settings.HL1_delineate_wave_front = obj.HL1_delineate_wave_front;
            settings.HL1_create_LF_movie = obj.HL1_create_LF_movie;       
            settings.HL1_create_frame_indexed_excited_regions_movie = obj.HL1_create_frame_indexed_excited_regions_movie;
            settings.H1_isochrones_map_step = obj.H1_isochrones_map_step;
            %            
            settings.NC_chNuc = obj.NC_chNuc;
            settings.NC_chCell = obj.NC_chCell;
            settings.NC_bckg_subtraction_proportion = obj.NC_bckg_subtraction_proportion;
            settings.NC_bckg_dilation_size = obj.NC_bckg_dilation_size;
            %
            settings.NC_cell_smoothing_radius = obj.NC_cell_smoothing_radius;
            settings.NC_cell_rg_std_factor_int = obj.NC_cell_rg_std_factor_int;
            settings.NC_cell_rg_std_factor_out = obj.NC_cell_rg_std_factor_out;
            settings.NC_cell_overpeak_ratio = obj.NC_cell_overpeak_ratio;
            settings.NC_nuc_scale = obj.NC_nuc_scale;
            settings.NC_nuc_threshold = obj.NC_nuc_threshold;
            settings.NC_nuc_smoothing = obj.NC_nuc_smoothing;
            settings.NC_nuc_min_area = obj.NC_nuc_min_area;
            settings.NC_nuc_rel_bg_scale = obj.NC_nuc_rel_bg_scale;
            settings.NC_nuc_breacking_distmap_smoothing_scale = obj.NC_nuc_breacking_distmap_smoothing_scale;
                                    
            settings.expar = obj.expar;
            
            settings.MPHG_chNuc = obj.MPHG_chNuc;
            settings.MPHG_chCell = obj.MPHG_chCell;
            settings.MPHG_chGran = obj.MPHG_chGran;
            settings.MPHG_nuc_scale = obj.MPHG_nuc_scale;
            settings.MPHG_nuc_threshold = obj.MPHG_nuc_threshold;
            settings.MPHG_nuc_smoothing_scale = obj.MPHG_nuc_smoothing_scale;
            settings.MPHG_nuc_min_area = obj.MPHG_nuc_min_area;
            settings.MPHG_nuc_rel_bg_scale = obj.MPHG_nuc_rel_bg_scale;
            settings.MPHG_cell_smoothing_radius = obj.MPHG_cell_smoothing_radius;
            settings.MPHG_cell_rg_std_factor_int = obj.MPHG_cell_rg_std_factor_int;
            settings.MPHG_cell_rg_std_factor_out = obj.MPHG_cell_rg_std_factor_out;
            settings.MPHG_cell_overpeak_ratio = obj.MPHG_cell_overpeak_ratio;   
            settings.MPHG_minimal_gran_size = obj.MPHG_minimal_gran_size;
            settings.MPHG_gran_overpeak_ratio = obj.MPHG_gran_overpeak_ratio;
            %
            settings.per_image_TCSPC_FLIM_irf = obj.per_image_TCSPC_FLIM_irf;
            settings.per_image_TCSPC_FLIM_irf_shift = obj.per_image_TCSPC_FLIM_irf_shift;
            settings.per_image_TCSPC_FLIM_rep_rate = obj.per_image_TCSPC_FLIM_rep_rate;
            settings.per_image_TCSPC_FLIM_irf_filename = obj.per_image_TCSPC_FLIM_irf_filename;            
            settings.per_image_TCSPC_FLIM_fit_model = obj.per_image_TCSPC_FLIM_fit_model;
            settings.per_image_TCSPC_FLIM_Tmin = obj.per_image_TCSPC_FLIM_Tmin;
            settings.per_image_TCSPC_FLIM_Tmax = obj.per_image_TCSPC_FLIM_Tmax;
            settings.per_image_TCSPC_FLIM_background_value = obj.per_image_TCSPC_FLIM_background_value;
            settings.per_image_TCSPC_FLIM_saturation_value = obj.per_image_TCSPC_FLIM_saturation_value;            
            settings.per_image_TCSPC_FLIM_tvb = obj.per_image_TCSPC_FLIM_tvb;
            settings.per_image_TCSPC_FLIM_tvb_filename = obj.per_image_TCSPC_FLIM_tvb_filename;  
            settings.per_image_TCSPC_FLIM_fixed_tauD = obj.per_image_TCSPC_FLIM_fixed_tauD;
            settings.per_image_TCSPC_FLIM_conv_irf_pp_69_70 = obj.per_image_TCSPC_FLIM_conv_irf_pp_69_70;
            settings.per_image_TCSPC_FLIM_irf_background = obj.per_image_TCSPC_FLIM_irf_background; 
            settings.per_image_TCSPC_FLIM_tvb_scaling = obj.per_image_TCSPC_FLIM_tvb_scaling;
            settings.per_image_TCSPC_FLIM_weights_resampling_factor = obj.per_image_TCSPC_FLIM_weights_resampling_factor;            
            settings.per_image_TCSPC_FLIM_Ref_lifetime = obj.per_image_TCSPC_FLIM_Ref_lifetime;
            settings.per_image_TCSPC_FLIM_averaging_sigma = obj.per_image_TCSPC_FLIM_averaging_sigma;
            %
            settings.FijiScriptsDirectory = obj.FijiScriptsDirectory;
            %
            settings.t_dependent_Nuclei_ratio_FRET_TIMESTEP = obj.t_dependent_Nuclei_ratio_FRET_TIMESTEP;
            settings.TrackMate_RADIUS = obj.TrackMate_RADIUS;
            settings.TrackMate_TARGET_CHANNEL = obj.TrackMate_TARGET_CHANNEL;
            settings.TrackMate_THRESHOLD = obj.TrackMate_THRESHOLD;
            settings.TrackMate_DO_MEDIAN_FILTERING = obj.TrackMate_DO_MEDIAN_FILTERING;
            settings.TrackMate_QUALITY = obj.TrackMate_QUALITY;
            settings.TrackMate_ALLOW_TRACK_SPLITTING = obj.TrackMate_ALLOW_TRACK_SPLITTING;
            settings.TrackMate_ALLOW_TRACK_MERGING = obj.TrackMate_ALLOW_TRACK_MERGING;
            settings.TrackMate_TRACK_DISPLACEMENT = obj.TrackMate_TRACK_DISPLACEMENT;                               
            settings.TrackMate_LINKING_MAX_DISTANCE = obj.TrackMate_LINKING_MAX_DISTANCE;
            settings.TrackMate_GAP_CLOSING_MAX_DISTANCE = obj.TrackMate_GAP_CLOSING_MAX_DISTANCE;
            settings.TrackMate_MAX_FRAME_GAP = obj.TrackMate_MAX_FRAME_GAP;            
                        
            settings.t_dependent_Nuclei_ratio_FRET_eps_A = obj.t_dependent_Nuclei_ratio_FRET_eps_A; 
            settings.t_dependent_Nuclei_ratio_FRET_eps_D = obj.t_dependent_Nuclei_ratio_FRET_eps_D; 
            settings.t_dependent_Nuclei_ratio_FRET_t_A  = obj.t_dependent_Nuclei_ratio_FRET_t_A ; 
            settings.t_dependent_Nuclei_ratio_FRET_t_D  = obj.t_dependent_Nuclei_ratio_FRET_t_D ; 
            settings.t_dependent_Nuclei_ratio_FRET_Q_A = obj.t_dependent_Nuclei_ratio_FRET_Q_A; 
            settings.t_dependent_Nuclei_ratio_FRET_Q_D = obj.t_dependent_Nuclei_ratio_FRET_Q_D; 
            settings.t_dependent_Nuclei_ratio_FRET_tau_A = obj.t_dependent_Nuclei_ratio_FRET_tau_A; 
            settings.t_dependent_Nuclei_ratio_FRET_tau_D = obj.t_dependent_Nuclei_ratio_FRET_tau_D; 
            settings.t_dependent_Nuclei_ratio_FRET_tau_FRET = obj.t_dependent_Nuclei_ratio_FRET_tau_FRET; 
            settings.t_dependent_Nuclei_ratio_FRET_K_AD = obj.t_dependent_Nuclei_ratio_FRET_K_AD;
            settings.t_dependent_Nuclei_ratio_FRET_K_DA = obj.t_dependent_Nuclei_ratio_FRET_K_DA;            
            settings.t_dependent_Nuclei_ratio_FRET_phi = obj.t_dependent_Nuclei_ratio_FRET_phi;
            settings.t_dependent_Nuclei_ratio_FRET_b_A = obj.t_dependent_Nuclei_ratio_FRET_b_A;  
            settings.t_dependent_Nuclei_ratio_FRET_b_D = obj.t_dependent_Nuclei_ratio_FRET_b_D; 
            %
            settings.t_dependent_Nuclei_ratio_FRET_donor_channel = obj.t_dependent_Nuclei_ratio_FRET_donor_channel;
            
            settings.ImageTiling_Ncols = obj.ImageTiling_Ncols;
            settings.ImageTiling_Nrows = obj.ImageTiling_Nrows;
            settings.ImageTiling_Ovlp_X = obj.ImageTiling_Ovlp_X;
            settings.ImageTiling_Ovlp_Y = obj.ImageTiling_Ovlp_Y;
            settings.ImageTiling_QT = obj.ImageTiling_QT;
            settings.ImageTiling_mode = obj.ImageTiling_mode;
            %
            settings.AI_Powered_2D_SMLM_Reconstruction_network = obj.AI_Powered_2D_SMLM_Reconstruction_network;
            settings.AI_Powered_2D_SMLM_Reconstruction_upscale_factor = obj.AI_Powered_2D_SMLM_Reconstruction_upscale_factor;
            settings.AI_Powered_2D_SMLM_Reconstruction_vicinity_half_width = obj.AI_Powered_2D_SMLM_Reconstruction_vicinity_half_width;
            settings.AI_Powered_2D_SMLM_Reconstruction_pixel_size = obj.AI_Powered_2D_SMLM_Reconstruction_pixel_size;
            settings.AI_Powered_2D_SMLM_Reconstruction_extraction_scale = obj.AI_Powered_2D_SMLM_Reconstruction_extraction_scale;
            settings.AI_Powered_2D_SMLM_Reconstruction_extraction_scale_ratio = obj.AI_Powered_2D_SMLM_Reconstruction_extraction_scale_ratio;
            settings.AI_Powered_2D_SMLM_Reconstruction_extraction_threshold = obj.AI_Powered_2D_SMLM_Reconstruction_extraction_threshold;
            settings.AI_Powered_2D_SMLM_Reconstruction_extraction_method = obj.AI_Powered_2D_SMLM_Reconstruction_extraction_method;
            settings.AI_Powered_2D_SMLM_Reconstruction_max_distance_to_spurious_pixl = obj.AI_Powered_2D_SMLM_Reconstruction_max_distance_to_spurious_pixl;
            settings.AI_Powered_2D_SMLM_Reconstruction_time_dependent_block_size = obj.AI_Powered_2D_SMLM_Reconstruction_time_dependent_block_size;
            settings.AI_Powered_2D_SMLM_Reconstruction_image_formation_method = obj.AI_Powered_2D_SMLM_Reconstruction_image_formation_method;
            settings.AI_Powered_2D_SMLM_Reconstruction_image_formation_scale = obj.AI_Powered_2D_SMLM_Reconstruction_image_formation_scale;
            settings.AI_Powered_2D_SMLM_Reconstruction_NA = obj.AI_Powered_2D_SMLM_Reconstruction_NA;
            settings.AI_Powered_2D_SMLM_Reconstruction_wavelength = obj.AI_Powered_2D_SMLM_Reconstruction_wavelength;
            settings.AI_Powered_2D_SMLM_Reconstruction_min_sigma = obj.AI_Powered_2D_SMLM_Reconstruction_min_sigma;
            settings.AI_Powered_2D_SMLM_Reconstruction_max_sigma = obj.AI_Powered_2D_SMLM_Reconstruction_max_sigma;            
            
            xml_write(fname,settings);
        end
%-------------------------------------------------------------------------%                        
        function load_settings(obj,fname,~)
             if exist(fname,'file') 
                [ settings, ~ ] = xml_read (fname);                                 
                obj.DefaultDirectory = settings.DefaultDirectory;  
                obj.IcyDirectory = settings.IcyDirectory;
                obj.microns_per_pixel = settings.microns_per_pixel;                
                %
                obj.send_analysis_output_to_Icy = settings.send_analysis_output_to_Icy;
                obj.save_analysis_output_as_OMEtiff = settings.save_analysis_output_as_OMEtiff;
                obj.save_analysis_output_as_xls = settings.save_analysis_output_as_xls;
                %                
                obj.send_original_to_Icy_on_show = settings.send_original_to_Icy_on_show;
                obj.ALYtools_always_on_top = settings.ALYtools_always_on_top;
                %
                % problem
                obj.problem = settings.problem;                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% granules
                obj.g_scale = settings.g_scale;
                obj.g_threshold = settings.g_threshold;
                obj.g_smoothing = settings.g_smoothing;
                obj.g_min_area = settings.g_min_area;
                obj.g_rel_bg_scale = settings.g_rel_bg_scale;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fungus
                obj.f_scale = settings.f_scale;
                obj.f_threshold = settings.f_threshold;
                obj.f_smoothing = settings.f_smoothing;
                obj.f_min_area = settings.f_min_area;
                obj.f_rel_bg_scale = settings.f_rel_bg_scale;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells
                obj.c_scale = settings.c_scale;
                obj.c_threshold = settings.c_threshold;
                obj.c_smoothing = settings.c_smoothing;
                obj.c_min_area = settings.c_min_area;
                obj.c_rel_bg_scale = settings.c_rel_bg_scale; 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells                    
                obj.GRANULE_TO_CELL_DISTANCE_THRESHOLD = settings.GRANULE_TO_CELL_DISTANCE_THRESHOLD;            
                %  
                % TTO
                obj.TTO_ref_channel = settings.TTO_ref_channel;
                % segmentation
                obj.TTO_threshold = settings.TTO_threshold;
                obj.TTO_smoothing_scale = settings.TTO_smoothing_scale;
                % analysis
                obj.TTO_Lmin = settings.TTO_Lmin;
                obj.TTO_Lmax = settings.TTO_Lmax;
                obj.TTO_Nmax = settings.TTO_Nmax;
                obj.TTO_fL = settings.TTO_fL;
                obj.TTO_fH = settings.TTO_fH;
                obj.TTO_df = settings.TTO_df;
                %
                % CIDR
                obj.CIDR_ref_channel = settings.CIDR_ref_channel;
                % segmentation
                obj.CIDR_ripple_scale = settings.CIDR_ripple_scale;
                obj.CIDR_ripple_threshold = settings.CIDR_ripple_threshold;
                obj.CIDR_rough_scale = settings.CIDR_rough_scale;
                obj.CIDR_rough_threshold = settings.CIDR_rough_threshold;   
                % analysis
                obj.CIDR_smoothing_scale = settings.CIDR_smoothing_scale;
                obj.CIDR_Nmax = settings.CIDR_Nmax;  
                %
                % PR
                obj.PR_ref_channel = settings.PR_ref_channel;
                obj.PR_coloc_channel = settings.PR_coloc_channel;
                % segmentation       
                obj.PR_K = settings.PR_K;
                obj.PR_S1 = settings.PR_S1;
                obj.PR_S2 = settings.PR_S2;
                obj.PR_a = settings.PR_a;
                obj.PR_t = settings.PR_t;
                obj.PR_mode = settings.PR_mode;
                obj.PR_min_size = settings.PR_min_size; 
                %                                
                % HL1
                % segmentation
                obj.HL1_sgm_minimal_radius = settings.HL1_sgm_minimal_radius;
                obj.HL1_sgm_stdT_threshold = settings.HL1_sgm_stdT_threshold;
                obj.HL1_sgm_I_threshold = settings.HL1_sgm_I_threshold;
                obj.HL1_sgm_mode = settings.HL1_sgm_mode;                
                %
                obj.HL1_ref_channel = settings.HL1_ref_channel;
                obj.HL1_TD_avreraging_size = settings.HL1_TD_avreraging_size;
                obj.HL1_df = settings.HL1_df;
                obj.HL1_min_std_T = settings.HL1_min_std_T; 
                obj.HL1_dynamic_amplitude_threshold = settings.HL1_dynamic_amplitude_threshold; 
                obj.HL1_calculate_phasemap = settings.HL1_calculate_phasemap;                
                obj.HL1_TD_smoothing_size = settings.HL1_TD_smoothing_size;
                obj.HL1_invert_input = settings.HL1_invert_input;
                obj.HL1_binning_radius = settings.HL1_binning_radius;                
                obj.HL1_MINPEAKDISTANCE = settings.HL1_MINPEAKDISTANCE;
                obj.HL1_MINPEAKHEIGHT_std_factor = settings.HL1_MINPEAKHEIGHT_std_factor;
                
                obj.HL1_delineate_wave_front = settings.HL1_delineate_wave_front;
                obj.HL1_create_LF_movie = settings.HL1_create_LF_movie;       
                obj.HL1_create_frame_indexed_excited_regions_movie = settings.HL1_create_frame_indexed_excited_regions_movie;
                obj.H1_isochrones_map_step = settings.H1_isochrones_map_step;
                %
                obj.NC_chNuc = settings.NC_chNuc;
                obj.NC_chCell = settings.NC_chCell;
                obj.NC_bckg_subtraction_proportion = settings.NC_bckg_subtraction_proportion;
                obj.NC_bckg_dilation_size = settings.NC_bckg_dilation_size;                
                %
                obj.NC_cell_smoothing_radius = settings.NC_cell_smoothing_radius;
                obj.NC_cell_rg_std_factor_int = settings.NC_cell_rg_std_factor_int;
                obj.NC_cell_rg_std_factor_out = settings.NC_cell_rg_std_factor_out;
                obj.NC_cell_overpeak_ratio = settings.NC_cell_overpeak_ratio;                
                obj.NC_nuc_scale = settings.NC_nuc_scale;
                obj.NC_nuc_threshold = settings.NC_nuc_threshold;
                obj.NC_nuc_smoothing = settings.NC_nuc_smoothing;
                obj.NC_nuc_min_area = settings.NC_nuc_min_area;
                obj.NC_nuc_rel_bg_scale = settings.NC_nuc_rel_bg_scale;
                obj.NC_nuc_breacking_distmap_smoothing_scale = settings.NC_nuc_breacking_distmap_smoothing_scale;
                %
                obj.MPHG_chNuc = settings.MPHG_chNuc;
                obj.MPHG_chCell = settings.MPHG_chCell;
                obj.MPHG_chGran = settings.MPHG_chGran;
                obj.MPHG_nuc_scale = settings.MPHG_nuc_scale;
                obj.MPHG_nuc_threshold = settings.MPHG_nuc_threshold;
                obj.MPHG_nuc_smoothing_scale = settings.MPHG_nuc_smoothing_scale;
                obj.MPHG_nuc_min_area = settings.MPHG_nuc_min_area;
                obj.MPHG_nuc_rel_bg_scale = settings.MPHG_nuc_rel_bg_scale;
                obj.MPHG_cell_smoothing_radius = settings.MPHG_cell_smoothing_radius;
                obj.MPHG_cell_rg_std_factor_int = settings.MPHG_cell_rg_std_factor_int;
                obj.MPHG_cell_rg_std_factor_out = settings.MPHG_cell_rg_std_factor_out;
                obj.MPHG_cell_overpeak_ratio = settings.MPHG_cell_overpeak_ratio;   
                obj.MPHG_minimal_gran_size = settings.MPHG_minimal_gran_size;
                obj.MPHG_gran_overpeak_ratio = settings.MPHG_gran_overpeak_ratio;
                                
                obj.expar = settings.expar;
                
                obj.per_image_TCSPC_FLIM_irf = settings.per_image_TCSPC_FLIM_irf;
                obj.per_image_TCSPC_FLIM_irf_shift = settings.per_image_TCSPC_FLIM_irf_shift;
                obj.per_image_TCSPC_FLIM_rep_rate = settings.per_image_TCSPC_FLIM_rep_rate;
                obj.per_image_TCSPC_FLIM_irf_filename = settings.per_image_TCSPC_FLIM_irf_filename;            
                obj.per_image_TCSPC_FLIM_fit_model = settings.per_image_TCSPC_FLIM_fit_model;
                obj.per_image_TCSPC_FLIM_Tmin = settings.per_image_TCSPC_FLIM_Tmin;
                obj.per_image_TCSPC_FLIM_Tmax = settings.per_image_TCSPC_FLIM_Tmax;
                obj.per_image_TCSPC_FLIM_background_value = settings.per_image_TCSPC_FLIM_background_value;
                obj.per_image_TCSPC_FLIM_saturation_value = settings.per_image_TCSPC_FLIM_saturation_value;                
                obj.per_image_TCSPC_FLIM_tvb = settings.per_image_TCSPC_FLIM_tvb;
                obj.per_image_TCSPC_FLIM_tvb_filename = settings.per_image_TCSPC_FLIM_tvb_filename; 
                obj.per_image_TCSPC_FLIM_fixed_tauD = settings.per_image_TCSPC_FLIM_fixed_tauD;
                obj.per_image_TCSPC_FLIM_conv_irf_pp_69_70 = settings.per_image_TCSPC_FLIM_conv_irf_pp_69_70;
                obj.per_image_TCSPC_FLIM_irf_background = settings.per_image_TCSPC_FLIM_irf_background;
                obj.per_image_TCSPC_FLIM_tvb_scaling = settings.per_image_TCSPC_FLIM_tvb_scaling;
                obj.per_image_TCSPC_FLIM_weights_resampling_factor = settings.per_image_TCSPC_FLIM_weights_resampling_factor;
                obj.per_image_TCSPC_FLIM_Ref_lifetime = settings.per_image_TCSPC_FLIM_Ref_lifetime;
                obj.per_image_TCSPC_FLIM_averaging_sigma = settings.per_image_TCSPC_FLIM_averaging_sigma;
                
                obj.FijiScriptsDirectory = settings.FijiScriptsDirectory;
                %
                obj.t_dependent_Nuclei_ratio_FRET_TIMESTEP = settings.t_dependent_Nuclei_ratio_FRET_TIMESTEP;
                obj.TrackMate_RADIUS = settings.TrackMate_RADIUS;
                obj.TrackMate_TARGET_CHANNEL = settings.TrackMate_TARGET_CHANNEL;
                obj.TrackMate_THRESHOLD = settings.TrackMate_THRESHOLD;
                obj.TrackMate_DO_MEDIAN_FILTERING = settings.TrackMate_DO_MEDIAN_FILTERING;
                obj.TrackMate_QUALITY = settings.TrackMate_QUALITY;
                obj.TrackMate_ALLOW_TRACK_SPLITTING = settings.TrackMate_ALLOW_TRACK_SPLITTING;
                obj.TrackMate_ALLOW_TRACK_MERGING = settings.TrackMate_ALLOW_TRACK_MERGING;
                obj.TrackMate_TRACK_DISPLACEMENT = settings.TrackMate_TRACK_DISPLACEMENT;
                obj.TrackMate_LINKING_MAX_DISTANCE = settings.TrackMate_LINKING_MAX_DISTANCE;
                obj.TrackMate_GAP_CLOSING_MAX_DISTANCE = settings.TrackMate_GAP_CLOSING_MAX_DISTANCE;
                obj.TrackMate_MAX_FRAME_GAP = settings.TrackMate_MAX_FRAME_GAP;            
                
                obj.t_dependent_Nuclei_ratio_FRET_eps_A = settings.t_dependent_Nuclei_ratio_FRET_eps_A; 
                obj.t_dependent_Nuclei_ratio_FRET_eps_D = settings.t_dependent_Nuclei_ratio_FRET_eps_D; 
                obj.t_dependent_Nuclei_ratio_FRET_t_A  = settings.t_dependent_Nuclei_ratio_FRET_t_A ; 
                obj.t_dependent_Nuclei_ratio_FRET_t_D  = settings.t_dependent_Nuclei_ratio_FRET_t_D ; 
                obj.t_dependent_Nuclei_ratio_FRET_Q_A = settings.t_dependent_Nuclei_ratio_FRET_Q_A; 
                obj.t_dependent_Nuclei_ratio_FRET_Q_D = settings.t_dependent_Nuclei_ratio_FRET_Q_D; 
                obj.t_dependent_Nuclei_ratio_FRET_tau_A = settings.t_dependent_Nuclei_ratio_FRET_tau_A; 
                obj.t_dependent_Nuclei_ratio_FRET_tau_D = settings.t_dependent_Nuclei_ratio_FRET_tau_D; 
                obj.t_dependent_Nuclei_ratio_FRET_tau_FRET = settings.t_dependent_Nuclei_ratio_FRET_tau_FRET; 
                obj.t_dependent_Nuclei_ratio_FRET_K_AD = settings.t_dependent_Nuclei_ratio_FRET_K_AD;
                obj.t_dependent_Nuclei_ratio_FRET_K_DA = settings.t_dependent_Nuclei_ratio_FRET_K_DA;
                
                obj.t_dependent_Nuclei_ratio_FRET_phi = settings.t_dependent_Nuclei_ratio_FRET_phi;
                obj.t_dependent_Nuclei_ratio_FRET_b_A = settings.t_dependent_Nuclei_ratio_FRET_b_A;  
                obj.t_dependent_Nuclei_ratio_FRET_b_D = settings.t_dependent_Nuclei_ratio_FRET_b_D;
                %
                obj.t_dependent_Nuclei_ratio_FRET_donor_channel = settings.t_dependent_Nuclei_ratio_FRET_donor_channel;
                %
                try
                    obj.ImageTiling_Ncols = settings.ImageTiling_Ncols;
                    obj.ImageTiling_Nrows = settings.ImageTiling_Nrows;
                    obj.ImageTiling_Ovlp_X = settings.ImageTiling_Ovlp_X;
                    obj.ImageTiling_Ovlp_Y = settings.ImageTiling_Ovlp_Y;
                    obj.ImageTiling_QT = settings.ImageTiling_QT;
                    obj.ImageTiling_mode = settings.ImageTiling_mode; 
                    %
                    obj.AI_Powered_2D_SMLM_Reconstruction_network = settings.AI_Powered_2D_SMLM_Reconstruction_network;
                    obj.AI_Powered_2D_SMLM_Reconstruction_upscale_factor = settings.AI_Powered_2D_SMLM_Reconstruction_upscale_factor;
                    obj.AI_Powered_2D_SMLM_Reconstruction_vicinity_half_width = settings.AI_Powered_2D_SMLM_Reconstruction_vicinity_half_width;
                    obj.AI_Powered_2D_SMLM_Reconstruction_pixel_size = settings.AI_Powered_2D_SMLM_Reconstruction_pixel_size;
                    obj.AI_Powered_2D_SMLM_Reconstruction_extraction_scale = settings.AI_Powered_2D_SMLM_Reconstruction_extraction_scale;
                    obj.AI_Powered_2D_SMLM_Reconstruction_extraction_scale_ratio = settings.AI_Powered_2D_SMLM_Reconstruction_extraction_scale_ratio;
                    obj.AI_Powered_2D_SMLM_Reconstruction_extraction_threshold = settings.AI_Powered_2D_SMLM_Reconstruction_extraction_threshold;
                    obj.AI_Powered_2D_SMLM_Reconstruction_extraction_method = settings.AI_Powered_2D_SMLM_Reconstruction_extraction_method;
                    obj.AI_Powered_2D_SMLM_Reconstruction_max_distance_to_spurious_pixl = settings.AI_Powered_2D_SMLM_Reconstruction_max_distance_to_spurious_pixl;
                    obj.AI_Powered_2D_SMLM_Reconstruction_time_dependent_block_size = settings.AI_Powered_2D_SMLM_Reconstruction_time_dependent_block_size;
                    obj.AI_Powered_2D_SMLM_Reconstruction_image_formation_method = settings.AI_Powered_2D_SMLM_Reconstruction_image_formation_method;
                    obj.AI_Powered_2D_SMLM_Reconstruction_image_formation_scale = settings.AI_Powered_2D_SMLM_Reconstruction_image_formation_scale;
                    obj.AI_Powered_2D_SMLM_Reconstruction_NA = settings.AI_Powered_2D_SMLM_Reconstruction_NA;
                    obj.AI_Powered_2D_SMLM_Reconstruction_wavelength = settings.AI_Powered_2D_SMLM_Reconstruction_wavelength;
                    obj.AI_Powered_2D_SMLM_Reconstruction_min_sigma = settings.AI_Powered_2D_SMLM_Reconstruction_min_sigma;
                    obj.AI_Powered_2D_SMLM_Reconstruction_max_sigma = settings.AI_Powered_2D_SMLM_Reconstruction_max_sigma;
                catch
                end  
                
             end
        end
%-------------------------------------------------------------------------%                
    function sgm = do_NC_Segmentation(obj,send_to_Icy,~)
            
            sgm = [];
            
            chNuc = obj.NC_chNuc;
            chCell = obj.NC_chCell;
            
            u1 = squeeze(double(obj.imgdata(:,:,chNuc))); % nukes
            u2 = squeeze(double(obj.imgdata(:,:,chCell))); % cells
            
            cell_smoothing_radius = obj.NC_cell_smoothing_radius;
            cell_rg_std_factor_int = obj.NC_cell_rg_std_factor_int;
            cell_rg_std_factor_out = obj.NC_cell_rg_std_factor_out;             
                        
            nuc_scale = obj.NC_nuc_scale;
            nuc_threshold = obj.NC_nuc_threshold;
            nuc_smoothing = obj.NC_nuc_smoothing;
            nuc_min_area = obj.NC_nuc_min_area;
            nuc_rel_bg_scale = obj.NC_nuc_rel_bg_scale; 

            nuc_breacking_distmap_smoothing_scale = obj.NC_nuc_breacking_distmap_smoothing_scale;
            
            hw = waitbar(0,'Segmenting cell nuclei and bodies, please wait');
            
            sgm_nukes = nth_segmentation(u1, ...
                nuc_scale, ...
                nuc_rel_bg_scale, ...
                nuc_threshold, ...
                nuc_smoothing, ...
                nuc_min_area);
            sgm_nukes(sgm_nukes~=0)=1;
            
                    waitbar(1/4,hw); drawnow;
            
            % segment cells and clumps
            [cnt,vls] = hist(u2(:),1:max(u2(:)));                        
            maxpeakval = find(cnt==max(cnt)); % ehm, better to apply findpeaks
            minval = min(u2(:));
            
            t = maxpeakval + obj.NC_cell_overpeak_ratio*(maxpeakval-minval);

            % this doesn't work - threshold too high            
            %             A = map(u2,0,1);
            %             min_u2 = min(u2(:));
            %             max_u2 = max(u2(:));
            %             level = graythresh(A);
            %             t1 = min_u2 + level*(max_u2-min_u2) ;
            %             disp([t t1]);            

            sgm_cells = u2>t; 
                        
            sgm_cells = sgm_cells | sgm_nukes; % should be OK...
            sgm_cells = segRegionGrowing(u2,sgm_cells,cell_rg_std_factor_int,cell_rg_std_factor_out);
            sgm_cells = bwmorph(imclose(sgm_cells,strel('disk',cell_smoothing_radius,0)),'clean');
            
                    waitbar(2/4,hw); drawnow;
                        
            sgm_cells_unbroken = sgm_cells;
                       
            % 
            % break nuclear clumps
            sgm_nukes = imfill(sgm_nukes,'holes');
            z2 = bwlabel(sgm_nukes);
            D = bwdist(~z2); %distance map                
               D = medfilt2(D,[nuc_breacking_distmap_smoothing_scale nuc_breacking_distmap_smoothing_scale]);
            D = -D;
            D(~z2) = -Inf;                                        
            L = watershed(D);                                                                                
            % remove background    
            stats = regionprops(L,'Area');    
            bckgind = find([stats.Area]==max([stats.Area]));
            L(L==bckgind) = 0;
            sgm_nukes = (L>0);
            %
            sgm_nukes = bwareaopen(sgm_nukes,nuc_min_area); % safety
            
                    waitbar(3/4,hw); drawnow;
            
            % break cell clumps
            z2 = bwmorph(sgm_nukes,'thicken',Inf);
            sep_lines = ~z2;
            sep_lines(~sgm_cells)=0;
            sgm_cells(sep_lines)=0;
                        
            % THIS MIGHT BE DANGEROUS - EXPERIMENTAL CORRECTION FOR MATLAB [bwmorph(.,'thicken',Inf)] - STARTS
            %{
            go along separation lines, for each pixel check # of different labels 
            in the 5x5 vicinity except 0, - if there is only one such label (internal "slit"), set that pixel to this label
            %}
            L_c = bwlabel(sgm_cells);                                    
            [sx,sy]=size(sep_lines);
            patch = zeros(sx,sy);
            for x=3:sx-2,
                for y=3:sy-2,
                    if 1==sep_lines(x,y)
                        rx = x-2:x+2;
                        ry = y-2:y+2;
                        s = L_c(rx,ry);
                        if 0~=sum(s(:)) 
                            u = unique(s(:));
                            if 2==max(size(u))
                                patch(x,y) = max(u);
                            end
                        end
                    end
                end
            end
            L_c = L_c + patch;
            sgm_cells = (L_c~=0);
            sgm_cells = imopen(sgm_cells,strel('disk',1));
            sgm_cells = bwareaopen(sgm_cells,cell_smoothing_radius*cell_smoothing_radius);            
            % THIS MIGHT BE DANGEROUS - EXPERIMENTAL CORRECTION FOR MATLAB [bwmorph(.,'thicken',Inf)] - ENDS
                                                            
            %                                                            
            % remove orphan pieces of cellular stuff..
            L_n = bwlabel(sgm_nukes);
            stats_n = regionprops(L_n,'Area','Centroid');    
            L_c = bwlabel(sgm_cells);                                    
            %
            for n = 1:numel(stats_n)
                xcn = fix(stats_n(n).Centroid(2));
                ycn = fix(stats_n(n).Centroid(1));
                cell_label = L_c(xcn,ycn);
                L_c(L_c==cell_label)=0;
            end            
            sgm_cells(L_c~=0)=0;            
            % remove orphan pieces of cellular stuff - ends            
                                                             
            sgm = cat(3,sgm_nukes,sgm_cells,sgm_cells_unbroken,u1,u2);
            
                   waitbar(4/5,hw); drawnow;                        
                   delete(hw); drawnow;
            
                if send_to_Icy                
                    try
                        icyvol(:,:,1,1,1) = u1;
                        icyvol(:,:,2,1,1) = sgm_nukes;
                        icyvol(:,:,3,1,1) = u2;
                        icyvol(:,:,4,1,1) = sgm_cells;                        
                        %
                        notification = [obj.current_filename ' - Segmentation: NC'];
                        if isempty(obj.h_Icy_segmentation_adjustment)
                            obj.h_Icy_segmentation_adjustment = icy_imshow(icyvol,notification);                    
                        else
                            icy_imshow(obj.h_Icy_segmentation_adjustment,icyvol,notification);                    
                        end
                    catch
                        errordlg('problem with Icy, - might be not running');
                    end                
                end                                   
    end                
%-------------------------------------------------------------------------%
        function sgm = do_TTO_Segmentation(obj,send_to_Icy,~)
                %
                u = squeeze(double(obj.imgdata(:,:,1,obj.TTO_ref_channel,:)));                
                % stuff mask definition - very crude
                threshold = quantile(u(:),obj.TTO_threshold);
                se = strel('disk',3,0);
                sgm = imclose(u>threshold,se);            
                %
                if send_to_Icy                
                    try
                        icyvol(:,:,1,1,1) = u;
                        icyvol(:,:,2,1,1) = sgm;
                        %
                        notification = [obj.current_filename ' - Segmentation: TTO'];
                        if isempty(obj.h_Icy_segmentation_adjustment)
                            obj.h_Icy_segmentation_adjustment = icy_imshow(icyvol,notification);                    
                        else
                            icy_imshow(obj.h_Icy_segmentation_adjustment,icyvol,notification);                    
                        end
                    catch
                        errordlg('problem with Icy, - might be not running');
                    end                
                end
        end
%-------------------------------------------------------------------------%                        
        function sgm = do_CIDR_Segmentation(obj,send_to_Icy,~)
            % 
            % ahem first one is phase contrast itlooks like
            [~,~,~,sC,~] = size(obj.imgdata);
            if 1==sC
                u = squeeze(double(obj.imgdata(:,:,1)));
            else
                u = squeeze(double(obj.imgdata(:,:,obj.CIDR_ref_channel)));
            end                        
    
                small_mask = obj.CIDR_ripple_scale; %1                         
                [gX,gY] = gsderiv(u,small_mask,1);
                %
                z1 = map(sqrt(gX.^2+gY.^2),0,1);
                %
                t = obj.CIDR_ripple_threshold; %0.05                               
                z = (z1>t); % threshold to get rippled patches
                                thick_size = obj.CIDR_rough_scale; % 6                        
                [gX,gY] = gsderiv(u,thick_size,1);
                z2 = map(sqrt(gX.^2+gY.^2),0,1);
                t2 = obj.CIDR_rough_threshold; % 0.3                              
                se = strel('disk',fix(thick_size/2),0);
                se2 = strel('disk',fix(thick_size),0);
                z2 = imdilate(imopen((z2>t2),se),se2);                 
                
                z = imclose(z &~ z2,se);
                % remove small holes - "size of maximal remaining hole"
                z = ~bwareaopen(~z,10*thick_size^2);
                % remove sticking parts with large "opening mask"
                se3 = strel('disk',fix(2*thick_size),0);
                z = imopen(z &~ z2,se3);
                % remove small stuff - "size-sieve"
                
                min_lin_size = thick_size*8; % line of that length...
                z = bwareaopen(z,min_lin_size^2);
                %
                % still many things hardcoded.. & performance not good
                %
                sgm = z;
                
            if send_to_Icy                
                try
                    icyvol(:,:,1,1,1) = u;
                    icyvol(:,:,2,1,1) = sgm;
                    
                    notification = [obj.current_filename ' - Segmentation: CIDR'];
                    if isempty(obj.h_Icy_segmentation_adjustment)
                        obj.h_Icy_segmentation_adjustment = icy_imshow(icyvol,notification);                    
                    else
                        icy_imshow(obj.h_Icy_segmentation_adjustment,icyvol,notification);                    
                    end
                catch
                    errordlg('problem with Icy, - might be not running');
                end                
            end
                        
        end % do_CIDR_segmentation
%-------------------------------------------------------------------------%
        function sgm = do_FJ_Segmentation(obj,send_to_Icy,~)            
            
            sgm = [];
                
            % u1    = squeeze(double(obj.imgdata(:,:,1,1,1)));
            u2    = squeeze(double(obj.imgdata(:,:,1,2,1))); % filaments & stuff
            u3    = squeeze(double(obj.imgdata(:,:,1,3,1))); % nuclear           
                                                          
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% nukes
            sgm_nukes = nth_segmentation(u3, ...
                200, ...
                2.5, ...
                0.05, ...
                5, ...
                20);
            sgm_nukes(sgm_nukes~=0)=1;
            
            % nonlinear GL mapping ("also known as histogram equalization") 
            x = u2;
            xmax = max(x(:)); %quantile(x(:),0.999);
            x_ = xmax/5;
            y = 2*xmax*(1./(1+exp(-x/x_))-1/2);
            %            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% clumps - starts            
            sgm_cells = nth_segmentation(y, ...
                    10, ...
                    2.5, ...
                    0.01, ...
                    2, ...
                    2);
            sgm_cells(sgm_cells~=0)=1;
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RG: relaxed
%                  std_factor_int = 0.8;
%                  std_factor_out = 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RG: restrictive
                 std_factor_int = 0.7;
                 std_factor_out = 0.7;

                sgm_cells = sgm_cells | sgm_nukes; % should be OK...
                sgm_cells = segRegionGrowing(y,sgm_cells,std_factor_int,std_factor_out);
                % 
                % we want to fill ditches..
                sgm_cells = imclose(sgm_cells,strel('disk',8,0));
                sgm_cells = imfill(sgm_cells,'holes');
                sgm_cells = imopen(sgm_cells,strel('disk',8,0));
                %
                % remove everything that is not clumps
                BAT = 1e4; %big area threshold
                sgm_cells = bwareaopen(sgm_cells,BAT);                                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% clumps - end                            
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% filaments               
               sigma1 = 1;
               sigma2 = 10;
               nmax = 8;
               n=64;
               fiber_texture = filaments_detection(u2,sigma1,sigma2,nmax,n);               
               S = 1;
               K  = 2.5;
               t = 0.045;
               nth1 = nonlinear_tophat(fiber_texture,S,K)-1;
               nth1(nth1<t)=0;
               sgm_filaments = bwmorph(nth1,'clean');                       
               
               % remove strong undesirable border "filaments" sometimes
               % delineating the border of clumps
               cells_internal_border = sgm_cells - imerode(sgm_cells,strel('disk',7,0));
               sgm_filaments = sgm_filaments.*(~cells_internal_border);
               
               % USEFUL - expanded filaments showing roughly where they are
               sgm_filaments_exp = imdilate(sgm_filaments,strel('disk',8,0));
               sgm_filaments_exp = sgm_filaments_exp & sgm_cells;
               %sgm_filaments = sgm_filaments_exp & sgm_cells &~ sgm_filaments;
                                          
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% junctions - start 
                %
                u = max(y(:))-y;                
                S = 12;
                str = zeros(size(u)); % strength
                %
                sigma = fix(S/2);
                [uxx,uxy,uyy] = gsderiv(u,sigma,2);
                %
                hw = waitbar(0,'Conducting junctions segmentation, please wait');
                for x=1:size(u,1)
                    if ~isempty(hw), waitbar(x/size(u,1),hw); drawnow, end;                            
                    for y=1:size(u,2)
                        H = [uxx(x,y) uxy(x,y); uxy(x,y) uyy(x,y)];
                        [~,D] = eig(H); 
                        upp = D(1,1); 
                        % uqq = D(2,2);                        
                        if upp < 0  
                            str(x,y) = abs(upp);
                        end;
                    end                    
                end
                if ~isempty(hw), delete(hw), drawnow; end;
                                
                sgm_junctions = map(str,0,1).*sgm_filaments_exp.*(~sgm_nukes);
                
                t = 0.07;
                sgm_junctions = sgm_junctions>t;
                % remove small garbage
                AT = 500; %area threshold
                sgm_junctions = bwareaopen(sgm_junctions,AT);
                
                % fill holes in "junctions" as prevention strategy before opening
                z = sgm_junctions;
                z = ~bwareaopen(~z,AT);
                % remove sticking parts with "opening mask"
                z = imopen(z,strel('disk',4,0));
                % remove small stuff - "size-sieve"
                z = bwareaopen(z,AT);
                sgm_junctions = z;
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% junctions - end
            
            %%%%%%%%%%
            
            %%%%%%%%%%
                                                
            sgm = cat(3,sgm_nukes,sgm_cells,sgm_filaments,sgm_junctions,u3,u2);
            
            if send_to_Icy                
                try
                    icyvol(:,:,1,1,1) = sgm_nukes;
                    icyvol(:,:,2,1,1) = sgm_cells;
                    icyvol(:,:,3,1,1) = sgm_filaments;
                    icyvol(:,:,4,1,1) = sgm_junctions;                    
                    icyvol(:,:,5,1,1) = u3;
                    icyvol(:,:,6,1,1) = u2;
                                        
                    notification = [obj.current_filename ' - Segmentation: Experimental'];
                    if isempty(obj.h_Icy_segmentation_adjustment)
                        obj.h_Icy_segmentation_adjustment = icy_imshow(icyvol,notification);                    
                    else
                        icy_imshow(obj.h_Icy_segmentation_adjustment,icyvol,notification);                    
                    end
                catch
                    errordlg('problem with Icy, - might be not running');
                end                
            end
                                
        end
        
%-------------------------------------------------------------------------%
        function sgm = do_FIBERTHICKNESS_Segmentation(obj,send_to_Icy,~) 
            
            sgm = [];
            %
            u = double(obj.imgdata);
            
% POSSIBLE ALTERNATIVE
% % % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% filaments               
% % % %                sigma1 = 2;
% % % %                sigma2 = 10;
% % % %                nmax = 10;
% % % %                n = 64;
% % % %                fiber_texture = filaments_detection(u,sigma1,sigma2,nmax,n);               
% % % %                
% % % %                sgm = fiber_texture;
% % % %                
% % % %                S = 3;
% % % %                K  = 2.5;
% % % %                t = 0.5;
% % % %                nth1 = nonlinear_tophat(u.*fiber_texture,S,K)-1;
% % % %                nth1(nth1<t)=0;
% % % %                sgm = bwmorph(nth1,'clean');                       
% % % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% filaments               
% POSSIBLE ALTERNATIVE            

                K = 2.5;
%                S1 = 6; % for Fred's thick simulated files
%                S2 = 12;
                S1 = 2;%FRED
                S2 = 7;%FRED
                %
                a = 0.5;
                t = 0.025; % threshold %FRED
                min_size = 64; %FRED
                 
                nth1 = nonlinear_tophat(u,S1,K)-1;
                nth1(nth1<0)=0;
                nth2 = nonlinear_tophat(u,S2,K)-1;
                nth2(nth2<0)=0;
                
                str1 = zeros(size(u));
                str2 = zeros(size(u));
                %
                sigma1 = fix(S1/2);
                sigma2 = fix(S2/2);
                [uxx1,uxy1,uyy1] = gsderiv(u,sigma1,2);
                [uxx2,uxy2,uyy2] = gsderiv(u,sigma2,2);
                %
                % dirty solution
                hw = waitbar(0,['Conducting Ridge segmentation, please wait']);
                for x=1:size(u,1)
                    if ~isempty(hw), waitbar(x/size(u,1),hw); drawnow, end;                            
                    for y=1:size(u,2)
                        H = [uxx1(x,y) uxy1(x,y); uxy1(x,y) uyy1(x,y)];
                        [~,D] = eig(H); 
                        upp = D(1,1); 
                        uqq = D(2,2);
                        %                                                        
                            % case 'Ridge'
                                if upp < 0 % && uqq~=0
                                    str1(x,y) = abs(upp);
                                    %str1(x,y) = uqq - upp;
                                end                                                                
                        %
                        H = [uxx2(x,y) uxy2(x,y); uxy2(x,y) uyy2(x,y)];
                        [~,D] = eig(H); 
                        upp = D(1,1); 
                        uqq = D(2,2);
                        %                                
                            % case 'Ridge'
                                if upp < 0 % && uqq~=0
                                    str2(x,y) = abs(upp);
                                    %str2(x,y) = uqq - upp;
                                end                                                                
                    end                    
                end
                if ~isempty(hw), delete(hw), drawnow; end;
                 %                             
                 % z = a*nth1.*map((str1),0,1) + (1-a)*nth2.*map((str2),0,1);
                 z = pixelwise_max(nth1.*map((str1),0,1),nth2.*map((str2),0,1));
                                                                  
                 z = z>t;
                 %
                 % remove small objects
                 L = bwlabel(z);
                 stats = regionprops(L,'Area');
                 idx = find([stats.Area] > min_size);
                 
                 sgm = ismember(L,idx);                
                 sgm(sgm>0)=1;
                 
                 sgm = imclose(sgm,strel('disk',fix(max(1,min(S1,S2)/2)),0)); 
                                                 
                if send_to_Icy                
                    try
                        icyvol(:,:,1,1,1) = u;
                        icyvol(:,:,2,1,1) = sgm;
                        %
                        notification = [obj.current_filename ' - Segmentation: PR'];
                        if isempty(obj.h_Icy_segmentation_adjustment)
                            obj.h_Icy_segmentation_adjustment = icy_imshow(icyvol,notification);                    
                        else
                            icy_imshow(obj.h_Icy_segmentation_adjustment,icyvol,notification);                    
                        end
                    catch
                        errordlg('problem with Icy, - might be not running');
                    end                
                end

        end                
%-------------------------------------------------------------------------%  also see the function below
        function [datas, captions, table_names, fig] = analyze_Experimental(obj,~,~)                                    
            %[datas, captions, table_names, fig] = obj.analyze_FJ; 
            %[datas, captions, table_names, fig] = obj.analyze_RPE;
            %[datas, captions, table_names, fig] = obj.analyze_FIBERTHICKNESS; 
            %[datas, captions, table_names, fig] = obj.analyze_OPT_Mouse_Lung;
            %[datas, captions, table_names, fig] = obj.analyze_MNEPR; 
            %[datas, captions, table_names, fig] = obj.analyze_DarkNuclei;
            %[datas, captions, table_names, fig] = obj.analyze_ImageTiling;
            [datas, captions, table_names, fig] = obj.analyze_AI_Powered_2D_SMLM_Reconstruction;            
        end 
%-------------------------------------------------------------------------%  also see the function above
        function sgm = do_Experimental_Segmentation(obj,send_to_Icy,~)                                    
            %sgm = obj.do_FJ_Segmentation(send_to_Icy); % filaments in junctions    
            %sgm = obj.do_RPE_Segmentation(send_to_Icy); % retinal pigment epithelial cells with lysosomes
            %sgm = obj.do_FIBERTHICKNESS_Segmentation(send_to_Icy); % special program to address fibers in SIM images
            %sgm = obj.do_OPT_Mouse_Lung_Segmentation(send_to_Icy);
            %sgm = obj.do_MNEPR_Segmentation(send_to_Icy);
            %sgm = obj.do_DarkNuclei_Segmentation(send_to_Icy);
            %sgm = obj.do_ImageTiling_Segmentation(send_to_Icy);
            sgm = obj.do_AI_Powered_2D_SMLM_Reconstruction_Segmentation(send_to_Icy);
        end
%-------------------------------------------------------------------------%          
function [datas, captions, table_names, fig] = analyze_FIBERTHICKNESS(obj,~,~) 
     datas = [];
     captions = [];
     table_names = 'ALYtools data';
     fig = [];
     %
     u = double(obj.imgdata);
     
     sgm = obj.do_FIBERTHICKNESS_Segmentation(false);
     
     sgm_skel = bwmorph(sgm,'skel',Inf);
     
     bp = bwmorph(sgm_skel,'branchpoints');     
     bp = imdilate(bp,strel('disk',1));
     sgm_skel(bp~=0)=0;
     sgm_skel = bwareaopen(sgm_skel,5); % HARDCODED... can do better %FRED
     
     %radius = 25;
     radius = 6; %FRED
     [sigmaIm, errIm]=estimate_fiber_thickness_via_fitting_FRED(sgm_skel,sgm,u,radius);
               
     icyvol(:,:,1,1,1) = u;
     icyvol(:,:,2,1,1) = sgm;
     icyvol(:,:,3,1,1) = sgm_skel;
     icyvol(:,:,4,1,1) = sigmaIm; % 2*sqrt(2*log(2))*sigmaIm*32; % FWHM, nanometers
     icyvol(:,:,5,1,1) = errIm;
     
     fig = icyvol;            
          
end
%-------------------------------------------------------------------------%          
 function [datas, captions, table_names, fig] = analyze_MPHG(obj,~,~) 
     datas = [];
     captions = [];
     table_names = 'ALYtools data';
     fig = [];
     
                %sgm = obj.do_MPHG_Segmentation(false);
                %             sgm = cat(3,u_cell,u_gran,u_nuc,sgm_cells,double(sgm_gran),sgm_nukes);
                sgm = obj.do_MPHG_2018_Segmentation(false);
                       
                if isempty(sgm), return, end
                
                u_cell = squeeze(sgm(:,:,1));
                u_gran = squeeze(sgm(:,:,2));
                u_nuc = squeeze(sgm(:,:,3));
                sgm_cell = squeeze(sgm(:,:,4));
                sgm_gran = squeeze(sgm(:,:,5));
                sgm_nuc = squeeze(sgm(:,:,6));
                                                                                
                % quantifiers estimating the effective distance of from lysososmes to nuclei
                nuc_dist_map = bwdist(sgm_nuc);
                cell_dist_map = bwdist(~sgm_cell);
                %icy_imshow(cell_dist_map);
                
                
                xcoords = zeros(size(sgm_nuc));
                xcoords(:,1)=1;
                xcoords = bwdist(xcoords)+1;
                ycoords = zeros(size(sgm_nuc));
                ycoords(1,:)=1;
                ycoords = bwdist(ycoords)+1;
                
                
                %%%%% CELL BY CELL QUANTIFICATION
                nuc_labs = bwlabel(sgm_nuc); 
                stats_n = regionprops(nuc_labs,'Area','Centroid');                                      
                cell_labs = bwlabel(sgm_cell);
                stats_c = regionprops(cell_labs,'Area','Centroid','PixelList','EquivDiameter');
                nuc_cell_lut = zeros(1,numel(stats_n));
                for n = 1:numel(stats_n)
                    nuc_cell_lut(n) = cell_labs(fix(stats_n(n).Centroid(2)),fix(stats_n(n).Centroid(1)));
                end
                %
                data = [];
                for n = 1:numel(stats_n)
                     cur_cell_lab = nuc_cell_lut(n);
                     if 0~=cur_cell_lab
                         cell_img = (cell_labs==cur_cell_lab)>0;
                         cur_nuclear_area = stats_n(n).Area;
                         cur_gran_area = sum(sum(sgm_gran.*cell_img));
                         cur_tot_gran_fluorescence = sum(sum(sgm_gran.*cell_img.*u_gran));
                         cur_tot_cell_fluorescence = sum(sum(cell_img.*u_cell));
                         %
                         cur_gran_distance_to_nucleus = nan;
                         cur_gran_distance_to_nucleus_weighted = nan;
                         cur_gran_distance_to_cell_border = nan;
                         cur_gran_distance_to_cell_border_weighted = nan;
                         cur_gran_Rg = nan;
                         cur_gran_Rg_w = nan;                         
                         if 0~=cur_gran_area                             
                            sample = nuc_dist_map(sgm_gran & cell_img);                             
                            cur_gran_distance_to_nucleus = mean(sample(:));
                            cur_gran_distance_to_nucleus_weighted = sum(sum(nuc_dist_map.*u_gran.*sgm_gran.*cell_img))/sum(sum(u_gran.*sgm_gran.*cell_img));
                            sample = cell_dist_map(sgm_gran & cell_img);
                            cur_gran_distance_to_cell_border = mean(sample(:));
                            cur_gran_distance_to_cell_border_weighted = sum(sum(cell_dist_map.*u_gran.*sgm_gran.*cell_img))/sum(sum(u_gran.*sgm_gran.*cell_img));
                            
                            gran_xC = sum(sum(xcoords.*cell_img.*sgm_gran))./sum(sum(cell_img.*sgm_gran));
                            gran_yC = sum(sum(ycoords.*cell_img.*sgm_gran))./sum(sum(cell_img.*sgm_gran));                                                        
                            res_x = (xcoords - ones(size(xcoords)).*gran_xC).*cell_img.*sgm_gran;
                            res_y = (ycoords - ones(size(ycoords)).*gran_yC).*cell_img.*sgm_gran;
                            cur_gran_Rg = sqrt( (sum(sum(res_x.*res_x)) + sum(sum(res_y.*res_y)) )/cur_gran_area );                            
                            cur_gran_Rg_w = sqrt( (sum(sum(res_x.*res_x.*u_gran)) + sum(sum(res_y.*res_y.*u_gran)) )/sum(sum(u_gran.*sgm_gran.*cell_img)) );                            
                            
                         end
                         %
                         cur_cell_area = stats_c(cur_cell_lab).Area;                                  
                         %
                         % cell gyration radius, etc..
                         pl = stats_c(cur_cell_lab).PixelList;
                         xvals = pl(:,1);
                         yvals = pl(:,2);
                         XC = stats_c(cur_cell_lab).Centroid(1);
                         YC = stats_c(cur_cell_lab).Centroid(2);
                         cur_cell_Rg = sqrt(sum((xvals-XC).*(xvals-XC)+(yvals-YC).*(yvals-YC))/length(xvals));
                         cur_cell_effective_diameter = stats_c(cur_cell_lab).EquivDiameter;
                         %
                         record = {obj.current_filename, ...
                         n, ...
                         stats_n(n).Centroid(1), ...
                         stats_n(n).Centroid(2), ...
                         cur_cell_lab, ...
                         cur_nuclear_area, ...
                         cur_gran_area, ...
                         cur_tot_gran_fluorescence, ...
                         cur_gran_distance_to_nucleus, ...
                         cur_gran_distance_to_nucleus_weighted, ...
                         cur_gran_distance_to_cell_border, ...
                         cur_gran_distance_to_cell_border_weighted, ...
                         cur_gran_Rg, ...
                         cur_gran_Rg_w, ...                                                  
                         cur_cell_area, ...
                         cur_cell_Rg, ...
                         cur_cell_effective_diameter, ...
                         cur_tot_cell_fluorescence};
                         %
                         data = [data; record];                                                                            
                     end
                end                                
                %%%%% CELL BY CELL QUANTIFICATION                
                                                                                                                                
                        icyvol(:,:,1,1,1) = u_nuc;
                        icyvol(:,:,2,1,1) = nuc_labs;
                        icyvol(:,:,3,1,1) = u_cell;
                        icyvol(:,:,4,1,1) = cell_labs;                        
                        icyvol(:,:,5,1,1) = u_gran;
                        icyvol(:,:,6,1,1) = sgm_gran; 
                
                fig = icyvol;            
                                  
                if ~isempty(data)
                    % check if it is possible to infer well plate info                    
                    res = parse_OpenFLIM_HCA_1(char(data(1,1)));
                    if ~isempty(res) && 5 == length(res)     
                        for k=1:size(data,1)
                            res = parse_OpenFLIM_HCA_1(char(data(k,1)));
                            well = res(1);
                            x = res(2);
                            y = res(3);
                            t = res(4);
                            z = res(5);
                            %
                            well = char(well);
                            LETTER = cellstr(well(1));
                            NUMBER = cellstr(well(2:length(well)));
                            %
                            currec = data(k,:);
                            fname = currec(1);
                            rest_of_rec = currec(2:length(currec));
                            rec = [fname LETTER NUMBER x y t z rest_of_rec];
                            datas = [datas; rec];                            
                        end
                        captions = {'filename','well_LETTER','wel_NUMBER','x','y','T','z','nuc_index','nuc_cX','nuc_cY','cell_index','nuclear_area', ...
                            'gran_area','gran_tot_fluorescence','gran_distance_to_nucleus','gran_distance_to_nucleus_weighted','gran_distance_to_cell_border','gran_distance_to_cell_border_weighted', ...
                            'gran_Rg','gran_Rg_w','cell_area','cell_Rg','cell_effective_diameter','cur_tot_cell_fluorescence'};
                        table_names = {'MPHG_CELL_BY_CELL'};
                        return;
                    else                    
                        datas = data;
                        captions = {'filename','nuc_index','nuc_cX','nuc_cY','cell_index','nuclear_area', ...
                            'gran_area','gran_tot_fluorescence','gran_distance_to_nucleus','gran_distance_to_nucleus_weighted','gran_distance_to_cell_border','gran_distance_to_cell_border_weighted', ...
                            'gran_Rg','gran_Rg_w','cell_area','cell_Rg','cell_effective_diameter','cur_tot_cell_fluorescence'};
                        table_names = {'MPHG_CELL_BY_CELL'};                                                                                                                                               
                    end
                end
 end     
%-------------------------------------------------------------------------%  
        function [datas, captions, table_names, fig] = analyze_FJ(obj,~,~)            
                datas = [];
                captions = [];
                table_names = 'ALYtools data';
                fig = [];
                     
                %[~,~,~,sC,~]=size(obj.imgdata);     
                
                sgm = obj.do_FJ_Segmentation(false);
                
                sgm_nukes = squeeze(sgm(:,:,1));
                sgm_cells = squeeze(sgm(:,:,2));
                sgm_filaments = squeeze(sgm(:,:,3));
                sgm_junctions = squeeze(sgm(:,:,4));
                nuke_ref = squeeze(sgm(:,:,5));
                filaments_ref = squeeze(sgm(:,:,6));                
                
                % ahem, only one param is calculated
                filaments_proj = sgm_junctions & sgm_filaments;
                FJI = sum(filaments_proj(:))/sum(sgm_junctions(:))*100; % % :)
                % disp([obj.current_filename ', FJI = ' num2str(FJI)]);
                
                icyvol(:,:,1,1,1) = sgm_nukes;
                icyvol(:,:,2,1,1) = sgm_cells;
                icyvol(:,:,3,1,1) = sgm_filaments;
                icyvol(:,:,4,1,1) = sgm_junctions;
                icyvol(:,:,5,1,1) = nuke_ref;
                icyvol(:,:,6,1,1) = filaments_ref;
                
                fig = icyvol;  
                
                % try to estimate fiber thickness via fitting
                [Npix,avr_width_val,std_width_val] = estimate_fiber_thickness_via_fitting(sgm_filaments,sgm_junctions,filaments_ref);
                %Npix=0;
                %avr_width_val=0;
                %std_width_val =0;
                % try to estimate fiber thickness via fitting
                
                datas = {obj.current_filename, FJI,Npix,avr_width_val,std_width_val};
                captions = {'filename','filament_junction','Npix','avr_width_val','std_width_val'};
                table_names = {'filament_junction'};                                                                                                
        end                
        
        
%-------------------------------------------------------------------------%  
        function [datas, captions, table_names, fig] = analyze_NC(obj,~,~)            
                datas = [];
                captions = [];
                table_names = 'ALYtools data';
                fig = [];
                     
                %[~,~,~,sC,~]=size(obj.imgdata);     
                    
                sgm = obj.do_NC_Segmentation(false);
                
                sgm_nukes = squeeze(sgm(:,:,1));
                sgm_cells = squeeze(sgm(:,:,2));
                sgm_cells_unbroken = squeeze(sgm(:,:,3));
                nuke_ref = squeeze(sgm(:,:,4));
                cell_ref = squeeze(sgm(:,:,5));  
                
                % disp('TODO - find nuc/cell correspondence, quantify morphology, etc..');
                
                L_n = bwlabel(sgm_nukes);
                stats_n = regionprops(L_n,cell_ref,'Area','Perimeter','EquivDiameter','Eccentricity','Centroid','MajorAxisLength','MinorAxisLength','Orientation','MeanIntensity');
                L_c = bwlabel(sgm_cells);
                stats_c = regionprops(L_c,cell_ref,'Area','Perimeter','EquivDiameter','Eccentricity','Centroid','MajorAxisLength','MinorAxisLength','Orientation','MeanIntensity');
                %                
                
                % excluded image
                cell_excl = L_c &~ L_n;
                %icy_imshow(cell_excl);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%% background -starts
                 z2 = bwmorph(sgm_cells_unbroken,'thicken',obj.NC_bckg_dilation_size);
                 %z2 = bwmorph(sgm_cells_unbroken,'thicken',Inf);
                 bckglabs = bwlabel(z2);
                 nuc_bck_lut = zeros(1,numel(stats_n));
                 nuc_bckds = zeros(1,numel(stats_n));                 
                 for n = 1:numel(stats_n)
                    nuc_bck_lut(n) = bckglabs(fix(stats_n(n).Centroid(2)),fix(stats_n(n).Centroid(1)));
                 end
                 bckglabs(sgm_cells_unbroken==1)=0;
                 for n = 1:numel(stats_n)
                    s = cell_ref(bckglabs==nuc_bck_lut(n));
                    nuc_bckds(n) = mean(s(:));
                 end                                  
                %%%%%%%%%%%%%%%%%%%%%%%%%% background -ends                
                
                DATA = [];
                
                hw = waitbar(0,'Performing NC quantification, please wait..');
                for n = 1:numel(stats_n)
                    if ~isempty(hw), waitbar(n/numel(stats_n),hw); drawnow, end;  
                    %
                    xc_n = fix(stats_n(n).Centroid(2));
                    yc_n = fix(stats_n(n).Centroid(1));                        
                    c = L_c(xc_n,yc_n);
                    if 0~=c % ??
                        %
                        % background of cell (or her clump) in cytoplasmic channel
                        s = cell_ref(bckglabs==nuc_bck_lut(n));
                        %
                        mean_bckgrnd = mean(s(:))*obj.NC_bckg_subtraction_proportion; % bckg attributed to the current cell
                        % disp([mean_bckgrnd mean_cytoplasm_intensity MeanIntensityCyt_n])
                        %
                        xc_c = fix(stats_c(c).Centroid(2));
                        yc_c = fix(stats_c(c).Centroid(1));                    
                        Area_c = stats_c(c).Area;
                        Perimeter_c = stats_c(c).Perimeter;
                        EquivDiameter_c = stats_c(c).EquivDiameter;
                        Eccentricity_c = stats_c(c).Eccentricity;
                        MajorAxisLength_c = stats_c(c).MajorAxisLength;
                        MinorAxisLength_c = stats_c(c).MinorAxisLength;
                        Orientation_c = stats_c(c).Orientation;
                        MeanIntensity_c = stats_c(c).MeanIntensity - mean_bckgrnd; % bckg corrected
                        %                    
                        Area_n = stats_n(n).Area;
                        Perimeter_n = stats_n(n).Perimeter;
                        EquivDiameter_n = stats_n(n).EquivDiameter;
                        Eccentricity_n = stats_n(n).Eccentricity;
                        MajorAxisLength_n = stats_n(n).MajorAxisLength;
                        MinorAxisLength_n = stats_n(n).MinorAxisLength;
                        Orientation_n = stats_n(n).Orientation;
                        MeanIntensityCyt_n = stats_n(n).MeanIntensity - mean_bckgrnd; % bckg corrected
                        % calculated
                        SF_n = Perimeter_n^2/4/pi/Area_n; 
                        SF_c = Perimeter_c^2/4/pi/Area_c; 
                        % cytoplasm intensity
                        sample = cell_ref(cell_excl & L_c==c);
                        mean_cytoplasm_intensity = mean(sample(:)) - mean_bckgrnd; % bckg corrected
                        tot_cytoplasm_intensity = sum(sample(:)) - mean_bckgrnd*numel(sample(:)); % bckg corrected
                        std_cytoplasm_intensity = std(sample(:));
                        %summary nuclear intensity in cytoplasm channel
                        tot_nuke_intensity  = MeanIntensityCyt_n*Area_n; % bckg corrected
                        nuc_cyt_ratio_tots = tot_nuke_intensity/tot_cytoplasm_intensity; % bckg corrected
                        nuc_cyt_ratio_means = MeanIntensityCyt_n/mean_cytoplasm_intensity; % bckg corrected

                        record = {obj.current_filename, ...
                            n, ...
                        xc_c, ...
                        yc_c, ...            
                        Area_c, ...
                        Perimeter_c, ...
                        EquivDiameter_c, ...
                        Eccentricity_c, ...
                        MajorAxisLength_c, ...
                        MinorAxisLength_c, ...
                        Orientation_c, ...
                        MeanIntensity_c, ...                    
                        xc_n, ...
                        yc_n, ...                    
                        Area_n, ...
                        Perimeter_n, ...
                        EquivDiameter_n, ...
                        Eccentricity_n, ...
                        MajorAxisLength_n, ...
                        MinorAxisLength_n, ...
                        Orientation_n, ...
                        MeanIntensityCyt_n, ...
                        SF_n, ...
                        SF_c, ...
                        mean_cytoplasm_intensity, ...
                        tot_cytoplasm_intensity, ...
                        std_cytoplasm_intensity, ...
                        tot_nuke_intensity, ...
                        nuc_cyt_ratio_tots, ...
                        nuc_cyt_ratio_means};

                        DATA = [DATA; record]; 
                    else
                        disp([n c]);
                    end
                    
                end            
                if ~isempty(hw), delete(hw), drawnow; end;
                
                datas = DATA;
                table_names = {''};
                captions = {'filename','nuc index', ...
                    'xc_c', ...
                    'yc_c', ...            
                    'Area_c', ...
                    'Perimeter_c', ...
                    'EquivDiameter_c', ...
                    'Eccentricity_c', ...
                    'MajorAxisLength_c', ...
                    'MinorAxisLength_c', ...
                    'Orientation_c', ...
                    'MeanIntensity_c', ...
                    'xc_n', ...
                    'yc_n', ...                    
                    'Area_n', ...
                    'Perimeter_n', ...
                    'EquivDiameter_n', ...
                    'Eccentricity_n', ...
                    'MajorAxisLength_n', ...
                    'MinorAxisLength_n', ...
                    'Orientation_n', ...
                    'MeanIntensity_n', ...
                    'SF_n', ...
                    'SF_c', ...
                    'mean_cytoplasm_intensity', ...
                    'tot_cytoplasm_intensity', ...
                    'std_cytoplasm_intensity', ...
                    'tot_nuke_intensity', ...
                    'nuc_cyt_ratio_tots', ...
                    'nuc_cyt_ratio_means'};                                                                                                                
                
               %fig
               icyvol(:,:,1,1,1) = nuke_ref;
               icyvol(:,:,2,1,1) = cell_ref;                              
               icyvol(:,:,3,1,1) = L_n;
               icyvol(:,:,4,1,1) = L_c;   
               %
               fig = icyvol;                                                
        end                
%-------------------------------------------------------------------------%
        function sgm = do_HL1_Segmentation(obj,send_to_Icy,~)
                 sgm = [];
                                 
            [sX,sY,sZ,sC,sT] = size(obj.imgdata);
            if sT < 16, return, end; % no point
            %
            data = [];
   
            threshold2 = obj.HL1_sgm_stdT_threshold;
            threshold1 = obj.HL1_sgm_I_threshold; % multiplier for golbal sum.img "std"        
            
            mode = obj.HL1_sgm_mode; 
            
            if 0 ~= obj.HL1_sgm_minimal_radius
                se = strel('disk',obj.HL1_sgm_minimal_radius,0);                                         
            end
            
            if 1==sC
                data = squeeze(obj.imgdata(:,:,1,1,:));
            else
                data = squeeze(obj.imgdata(:,:,1,obj.HL1_ref_channel,:));
            end
            
            %
            sgm1 = zeros(sX,sY);
            sgm2 = zeros(sX,sY);
            
            if ~strcmp(mode,'Intensity std only')            
                hw = waitbar(0,'Creating background mask - please wait');
                tQ = zeros(sX,sY); % time qualifier
                for x=1:sX,
                    if ~isempty(hw), waitbar(x/sX,hw); drawnow, end;   
                    for y=1:sY,                    
                            s = squeeze(double(data(x,y,:)));
                            std_s = std(s(:));
                            tQ(x,y) = std_s;
                    end
                end
                if ~isempty(hw), delete(hw), drawnow; end;
                z = bwmorph(tQ>threshold2,'clean');
                if 0 ~= obj.HL1_sgm_minimal_radius                
                    sgm2 = imclose(z,se);                
                else
                    sgm2 = z;                
                end
            end

            sumorg = sum(data,3);            
            
            if ~strcmp(mode,'Time std only')            
                z = sumorg>(threshold1*std(sumorg(:)));                        
                z = bwmorph(z,'clean');
                if 0 ~= obj.HL1_sgm_minimal_radius                
                    sgm1 = imclose(z,se);                
                else
                    sgm1 = z;                
                end

            end

            switch mode
                case 'Intensity std only'
                    sgm = sgm1;
                case 'Time std only'
                    sgm = sgm2;
                case 'Intensity AND Time std'
                    sgm = sgm1 | sgm2;
                case 'Intensity OR Time std'
                    sgm = sgm1 & sgm2;
            end
                        
            %leave biggest object
            if 0~=sum(sgm(:))
                [L, num] = bwlabel(sgm);
                count_pixels_per_obj = sum(bsxfun(@eq,L(:),1:num));
                [~,ind] = max(count_pixels_per_obj);
                sgm = (L==ind);
            end
            
            if obj.HL1_sgm_fill_holes
                sgm = imfill(sgm,'holes');
            end
            
            %draw 1-pix frame;
            sgm(1:sX,sY)=0;
            sgm(1:sX,1)=0;
            sgm(sX,1:sY)=0;
            sgm(1,1:sY)=0;
                        
            if send_to_Icy                
                try
                    icyvol(:,:,1,1,1) = sumorg;
                    icyvol(:,:,2,1,1) = sgm;                    
                    
                    notification = [obj.current_filename ' - Segmentation: HL1'];
                    if isempty(obj.h_Icy_segmentation_adjustment)
                        obj.h_Icy_segmentation_adjustment = icy_imshow(icyvol,notification);                    
                    else
                        icy_imshow(obj.h_Icy_segmentation_adjustment,icyvol,notification);                    
                    end
                catch
                    errordlg('problem with Icy, - might be not running');
                end                
            end
        end
%-------------------------------------------------------------------------%
        function sgm = do_PR_Segmentation(obj,send_to_Icy,~)
                %
                u = squeeze(double(obj.imgdata(:,:,1,obj.PR_ref_channel,:)));
                %                          
                K = obj.PR_K;
                S1 = obj.PR_S1;
                S2 = obj.PR_S2;
                a = obj.PR_a;
                t = obj.PR_t;
                mode = obj.PR_mode;
                min_size = obj.PR_min_size;
                 
                nth1 = nonlinear_tophat(u,S1,K)-1;
                nth1(nth1<0)=0;
                nth2 = nonlinear_tophat(u,S2,K)-1;
                nth2(nth2<0)=0;
                
                str1 = zeros(size(u));
                str2 = zeros(size(u));
                %
                sigma1 = fix(S1/2);
                sigma2 = fix(S2/2);
                [uxx1,uxy1,uyy1] = gsderiv(u,sigma1,2);
                [uxx2,uxy2,uyy2] = gsderiv(u,sigma2,2);
                %
                % dirty solution
                hw = waitbar(0,['Conducting ' obj.PR_mode ' segmentation, please wait']);
                for x=1:size(u,1)
                    if ~isempty(hw), waitbar(x/size(u,1),hw); drawnow, end;                            
                     for y=1:size(u,2)
                        H = [uxx1(x,y) uxy1(x,y); uxy1(x,y) uyy1(x,y)];
                        [~,D] = eig(H); 
                        upp = D(1,1); 
                        uqq = D(2,2);
                        %
                        switch mode 
                            case 'Ridge'
                                if upp < 0
                                    str1(x,y) = abs(upp);
                                end
                            case 'Peak'
                            if upp < 0 && uqq < 0
                                str1(x,y) = sqrt(upp*uqq);
                            end             
                        end
                        %
                        H = [uxx2(x,y) uxy2(x,y); uxy2(x,y) uyy2(x,y)];
                        [~,D] = eig(H); 
                        upp = D(1,1); 
                        uqq = D(2,2);
                        %
                        switch mode 
                            case 'Ridge'
                                if upp < 0
                                    str2(x,y) = abs(upp);
                                end
                            case 'Peak'
                            if upp < 0 && uqq < 0
                                str2(x,y) = sqrt(upp*uqq);
                            end             
                        end                                                         
                    end                    
                end
                if ~isempty(hw), delete(hw), drawnow; end;
                %                             
                z = a*nth1.*map((str1),0,1) + (1-a)*nth2.*map((str2),0,1);
                z = imclose(z>t,strel('disk',fix(max(1,min(S1,S2)/2)),0));
                %
                % remove small objects
                L = bwlabel(z);
                stats = regionprops(L,'Area');
                idx = find([stats.Area] > min_size);
                
                sgm = ismember(L,idx);                
                sgm(sgm>0)=1;
                
                if send_to_Icy                
                    try
                        icyvol(:,:,1,1,1) = u;
                        icyvol(:,:,2,1,1) = sgm;
                        %
                        notification = [obj.current_filename ' - Segmentation: PR'];
                        if isempty(obj.h_Icy_segmentation_adjustment)
                            obj.h_Icy_segmentation_adjustment = icy_imshow(icyvol,notification);                    
                        else
                            icy_imshow(obj.h_Icy_segmentation_adjustment,icyvol,notification);                    
                        end
                    catch
                        errordlg('problem with Icy, - might be not running');
                    end                
                end
        end
%-------------------------------------------------------------------------%                
    function sgm = do_MPHG_Segmentation(obj,send_to_Icy,~)
            
            sgm = [];
                                   
            u_nuc = sum(double(obj.imgdata(:,:,:,obj.MPHG_chNuc,1)),3); 
            u_cell = sum(double(obj.imgdata(:,:,:,obj.MPHG_chCell,1)),3); 
            u_gran = sum(double(obj.imgdata(:,:,:,obj.MPHG_chGran,1)),3); 
             
            % segment nukes            
            nuc_scale = fix(obj.u2pix(obj.MPHG_nuc_scale));
            nuc_threshold = obj.MPHG_nuc_threshold;
            nuc_smoothing = fix(obj.u2pix(obj.MPHG_nuc_smoothing_scale));
            nuc_min_area = fix(obj.u22pix2(obj.MPHG_nuc_min_area));
            nuc_rel_bg_scale = obj.MPHG_nuc_rel_bg_scale;
            
            % segment cells and clumps            
            cell_smoothing_radius = fix(obj.u2pix(obj.MPHG_cell_smoothing_radius));
            cell_rg_std_factor_int = obj.MPHG_cell_rg_std_factor_int;
            cell_rg_std_factor_out = obj.MPHG_cell_rg_std_factor_out;
            cell_overpeak_ratio = obj.MPHG_cell_overpeak_ratio;               
             
            % segment granules
            minimal_gran_size = fix(obj.u2pix(obj.MPHG_minimal_gran_size));
            gran_overpeak_ratio = obj.MPHG_gran_overpeak_ratio;
                        
            sgm_nukes = nth_segmentation(u_nuc, ...
                nuc_scale, ...
                nuc_rel_bg_scale, ...
                nuc_threshold, ...
                nuc_smoothing, ...
                nuc_min_area);
            sgm_nukes(sgm_nukes~=0)=1;                                                
            % segment nukes                         
                                                                                       
            % cell image histo analysis, to find proper threshold
            minval = min(u_cell(:));
            
            [cnt,vls] = hist(u_cell(:),100);
            [pks,locs]= findpeaks(cnt,'MINPEAKDISTANCE',10);
            %
            candidate_index = min(locs);
            %
            if cnt(1) > cnt(candidate_index)
                candidate_index = 1;
            end
            %
            if candidate_index > 15
                t = vls(candidate_index); % or, one can try min on inverted histo..
            else 
                t = vls(candidate_index) + cell_overpeak_ratio*(vls(candidate_index)-minval);
            end            
            % cell image histo analysis, to find proper threshold
            
            sgm_cells = u_cell > t; 
                        
            sgm_cells = sgm_cells | sgm_nukes; % should be OK...
            sgm_cells = segRegionGrowing(u_cell,sgm_cells,cell_rg_std_factor_int,cell_rg_std_factor_out);
            sgm_cells = bwmorph(imclose(sgm_cells,strel('disk',cell_smoothing_radius,0)),'clean');            
            %
            % .. that was pre-segmentation
            
            % watershed-based segmentation          
            I = u_cell;            
            % morphological smoothing - worse than simlpe Gaussian

            %             H = 20; % GL height
            %             R = 8; % radius;
            %             SE = strel('ball',R,H,0);
            %             smthd2 = imdilate(imerode(I,SE),SE);
            
            smthd2 = gsderiv(I,cell_smoothing_radius,0);

            % ONE (attempt to make it smarter.. failed?) see below...
            % mask_em = imextendedmax(smthd2, 30);
            %             mask_em = nth_segmentation(smthd2, ...
            %                 60, ...
            %                 2.5, ...
            %                 2.5, ...
            %                 4, ...
            %                 100);
            %             mask_em(mask_em~=0)=1;                                                
            % mask_em = (mask_em & sgm_cells) | sgm_nukes;                                   

            % TWO
            mask_em = sgm_nukes;
            
            % main course..
            I_c = imcomplement(smthd2);
            I_mod = imimposemin(I_c, ~sgm_cells | mask_em);
            %icy_imshow(~sgm_cells | mask_em);
            % ???
                        
            L = watershed(I_mod);
            L(~sgm_cells)=0;

            L = bwlabel(bwareaopen(L,4*nuc_min_area));

            sgm_cells1 = (L>0);

            % remove orphan pieces of cellular stuff..
            L_n = bwlabel(sgm_nukes);
            stats_n = regionprops(L_n,'Area','Centroid');    
            L_c = bwlabel(sgm_cells1);
            %
            for n = 1:numel(stats_n)
                xcn = fix(stats_n(n).Centroid(2));
                ycn = fix(stats_n(n).Centroid(1));
                cell_label = L_c(xcn,ycn);
                L_c(L_c==cell_label)=0;
            end            
            sgm_cells1(L_c~=0)=0;
            % remove orphan pieces of cellular stuff - ends                        
                                                
            % icy_imshow(label2rgb(bwlabel(sgm_cells1)));
            sgm_cells = sgm_cells1; % first iteration.. (and last.. Feb 2)
            % segment cells and clumps    
            
            % segment pathogen stuff
                 % histogram analysis to find threshold t
                    z2 = u_gran(sgm_cells~=0); % discarding empty spaces                                                                                         
                    minval = min(z2(:));
                    
                    % h = 2?IQR?n(-1/3) N - (max-min)/h %Freedman-Diaconis
%                     s = z2(:);
%                     h = 2*iqr(s)*numel(s).^(-1/3);
%                     Nbins = ceil( range(s)/h ); % damn too many
                    %                                        
                    Nbins = 100;
                    [cnt,vls] = hist(z2(:),Nbins);                    
                    
                    [pks,locs]= findpeaks(cnt,'MINPEAKDISTANCE',10);
                    %
                    candidate_index = min(locs);
                    %
                    if cnt(1) > cnt(candidate_index)
                        candidate_index = 1;
                    end
                    %
                    if candidate_index > 15; % ?? can happen if very crowded only...
                        t = vls(candidate_index); % or, one can try min on inverted histo..
                    else 
                        t = vls(candidate_index) + gran_overpeak_ratio*(vls(candidate_index)-minval);
                    end                             
                 % histogram analysis to find threshold t - ends
                 sgm_gran = u_gran > t; 
                 sgm_gran = imopen(sgm_gran,strel('disk',1));
                 sgm_gran = bwareaopen(sgm_gran,minimal_gran_size); 
                 
                 sgm_gran = sgm_gran & sgm_cells; % &~ sgm_nukes;                                               
            % segment pathogen stuff                 
                                                                    
            sgm = cat(3,u_cell,u_gran,u_nuc,sgm_cells,double(sgm_gran),sgm_nukes);
            
                if send_to_Icy                
                    try
                        icyvol(:,:,1,1,1) = u_cell;
                        icyvol(:,:,2,1,1) = u_gran;
                        icyvol(:,:,3,1,1) = u_nuc;                        
                        icyvol(:,:,4,1,1) = sgm_cells;
                        icyvol(:,:,5,1,1) = sgm_nukes;
                        icyvol(:,:,6,1,1) = double(sgm_gran);
                        %
                        notification = [obj.current_filename ' - Segmentation: MPHG'];
                        if isempty(obj.h_Icy_segmentation_adjustment)
                            obj.h_Icy_segmentation_adjustment = icy_imshow(icyvol,notification);
                        else
                            icy_imshow(obj.h_Icy_segmentation_adjustment,icyvol,notification);
                        end
                    catch
                        errordlg('problem with Icy, - might be not running');
                    end                
                end                                   
    end                
%-------------------------------------------------------------------------%                
    function sgm = do_RPE_Segmentation(obj,send_to_Icy,~)
            
            sgm = [];
            
            chNuc = 1;
            chPhaloidin = 5;
            chLysosomes = 4;
            
            un = squeeze(double(obj.imgdata(:,:,chNuc))); % nukes
            up = squeeze(double(obj.imgdata(:,:,chPhaloidin))); % fibers
            ul = squeeze(double(obj.imgdata(:,:,chLysosomes))); % lysosomes
            
            sgm_nukes = zeros(size(un));
            sgm_phal = zeros(size(un));
            sgm_lys = zeros(size(un));
   
            %%%%%%%%%%%%%%
           hw = waitbar(0,'Segmenting nuclei, cell bodies and lysosomes, - please wait');            
            %%%%%%%%%%%%%%                        
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% nuclei
            nuc_scale = 70;
            nuc_threshold = 0.05;
            nuc_smoothing = 16;
            nuc_min_area = 900;
            nuc_rel_bg_scale = 2.5; 
            %
            z2 = gsderiv(un,2,0);
            sgm_nukes = nth_segmentation(z2, ...
                nuc_scale, ...
                nuc_rel_bg_scale, ...
                nuc_threshold, ...
                nuc_smoothing, ...
                nuc_min_area);
            sgm_nukes(sgm_nukes~=0)=1;
            sgm_nukes = imfill(sgm_nukes,'holes');            
                       
            % break nuclear clumps
            nuc_breacking_distmap_smoothing_scale = 16;            
            z2 = bwlabel(sgm_nukes);
            D = bwdist(~z2); %distance map                
               D = medfilt2(D,[nuc_breacking_distmap_smoothing_scale nuc_breacking_distmap_smoothing_scale]);
            D = -D;
            D(~z2) = -Inf;                                        
            L = watershed(D);                                                                                
            % remove background    
            stats = regionprops(L,'Area');    
            bckgind = find([stats.Area]==max([stats.Area]));
            L(L==bckgind) = 0;
            sgm_nukes = (L>0);                                                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% nuclei            
           waitbar(1/5,hw); drawnow;
                        
            %%%%%%%%%%%%%%%%%% filaments                        
            
                    %%%%%%%%%%%%%%%%
        %             smoothing_radius = 4;
        %             min_nuc_size = 625;            
        %             phal_thresh = 24;
        %             sgm_phal = up>phal_thresh;    
        %             sgm_phal = imopen(sgm_phal,strel('disk',smoothing_radius));
        %             sgm_phal = bwareaopen(sgm_phal,min_nuc_size);             
                    %%%%%%%%%%%%%%%%
            
               sigma1 = 2;
               sigma2 = 10;
               nmax = 8;
               n=128;
               fiber_texture = filaments_detection(up,sigma1,sigma2,nmax,n);               
               S = 2;
               K  = 2.5;
               t = 0.08;
               nth1 = nonlinear_tophat(fiber_texture,S,K)-1;
               nth1(nth1<t)=0;
               sgm_phal = bwmorph(nth1,'clean');   
               
               % USEFUL - expanded filaments showing roughly where they are
               sgm_phal_exp = imdilate(sgm_phal,strel('disk',6,0));
               min_size = nuc_min_area*1.5;
               sgm_phal_exp = bwareaopen(sgm_phal_exp,min_size);
               % should fill in small holes... !!!!!
               % remove small holes - "size of maximal remaining hole"
               sgm_phal_exp = ~bwareaopen(~sgm_phal_exp,min_size);               
               %
               sgm_phal = sgm_phal & sgm_phal_exp;
               %
               % skeletonize phalloidin
               % sgm_phal = bwmorph(sgm_phal,'skel',Inf); % thin looks
               % better
               sgm_phal = bwmorph(sgm_phal,'thin',Inf);
               
            %%%%%%%%%%%%%%%%%% filaments                                       
           waitbar(2/5,hw); drawnow;                        
                        
            %%%%%%%%%%%%%%%%%% LYSOSOSMES
                 % histogram analysis to find threshold t
                 
                 % first version - was working not that bad..
                 %z2 = ul(sgm_phal_exp~=0); % discarding empty spaces
                 %sgm_lys = ul > quantile(z2(:),0.7);
                                  
                    z2 = ul(sgm_phal_exp~=0); % discarding empty spaces                                                                                         
                    minval = min(z2(:));
                    [cnt,vls] = hist(z2(:),100);
                    
                    % there is one gorgeous method optimizing #bins
                    % http://176.32.89.45/~hideaki/res/histogram.html -
                    % Shimazaki-Shinomoto (Matlab implementation a well?..                    
                    % http://www.mathworks.com/matlabcentral/fileexchange/24913-histogram-binwidth-optimization
                    
                    % also:
                    % h = 2?IQR?n(-1/3) N - (max-min)/h %Freedman-Diaconis
                    
                    %[pks,locs]= findpeaks(cnt);
                    [pks,locs]= findpeaks(cnt,'MINPEAKDISTANCE',10);
                    %
                    candidate_index = min(locs);
                    %
                    if cnt(1) > cnt(candidate_index)
                        candidate_index = 1;
                    end
                    %
                    lys_overpeak_ratio = 1.5;
                    if candidate_index > 15; % ?? can happen if very crowded only...
                        t = vls(candidate_index); % or, one can try min on inverted histo..
                    else 
                        t = vls(candidate_index) + lys_overpeak_ratio*(vls(candidate_index)-minval);
                    end                             
                 % histogram analysis to find threshold t - ends
                 sgm_lys = ul > t; 
                                                                    
                 sgm_lys = imopen(sgm_lys,strel('disk',1));
                 S1 = 3;
                 sgm_lys = bwareaopen(sgm_lys,S1*S1);                  
                 sgm_lys = sgm_lys & sgm_phal_exp &~ sgm_nukes;
            %%%%%%%%%%%%%%%%%% LYSOSOSMES       
           waitbar(3/5,hw); drawnow;                        
                        
            %%%%%%%%%%%%%%%%%% CELLS
            dist_from_nuc = 130; %100
            z2 = imdilate(sgm_nukes,strel('disk',dist_from_nuc,0));
            sgm_cells = sgm_phal_exp & z2; 
            
            % break cell clumps
            z2 = bwmorph(sgm_nukes,'thicken',Inf);
            sep_lines = ~z2;
            sep_lines(~sgm_cells)=0;
            sgm_cells(sep_lines)=0;
            
            % remove orphan pieces of cellular stuff..
            L_n = bwlabel(sgm_nukes);
            stats_n = regionprops(L_n,'Area','Centroid');    
            L_c = bwlabel(sgm_cells);                                    
            %
            for n = 1:numel(stats_n)
                xcn = fix(stats_n(n).Centroid(2));
                ycn = fix(stats_n(n).Centroid(1));
                cell_label = L_c(xcn,ycn);
                L_c(L_c==cell_label)=0;
            end            
            sgm_cells(L_c~=0)=0;          
            % remove orphan pieces of cellular stuff - ends            
            %%%%%%%%%%%%%%%%%% CELLS            
            waitbar(4/5,hw); drawnow;            
                                    
            sgm = cat(3,un,sgm_nukes,up,sgm_phal,ul,sgm_lys,sgm_cells);
                                            
                if send_to_Icy                
                    try
                        icyvol(:,:,1,1,1) = un;
                        icyvol(:,:,2,1,1) = sgm_nukes;
                        icyvol(:,:,3,1,1) = up;
                        icyvol(:,:,4,1,1) = sgm_phal;                        
                        icyvol(:,:,5,1,1) = ul;
                        icyvol(:,:,6,1,1) = sgm_lys;                                                
                        icyvol(:,:,7,1,1) = sgm_cells;                                                                        
                        %
                        notification = [obj.current_filename ' - Segmentation: RPE'];
                        if isempty(obj.h_Icy_segmentation_adjustment)
                            obj.h_Icy_segmentation_adjustment = icy_imshow(icyvol,notification);                    
                        else
                            icy_imshow(obj.h_Icy_segmentation_adjustment,icyvol,notification);                    
                        end
                    catch
                        errordlg('problem with Icy, - might be not running');
                    end                
                end 

                waitbar(5/5,hw); drawnow;            
                if ~isempty(hw), delete(hw), drawnow; end;
    end     
%-------------------------------------------------------------------------%          
 function [datas, captions, table_names, fig] = analyze_RPE(obj,~,~) 
     
%{
The purpose of the experiment was to identify effect of Rab38 on lysosomes in retinal pigment epithelial (RPE) cells.   
I isolated primary RPE cells from Rab38 mutant(chocolate) and control WT mice, 
cultured the cells on coverslips, 
fixed with 4%PFA and stained with anti-LAMP1(+Alexa 568-secondary), DAPI and phalloidin-647. I used confocal SF5 microscope to collect the images.
 
Lysosomes are small red punctate structures positive for LAMP1.       
     
- Is there a difference in size between cht and WT lysosomes? So it would be great to measure the size of lysosomes (diameter or area)
 
- Is there a difference in intracellular distribution of lysosomes in cht and WT cells? 
So probably the parameter is  �distance from the nucleus�? It seems that in WT lysosomes are perinuclear, and in cht lysosomes are more diffused.
 
- What is the difference between total lysosomal �mass� (possibly �ratio of LAMP1 signal to phalloidin�)
 
- Is it possible to see the difference in the nucleus size (diameter or area) or granularity between Rab38 and WT mutant? 
Nucleus is stained with the DAPI (blue). It is an indirect parameter of cell wellbeing which might be compromised when lysosomes are not digesting material properly.     
     
1. - measure granule size - try Fourier?
2. - distance from nucleus - use nuclei distance map (break confluent nukes beforehand..) use distance weighted with intensity (sum(dist_to_nuc.*U))/sum(U)    
3. - lysosomal mass - sum(UL)/sum(Unuc) or/and sum(UL)/sum(Unuc) or ratio of areas of lysossomes to nuclei or to phaloidin
4. - nuclear size srtaightforward, but nucleus granularity not easy as there is saturation       
%}
   
                sgm = obj.do_RPE_Segmentation(false);
                
                un = squeeze(sgm(:,:,1));
                sgm_nukes = squeeze(sgm(:,:,2));
                up = squeeze(sgm(:,:,3));
                sgm_phal = squeeze(sgm(:,:,4));
                ul = squeeze(sgm(:,:,5));
                sgm_lys = squeeze(sgm(:,:,6));                
                %
                sgm_cells = squeeze(sgm(:,:,7));
                                                
                % per field of view parameters
                total_phalloidin_length = sum(sgm_phal(:)); % fiber length
                total_nuclear_area = sum(sgm_nukes(:));  
                total_cell_area = sum(sgm_cells(:));  
                total_lysosomes_area = sum(sgm_lys(:));                  
                
                % quantifiers estimating the effective distance of from lysososmes to nuclei
                nuc_dist_map = bwdist(sgm_nukes);
                sample = nuc_dist_map(sgm_lys~=0);
                average_lysosome_distance_to_nucleus = mean(sample(:));
                % same thing lysosome intensity weighted
                average_lysosome_distance_to_nucleus_weighted = sum(sum(nuc_dist_map.*ul.*sgm_lys))/sum(sum(ul.*sgm_lys));
                
                %%%%% CELL BY CELL QUANTIFICATION
                 nuc_labs = bwlabel(sgm_nukes); 
                 stats_n = regionprops(nuc_labs,'Area','Centroid');                                      
                 cell_labs = bwlabel(sgm_cells);
                 stats_c = regionprops(cell_labs,'Area','Centroid','PixelList','EquivDiameter');
                 nuc_cell_lut = zeros(1,numel(stats_n));
                 for n = 1:numel(stats_n)
                    nuc_cell_lut(n) = cell_labs(fix(stats_n(n).Centroid(2)),fix(stats_n(n).Centroid(1)));
                 end
                 %
                 data = [];
                 for n = 1:numel(stats_n)
                     cur_cell_lab = nuc_cell_lut(n);
                     if 0~=cur_cell_lab
                         cell_img = (cell_labs==cur_cell_lab)>0;
                         cur_phalloidin_length = sum(sum(sgm_phal.*cell_img));
                         cur_nuclear_area = stats_n(n).Area;
                         cur_lysosomes_area = sum(sum(sgm_lys.*cell_img));
                         sample = nuc_dist_map(sgm_lys & cell_img);
                         cur_lysosome_distance_to_nucleus = mean(sample(:));
                         cur_lysosome_distance_to_nucleus_weighted = sum(sum(nuc_dist_map.*ul.*sgm_lys.*cell_img))/sum(sum(ul.*sgm_lys.*cell_img));
                         cur_cell_area = stats_c(cur_cell_lab).Area;         
                         %
                         % cell gyration radius, etc..
                         pl = stats_c(cur_cell_lab).PixelList;
                         xvals = pl(:,1);
                         yvals = pl(:,2);
                         XC = stats_c(cur_cell_lab).Centroid(1);
                         YC = stats_c(cur_cell_lab).Centroid(2);
                         cur_cell_Rg = sqrt(sum((xvals-XC).*(xvals-XC)+(yvals-YC).*(yvals-YC))/length(xvals));
                         cur_cell_effective_diameter = stats_c(cur_cell_lab).EquivDiameter;
                         %
                         record = {obj.current_filename, ...
                         n, ...
                         stats_n(n).Centroid(1), ...
                         stats_n(n).Centroid(2), ...            
                         cur_cell_lab, ...                         
                         cur_nuclear_area, ...
                         cur_phalloidin_length, ...
                         cur_lysosomes_area, ...
                         cur_lysosome_distance_to_nucleus, ...
                         cur_lysosome_distance_to_nucleus_weighted, ...
                         cur_cell_area, ...
                         cur_cell_Rg, ...
                         cur_cell_effective_diameter};                        
                         %
                         data = [data; record];                                                                            
                     end
                 end                                
                %%%%% CELL BY CELL QUANTIFICATION                
                                                                                                                                
                        icyvol(:,:,1,1,1) = un;
                        icyvol(:,:,2,1,1) = nuc_labs;
                        icyvol(:,:,3,1,1) = up;
                        icyvol(:,:,4,1,1) = sgm_phal;                        
                        icyvol(:,:,5,1,1) = ul;
                        icyvol(:,:,6,1,1) = sgm_lys; 
                        icyvol(:,:,7,1,1) = cell_labs;
                
                fig = icyvol;            
                
                %  it looks cumbersome to save cell-by-cell AND per-FOV
                %  data, especially in a batch.. 
                %  so maybe think about OR, for now via comment/uncomment
                
                datas = data;
                captions = {'filename','nuc_index','nuc_cX','nuc_cY','cell_index','nuclear_area','phalloidin_length', ...
                    'lysosomes_area','lysosome_distance_to_nucleus','lysosome_distance_to_nucleus_weighted','cur_cell_area','cur_cell_Rg','cur_cell_effective_diameter'};
                table_names = {'RPE_CELL_BY_CELL'};                                                                                                               
                                
%                  datas = {obj.current_filename,numel(stats_n),total_nuclear_area,total_cell_area,total_phalloidin_length,total_lysosomes_area,average_lysosome_distance_to_nucleus,average_lysosome_distance_to_nucleus_weighted};
%                  captions = {'filename','#nuc','total_nuclear_area','total_cell_area','total_phalloidin_length','total_lysosomes_area','average_lysosome_distance_to_nucleus','average_lysosome_distance_to_nucleus_weighted'};
%                  table_names = {'RPE_per_FOV'};                                                                                                               
 end
%-------------------------------------------------------------------------%          
function sgm = do_Sparks_Segmentation(obj,send_to_Icy,~)                     
         
     sgm = ones(size(obj.imgdata,1),size(obj.imgdata,2));
            
     if send_to_Icy                
        try
            icyvol(:,:,1,1,1) = mask;
            %
            notification = [obj.current_filename ' - Segmentation: Sparks'];
            if isempty(obj.h_Icy_segmentation_adjustment)
                obj.h_Icy_segmentation_adjustment = icy_imshow(icyvol,notification);                    
            else
                icy_imshow(obj.h_Icy_segmentation_adjustment,icyvol,notification);                    
            end
        catch
                errordlg('problem with Icy, - might be not running');
        end                
     end             
    
end
%-------------------------------------------------------------------------%          
 function [datas, captions, table_names, fig] = analyze_Sparks(obj,~,~) 
     datas = [];
     captions = [];
     table_names = 'ALYtools data';
     fig = [];     
     
     mask = obj.do_Sparks_Segmentation(false);                          
     %
     Parameters.PreProcessedData = squeeze(double(obj.imgdata));
     Parameters.RegionOfInterestMask = mask;
     % Parameters.ThresholdHigh = 420;
     %
     % calculate threshold via the Threshold parameter by 3.5 of sigma,
     % assuming symetrical distribution
                    [sX,sY,sZ,sC,sT] = size(obj.imgdata);
                    mask5 = repmat(mask,[1,1,1,1,sT]);
                    z = double(obj.imgdata(mask5~=0)); % discarding empty spaces                                                                                         
                    minval = min(z(:));
                    [cnt,vls] = hist(z(:),100);
                    %
                    [pks,locs]= findpeaks(cnt,'MINPEAKDISTANCE',10);
                    %
                    candidate_index = min(locs); % the dimmest peak - bckgrnd
                    %
                    if cnt(1) > cnt(candidate_index)
                        candidate_index = 1;
                    end
                    %
                    overpeak_ratio = 0.477; % percentage of the value 
                    t = vls(candidate_index) + overpeak_ratio*(vls(candidate_index)-minval);
% Suggested Threshold for analysis and promt for input
prompt = {'Threshold (Suggested Value):'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {num2str(t)};
t1 = inputdlg(prompt,dlg_title,num_lines,defaultans);
t = str2num(cell2mat(t1(1)));
%
     Parameters.ThresholdHigh = t;                              
     %
     res = SparkDetection(Parameters);
     %
     D = res.CurrentSparkMap;
     fig = reshape(D,[size(D,1) size(D,2) 1 1 size(D,3)]);
 end
%-------------------------------------------------------------------------%            
    function sgm = do_per_image_TCSPC_FLIM_Segmentation(obj,send_to_Icy,~)     
                    
            sgm = [];
            %
            for k=1:numel(obj.M_imgdata)
                    obj.imgdata = squeeze(sum(obj.M_imgdata{k},5));
                    sgm{k} = obj.do_PR_Segmentation(false);
                    %
                    if send_to_Icy
                            icyvol(:,:,1,1,1) = obj.imgdata;
                            icyvol(:,:,2,1,1) = sgm{k};
                            icy_imshow(icyvol,char(obj.M_filenames{k}));
                    end
            end                                    
    end            
%-------------------------------------------------------------------------%  
    function [datas, captions, table_names, fig] = analyze_per_image_TCSPC_FLIM(obj,~,~)            
                datas = [];
                captions = [];
                table_names = 'ALYtools data';
                fig = [];                
                % 
                global Tp;
                global DT;
                global t;
                global IRF;
                global fitting_mask;
                global y;
                %
                global E_avr;
                global tau_FRET_avr;
                global E;
                %
                % for double_sk2
                global E_avr_1;
                global tau_FRET_avr_1;
                global E_avr_2;
                global tau_FRET_avr_2;                
                %
                global fixed_tauD;
                
                E = linspace(1e-4,1-1e-4,1000);
                %
                % %polynomial "p"
                load('eta_to_Eavr_coefs.mat');
                
                mega = 1e6;
                nano = 1e-9;
                pico = 1e-12;
                sec = 1;
                Hz = 1;

                IRF = obj.per_image_TCSPC_FLIM_irf;
                %
                if isempty(IRF)
                    errordlg('no IRF loaded, can not contiue')
                    return;
                end
                Nbins = length(IRF);

                f = obj.per_image_TCSPC_FLIM_rep_rate*mega*Hz;
                Tp = 1/f/(pico*sec);
                DT = Tp/Nbins;
                t = DT*(0:Nbins - 1);
                %
                shift = 0;
                if 0~=obj.per_image_TCSPC_FLIM_irf_shift
                    shift = shift + obj.per_image_TCSPC_FLIM_irf_shift/DT;
                end                
                IRF = shift_curve(IRF,shift);
                %
                IRF = IRF/sum(IRF);                
                %
                IRF = IRF - obj.per_image_TCSPC_FLIM_irf_background;
                IRF(IRF<0)=0;
                %                                                                
                bckg = obj.per_image_TCSPC_FLIM_background_value;
                %
                tvb = obj.per_image_TCSPC_FLIM_tvb*obj.per_image_TCSPC_FLIM_tvb_scaling;
                
                if ~isempty(tvb) % sizes should match
                    if length(tvb) ~= length(IRF)
                        errordlg('length of tvb and IRF dont match, can not continue');
                        return;
                    end
                else
                    tvb = zeros(size(IRF));                    
                end
                %
                include_background_to_model = true;
                %
                fitting_mask = ones(size(t));
                %                
                i_min = round(obj.per_image_TCSPC_FLIM_Tmin/DT);
                i_max = round(obj.per_image_TCSPC_FLIM_Tmax/DT);   
                % NB - IF ONE TRIES TO CROP IRF THE SAME WAY - IT DOESN'T LOOK GOOD
                if 0~= i_min
                    fitting_mask(1:i_min)=0;
                end
                if 0~=i_max % TO DO
                    fitting_mask(i_max:Nbins)=0;
                end
                %
                Ndecays = numel(obj.M_imgdata);                
                % to remove left and right dead zones ..   
                FITTING_MASK = repmat(fitting_mask,[1 Ndecays]);                                                
                %
                y = sparse([]);
                %
                fig = cell(1,numel(obj.M_imgdata));
                %
                raw_decays = zeros(Nbins,Ndecays); % for export
                bckg_corrected_decays = zeros(Nbins,Ndecays); % for export
                %
                if obj.per_image_TCSPC_FLIM_nonimaging % cuvette data
                    %
                   for k=1:numel(obj.M_imgdata)
                        %
                        D = squeeze(obj.M_imgdata{k});
                        %
                        if Nbins ~= numel(D)
                            errordlg('data and IRF are not compatible, can not continue');
                            return;                    
                        end
                        % 
                        raw_decays(:,k) = D;                        
                        D = D - (bckg + tvb);
                        D(D<0)=0;
                        bckg_corrected_decays(:,k) = D;                        
                        y = [y D'];                                                
                   end   
                   %
                   numbers_of_pixels = ones(1,numel(obj.M_imgdata));
                   %
                else % imaging
                    %
                    sgm = obj.do_per_image_TCSPC_FLIM_Segmentation(false);
                    %
                    for k=1:numel(obj.M_imgdata)
                        u = squeeze(obj.M_imgdata{k});                        
                        if Nbins ~= size(u,3)
                            errordlg('data and IRF are not compatible, can not continue');
                            return;                    
                        end
                    end
                    %
                    numbers_of_pixels = zeros(1,numel(obj.M_imgdata));
                    %                                                                                
                    for k=1:numel(obj.M_imgdata)
                        %
                        u = double(squeeze(obj.M_imgdata{k})); 
                        %
                        %averaging
                        if obj.per_image_TCSPC_FLIM_averaging_sigma > 0
                            sigma = obj.per_image_TCSPC_FLIM_averaging_sigma;
                            for mm=1:size(u,3)
                                uu = squeeze(u(:,:,mm));
                                [I_avr] = gsderiv(uu,sigma,0);
                                u(:,:,mm) = I_avr;
                            end
                        end
                        %                                                                                                                       
                        mask = sgm{k};       
                        Npixs_k = sum(mask(:));
                        numbers_of_pixels(k) = Npixs_k;
                        %
                        u = double(u).*repmat(double(mask),[1 1 256]);
                        u = squeeze(sum(u,1));
                        D = squeeze(sum(u,1)); % raw total decay
                        %
                        raw_decays(:,k) = D;                        
                        D = D - Npixs_k*(bckg + tvb');
                        D(D<0)=0;
                        bckg_corrected_decays(:,k) = D;                                                
                        %
                        y = [y D];                        
                        %
                        % take care of figs                    
                        icyvol(:,:,1,1,1) = squeeze(sum(obj.M_imgdata{k},5)); % intensity
                        icyvol(:,:,2,1,1) = mask; % intensity
                        fig{k} = icyvol;
                        %
                    end
                    %
                end % if obj.per_image_TCSPC_FLIM_nonimaging % cuvette data
                %
                y(y<0)=0;                                           
                %
                % to correct the bug found by Sean
                RF = obj.per_image_TCSPC_FLIM_weights_resampling_factor;
                y_ = resample(resample(full(y),1,RF),RF,1);
                w = 1./sqrt(y_);                
                w(isinf(w))=0;
                %
                w(~FITTING_MASK) = 0;
                %
                % take care of saturation ..
                w = w.*( y < obj.per_image_TCSPC_FLIM_saturation_value ); 
                %
                options = optimset('Display','iter','DerivativeCheck','on');
                %
                D = [];
                %
                chi2_index = 7; % 2exp
                %
                switch  char(obj.per_image_TCSPC_FLIM_fit_model)
                    
                    case '1exp'
                        
                        chi2_index = 3;
                        
                        alphainit = 1000; 
                        %                            
                        lb = 0;
                        ub = Inf;
                        % 
                        tic
                            if obj.per_image_TCSPC_FLIM_conv_irf_pp_69_70
                                [alpha,c,wresid,resid_norm,y_est,Regression] = ...
                                    varpro(y',w',alphainit,Ndecays,@adaex_1exp_FLIM_pp_69_70,lb,ub,options);
                            else
                                [alpha,c,wresid,resid_norm,y_est,Regression] = ...
                                    varpro(y',w',alphainit,Ndecays,@adaex_1exp_FLIM,lb,ub,options);
                            end
                        toc
                            %
                        tau = alpha(1);                        
                        %
                        chi2 = zeros(size(c));
                        %
                        % stuff the output - starts
                        for m=1:numel(obj.M_imgdata)                            
                            record = {char(obj.M_filenames{m}), ...
                                        tau ...,
                                        chi2(m)...,
                                        c(m),numbers_of_pixels(m)};
                            D = [D; record];    
                        end
                        %
                        datas = D;
                        captions = {'filename','tau','chi2','Nphot','Npixs'};
                        table_names = {'1exp'};
                        % stuff the output - ends                        
                                            
                    case '2exp'
                        chi2_index = 7;
                        %
                        if obj.per_image_TCSPC_FLIM_fixed_tauD > 0 %fixed
                        
                            fixed_tauD = obj.per_image_TCSPC_FLIM_fixed_tauD;
                            %
                            alphainit = 500; 
                            %
                            lb = 0; % DT; % [] % :)
                            ub = Inf; % []
                            % 
                            tic
                                if obj.per_image_TCSPC_FLIM_conv_irf_pp_69_70                            
                                    [alpha,c,wresid,resid_norm,y_est,Regression] = ...
                                        varpro(y',w',alphainit,Ndecays*2,@adaex_2exp_FLIMFRET_fixed_tauD_pp_69_70,lb,ub,options);
                                else
                                    [alpha,c,wresid,resid_norm,y_est,Regression] = ...
                                        varpro(y',w',alphainit,Ndecays*2,@adaex_2exp_FLIMFRET_fixed_tauD,lb,ub,options);
                                end
                            toc
                            %
                            tauD_2exp = fixed_tauD;
                            tauFRET_2exp = alpha(1);
                        
                        else % non-fixed
                            
                            alphainit = [3000; 500]; 
                            %
                            % lb = [1000   DT]'; % [] % :)
                            lb = [100   0]'; % [] % :)
                            ub = [Inf   Inf]'; % []
                            Nnonlins  = length(alphainit);
                            % 
                            tic
                                if obj.per_image_TCSPC_FLIM_conv_irf_pp_69_70                                                        
                                    [alpha,c,wresid,resid_norm,y_est,Regression] = ...
                                        varpro(y',w',alphainit,Ndecays*2,@adaex_2exp_FLIMFRET_pp_69_70,lb,ub,options);
                                else
                                    [alpha,c,wresid,resid_norm,y_est,Regression] = ...
                                        varpro(y',w',alphainit,Ndecays*2,@adaex_2exp_FLIMFRET,lb,ub,options);                                
                                end
                            toc
                            %
                            tauD_2exp = alpha(1);
                            tauFRET_2exp = alpha(2);
                        end
                        
                        n = length(c);
                        c1=c(1:2:n);
                        c2=c(2:2:n);
                        Nphot = c1+c2;
                        %
                        mf_2exp = 1 - (c1/tauD_2exp)./(c1/tauD_2exp + c2/tauFRET_2exp);
                        FRET_efficiency_2exp = 1 - tauFRET_2exp/tauD_2exp;
                        eta_2exp = (1/FRET_efficiency_2exp-1)^(1/6);
                        %
                        chi2 = zeros(size(mf_2exp));
                        %
                        % stuff the output - starts
                        for m=1:numel(obj.M_imgdata)                            
                            record = {char(obj.M_filenames{m}), ...
                                        tauD_2exp ...,
                                        tauFRET_2exp ...,
                                        mf_2exp(m) ..., % FRET molar fraction
                                        FRET_efficiency_2exp ...,
                                        eta_2exp, ...
                                        chi2(m),Nphot(m),numbers_of_pixels(m)};
                            D = [D; record];                             
                        end         
                        %
                        datas = D;
                        captions = {'filename','tau_D','tau_FRET','beta_FRET','FRET_efficiency','r_div_R0','chi2','Nphot','Npixs'};
                        table_names = {'2exp'};
                        % stuff the output - ends                        
                                                
                    case 'sk2'
                        chi2_index = 7;
                        
                        if obj.per_image_TCSPC_FLIM_fixed_tauD > 0 %fixed

                            fixed_tauD = obj.per_image_TCSPC_FLIM_fixed_tauD;
                            
                            alphainit = 1.2;  
                            lb = 0.05;
                            ub = 2;
                            % 
                            tic
                                if obj.per_image_TCSPC_FLIM_conv_irf_pp_69_70
                                    [alpha,c,wresid,resid_norm,y_est,Regression] = ... 
                                        varpro(y',w',alphainit,Ndecays*2,@adaex_sk2_FLIMFRET_fixed_tauD_pp_69_70,lb,ub,options);                                
                                else
                                    [alpha,c,wresid,resid_norm,y_est,Regression] = ... 
                                        varpro(y',w',alphainit,Ndecays*2,@adaex_sk2_FLIMFRET_fixed_tauD,lb,ub,options);                                    
                                end
                            toc
                            %
                            tauD_sk2 = fixed_tauD;
                            eta_f = alpha(1);
                            
                        else % non-fixed
                            %
                            alphainit = [5000; 1.2];  
                            % lb = [DT 0.05]';
                            lb = [0 0.05]';
                            ub = [Inf 2]';
                            % 
                            tic
                                if obj.per_image_TCSPC_FLIM_conv_irf_pp_69_70
                                    [alpha,c,wresid,resid_norm,y_est,Regression] = ... 
                                        varpro(y',w',alphainit,Ndecays*2,@adaex_sk2_FLIMFRET_pp_69_70,lb,ub,options);
                                else
                                    [alpha,c,wresid,resid_norm,y_est,Regression] = ...                                     
                                        varpro(y',w',alphainit,Ndecays*2,@adaex_sk2_FLIMFRET,lb,ub,options);                                        
                                end
                            toc
                            %
                            tauD_sk2 = alpha(1);
                            eta_f = alpha(2);
   
                        end
                        %
                        n = length(c);  
                        c1=c(1:2:n);
                        c2=c(2:2:n);
                        Nphot = c1+c2;
                            
                        mf_sk2 = (c1/tau_FRET_avr)./(c1/tau_FRET_avr + c2/tauD_sk2);
                        Eavr_f = polyval(p,eta_f);
                        %
                        chi2 = zeros(size(mf_sk2));
                        %                        
                        % stuff the output - starts
                        for m=1:numel(obj.M_imgdata)                            
                            record = {char(obj.M_filenames{m}), ...
                                        tauD_sk2 ...,
                                        tau_FRET_avr ...,
                                        mf_sk2(m) ..., % FRET molar fraction
                                        Eavr_f ...,
                                        eta_f, ...
                                        chi2(m),Nphot(m),numbers_of_pixels(m)};
                            D = [D; record];                             
                        end         
                        %
                        datas = D;
                        captions = {'filename','tau_D','tau_FRET','beta_FRET','FRET_efficiency','r_div_R0','chi2','Nphot','Npixs'};
                        table_names = {'sk2'};
                        % stuff the output - ends                        
                                                
                    case 'sk2 only'
                        chi2_index = 6;
                        
                        if obj.per_image_TCSPC_FLIM_fixed_tauD > 0 %fixed
                            
                            fixed_tauD = obj.per_image_TCSPC_FLIM_fixed_tauD;
                            
                            alphainit = 1.2;  
                            lb = 0.05;
                            ub = 2;
                            % 
                            tic
                                if obj.per_image_TCSPC_FLIM_conv_irf_pp_69_70
                                    [alpha,c,wresid,resid_norm,y_est,Regression] = ... 
                                        varpro(y',w',alphainit,Ndecays,@adaex_sk2_only_FLIMFRET_fixed_tauD_pp_69_70,lb,ub,options);                               
                                else
                                    [alpha,c,wresid,resid_norm,y_est,Regression] = ... 
                                        varpro(y',w',alphainit,Ndecays,@adaex_sk2_only_FLIMFRET_fixed_tauD,lb,ub,options);                                                               
                                end
                            toc
                            %
                            tauD_sk2 = fixed_tauD;
                            % no betas
                            eta_f = alpha(1);
                            
                        else % non-fixed
                            
                            alphainit = [5000; 1.2];  
                            % lb = [DT 0.05]';
                            lb = [0 0.05]';
                            ub = [Inf 2]';
                            % 
                            tic
                                if obj.per_image_TCSPC_FLIM_conv_irf_pp_69_70
                                    [alpha,c,wresid,resid_norm,y_est,Regression] = ... 
                                        varpro(y',w',alphainit,Ndecays,@adaex_sk2_only_FLIMFRET_pp_69_70,lb,ub,options);                                    
                                else
                                    [alpha,c,wresid,resid_norm,y_est,Regression] = ... 
                                        varpro(y',w',alphainit,Ndecays,@adaex_sk2_only_FLIMFRET,lb,ub,options);                                                                        
                                end
                            toc
                            %
                            tauD_sk2 = alpha(1);
                            % no betas
                            eta_f = alpha(2);
                        end
                        %
                        Eavr_f = polyval(p,eta_f);
                        %
                        chi2 = zeros(1,numel(obj.M_imgdata));
                        %
                        Nphot = c;
                        %
                        % stuff the output - starts
                        for m=1:numel(obj.M_imgdata)                            
                            record = {char(obj.M_filenames{m}), ...
                                        tauD_sk2 ...,
                                        tau_FRET_avr ...,
                                        Eavr_f ...,
                                        eta_f, ...
                                        chi2(m),Nphot(m),numbers_of_pixels(m)};
                            D = [D; record];                             
                        end         
                        %
                        datas = D;
                        captions = {'filename','tau_D','tau_FRET','FRET_efficiency','r_div_R0','chi2','Nphot','Npixs'};
                        table_names = {'sk2 only'};
                        % stuff the output - ends                                                
                            
                    case 'double sk2'
                        chi2_index = 10;
                        
                         if obj.per_image_TCSPC_FLIM_fixed_tauD > 0 %fixed
                             
                            fixed_tauD = obj.per_image_TCSPC_FLIM_fixed_tauD;
                            
                            alphainit = [0.8; 1.2];  
                            lb = [0.05 0.05]';
                            ub = [2 2]';
                            %
                            tic
                                if obj.per_image_TCSPC_FLIM_conv_irf_pp_69_70
                                    [alpha,c,wresid,resid_norm,y_est,Regression] = ...
                                         varpro(y',w',alphainit,Ndecays*2,@adaex_sk2_sk2_FLIMFRET_fixed_tauD_pp_69_70,lb,ub,options);
                                else
                                    [alpha,c,wresid,resid_norm,y_est,Regression] = ...
                                         varpro(y',w',alphainit,Ndecays*2,@adaex_sk2_sk2_FLIMFRET_fixed_tauD,lb,ub,options);
                                end
                            toc
                            %
                            tauD_22 = fixed_tauD;
                            eta_f_1 = alpha(1); 
                            eta_f_2 = alpha(2);
                            
                        else % non-fixed                       
                            %
                            alphainit = [2700; 0.8; 1.2];  
                            lb = [100 0.05 0.05]';
                            ub = [Inf 2 2]';
                            %
                            tic
                                if obj.per_image_TCSPC_FLIM_conv_irf_pp_69_70
                                    [alpha,c,wresid,resid_norm,y_est,Regression] = ...
                                         varpro(y',w',alphainit,Ndecays*2,@adaex_sk2_sk2_FLIMFRET_pp_69_70,lb,ub,options);                                    
                                else
                                    [alpha,c,wresid,resid_norm,y_est,Regression] = ...
                                         varpro(y',w',alphainit,Ndecays*2,@adaex_sk2_sk2_FLIMFRET,lb,ub,options);                                                                        
                                end
                            toc
                            %
                            tauD_22 = alpha(1);
                            eta_f_1 = alpha(2); 
                            eta_f_2 = alpha(3);                           
                         end
                        %
                        n = length(c);
                        c1=c(1:2:n);
                        c2=c(2:2:n);
                        pf = c1./(c1+c2);
                        Nphot = c1+c2;
                        %
                        mf_sk2_22 = (c1/tau_FRET_avr_1)./(c1/tau_FRET_avr_1 + c2/tau_FRET_avr_2);
      
                        Eavr_f_1 = polyval(p,eta_f_1);
                        Eavr_f_2 = polyval(p,eta_f_2);                                

                        chi2 = zeros(size(mf_sk2_22));                       
                        
                        % stuff the output - starts
                        for m=1:numel(obj.M_imgdata)                            
                            record = {char(obj.M_filenames{m}), ...
                                        tauD_22 ...,
                                        tau_FRET_avr_1 ...,
                                        Eavr_f_1 ...,
                                        eta_f_1 ...,
                                        tau_FRET_avr_2 ...,
                                        Eavr_f_2 ...,
                                        eta_f_2 ...,                                        
                                        mf_sk2_22(m), ...
                                        chi2(m),Nphot(m),numbers_of_pixels(m)}; % FRET molar fraction of the first term
                            D = [D; record];                             
                        end         
                        %
                        datas = D;
                        captions = {'filename','tau_D','tau_FRET_1','FRET_efficiency_1','r_div_R0_1', ...
                                                       'tau_FRET_2','FRET_efficiency_2','r_div_R0_2', ...
                                                       'beta_1','chi2','Nphot','Npixs'};
                        table_names = {'double sk2'};
                        % stuff the output - ends                                                                                                
                end
                %                
                % chi squared                 
                chi2s = zeros(Ndecays,1);                
                n_activebins = numel(fitting_mask(fitting_mask~=0));
                r = y-y_est';
                norm_res = r./sqrt(y_est');                                    
                for k=1:Ndecays
                        RANGE = ((k-1)*Nbins+1):(k*Nbins);
                        s = norm_res(RANGE).*fitting_mask;
                        chi2s(k) = sum(s.*s)/n_activebins;
                end                
                datas(:,chi2_index) = num2cell(chi2s);
                % chi squared
                %
                % averages of parameters
                ncols = size(datas,2);
                avrrec = [{'average'} num2cell(mean(cell2mat(datas(:,2:ncols)),1))];
                datas = [datas; avrrec];
                %                
                res.t = t;
                res.y = y;
                res.y_est = y_est;
                res.fitting_mask = fitting_mask;
                res.data_controller = obj;
                res.Ndecays = Ndecays;
                res.Nbins = Nbins;
                res.datas = datas;
                res.captions = captions;
                res.table_names = table_names;
                res.IRF = IRF;
                res.raw_decays = raw_decays;
                res.bckg_corrected_decays = bckg_corrected_decays;
                
                per_image_TCSPC_FLIM_results_panel(res);
    end        
%-------------------------------------------------------------------------%          
function load_irf(obj,~,~)    
             obj.per_image_TCSPC_FLIM_irf_filename = [];
             [fname, fpath] = uigetfile({'*.sdt;*.csv;*.xls;*.txt'},'Select IRF file',obj.DefaultDirectory);
             if fpath == 0, return, end;
             filespec = fullfile(fpath,fname);
             %
             extension = fname(length(fname)-2:length(fname));
             %
             if strcmp(extension,'sdt')
                 try
                    [~,~,I] = bfopen_v(filespec);
                    [sX,sY,sZ,sC,sT] = size(I);
                    I = reshape(I,[sX,sY,sZ,sT,sC]);
                    I = double(squeeze(I(:,:,1,:,1)));
                    irf = squeeze(sum(sum(I,1),2));
                    irf = irf/sum(irf(:));
                    % figure;semilogy(irf);
                    obj.per_image_TCSPC_FLIM_irf = irf;
                    obj.per_image_TCSPC_FLIM_irf_filename = filespec;
                 catch
                    errordlg('Error while trying to load sdt IRF');
                 end
             elseif strcmp(extension,'csv') || strcmp(extension,'xls')
                 try
                    [NUM,TXT,RAW]=xlsread(filespec);
                    irf = NUM(:,2);
                    irf = irf/sum(irf(:));                    
                    obj.per_image_TCSPC_FLIM_irf = irf;
                    obj.per_image_TCSPC_FLIM_irf_filename = filespec;
                 catch
                    errordlg('Error while trying to load csv IRF');
                 end
             elseif strcmp(extension,'txt')
                 try
                   data = importdata(filespec);
                    if isnumeric(data)
                        data_v = data(:,2);
                    else
                        nrows = size(data.data,1);
                        if nrows - 256 < 20 % hack
                            L = 256;
                        end
                         offset = nrows - L + 1 ;
                         data_v = data.data(offset:nrows,2);
                    end                                                            
                    irf = data_v;
                    irf = irf/sum(irf(:));                    
                    obj.per_image_TCSPC_FLIM_irf = irf;
                    obj.per_image_TCSPC_FLIM_irf_filename = filespec;
                 catch
                    errordlg('Error while trying to load txt IRF');
                 end    
             end                 
end
%-------------------------------------------------------------------------%          
function load_tvb(obj,~,~)    
             obj.per_image_TCSPC_FLIM_tvb_filename = [];
             [fname, fpath] = uigetfile({'*.xls;*.xlsx;*.csv;'},'Select tvb file',obj.DefaultDirectory);
             if fpath == 0, return, end;
             filespec = fullfile(fpath,fname);
             try
                [NUM,TXT,RAW] = xlsread(filespec);                                
                obj.per_image_TCSPC_FLIM_tvb = NUM(:,2);
                obj.per_image_TCSPC_FLIM_tvb_filename = filespec;                
             catch
                errordlg('Error while trying to load tvb');
             end
end
%-------------------------------------------------------------------------%          
function phasors = per_image_TCSPC_FLIM_get_phasors(obj,~,~)
    
                phasors = [];
    
                mega = 1e6;
                pico = 1e-12;
                sec = 1;
                Hz = 1;

                IRF = obj.per_image_TCSPC_FLIM_irf;
                %
                if isempty(IRF)
                    errordlg('no IRF loaded, can not contiue')
                    return;
                end
                Nbins = length(IRF);

                f = obj.per_image_TCSPC_FLIM_rep_rate*mega*Hz;
                Tp = 1/f/(pico*sec);
                DT = Tp/Nbins;
                t = DT*(0:Nbins - 1);
                %
                shift = 0;
                if 0~=obj.per_image_TCSPC_FLIM_irf_shift
                    shift = shift + obj.per_image_TCSPC_FLIM_irf_shift/DT;
                end                
                IRF = shift_curve(IRF,shift);
                %
                IRF = IRF/sum(IRF);                
                %
                IRF = IRF - obj.per_image_TCSPC_FLIM_irf_background;
                IRF(IRF<0)=0;
                IRF = IRF/sum(IRF); 
                %                                                                
                bckg = obj.per_image_TCSPC_FLIM_background_value;
                %
                tvb = obj.per_image_TCSPC_FLIM_tvb*obj.per_image_TCSPC_FLIM_tvb_scaling;
                
                if ~isempty(tvb) % sizes should match
                    if length(tvb) ~= length(IRF)
                        errordlg('length of tvb and IRF dont match, can not continue');
                        return;
                    end
                else
                    tvb = zeros(size(IRF));                    
                end
                %
                fitting_mask = ones(size(t));
                %                
                i_min = round(obj.per_image_TCSPC_FLIM_Tmin/DT);
                i_max = round(obj.per_image_TCSPC_FLIM_Tmax/DT);   
                % NB - IF ONE TRIES TO CROP IRF THE SAME WAY - IT DOESN'T LOOK GOOD
                if 0~= i_min
                    fitting_mask(1:i_min)=0;
                end
                if 0~=i_max % TO DO
                    fitting_mask(i_max:Nbins)=0;
                end
                    %reference frequency
                    W = 2*pi/Tp; 
                    %
                    % G,S phasor for the tau_R reference IRF
                    G_irf = sum(IRF.*cos(W*t).*fitting_mask)/sum(IRF.*fitting_mask);
                    S_irf = sum(IRF.*sin(W*t).*fitting_mask)/sum(IRF.*fitting_mask);
                    %                                                                                                                
                if obj.per_image_TCSPC_FLIM_nonimaging % cuvette data
                    return;
                else % imaging
                    %
                    if ~isempty(obj.M_sgm)
                        sgm = obj.M_sgm;
                    else
                        sgm = obj.do_per_image_TCSPC_FLIM_Segmentation(false);
                    end
                    %
                    for k=1:numel(obj.M_imgdata)
                        u = squeeze(obj.M_imgdata{k});                        
                        if Nbins ~= size(u,3)
                            errordlg('data and IRF are not compatible, can not continue');
                            return;                    
                        end
                    end
                    %
                    hw = waitbar(0,'Collecting phasors, please wait');                    
                    for k=1:numel(obj.M_imgdata)
                        %
                        if ~isempty(hw), waitbar(k/numel(obj.M_imgdata),hw); drawnow, end;                            
                        %
                        u = double(squeeze(obj.M_imgdata{k}));
                                               
                        %averaging
                        if obj.per_image_TCSPC_FLIM_averaging_sigma > 0
                            sigma = obj.per_image_TCSPC_FLIM_averaging_sigma;
                            for m=1:size(u,3)
                                I = squeeze(u(:,:,m));
                                [I_avr] = gsderiv(I,sigma,0);
                                u(:,:,m) = I_avr;
                            end
                        end
                        %                       
                        mask = sgm{k};       
                        %
                        phasors_k = zeros(sum(mask(:)),2);
                        index = 0;
                        for xx=1:size(u,1)
                            for yy=1:size(u,2)
                                if 0 ~= mask(xx,yy)
                                    D = squeeze(u(xx,yy,:)); %decay                                    
                                    D = D - (bckg + tvb);
                                    D(D<0)=0;
                                    %
                                    I_tot_m = sum(D'.*fitting_mask);
                                    G_m = sum(D'.*cos(W*t).*fitting_mask)/I_tot_m;
                                    S_m = sum(D'.*sin(W*t).*fitting_mask)/I_tot_m;
                                    %
                                    % irf compensation
                                    tau_R = obj.per_image_TCSPC_FLIM_Ref_lifetime;
                                    L_I = (G_m - 1i*S_m)/(G_irf - 1i*S_irf)/(1+1i*W*tau_R);
                                    G =   real( L_I );
                                    S = - imag( L_I );
                                    index = index+1;
                                    phasors_k(index,:)=[G S];
                                end
                            end
                        end
                        phasors = [phasors; phasors_k];
                    end
                    if ~isempty(hw), delete(hw), drawnow; end;
                end % if obj.per_image_TCSPC_FLIM_nonimaging % cuvette data        
end
%-------------------------------------------------------------------------%
    function [datas, captions, table_names, fig] = analyze_per_image_TCSPC_FLIM_PHASOR(obj,~,~)            
                datas = [];
                captions = [];
                table_names = 'ALYtools data';
                fig = [];
                %
                mega = 1e6;
                pico = 1e-12;
                sec = 1;
                Hz = 1;
                %    
                f = obj.per_image_TCSPC_FLIM_rep_rate*mega*Hz;
                Tp = 1/f/(pico*sec);      
                %
                phasors = obj.per_image_TCSPC_FLIM_get_phasors;
                %
                res{1}=phasors;
                res{2}=Tp;
                per_image_TCSPC_FLIM_PHASOR_results_panel(res);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XYF = do_AI_Powered_2D_SMLM_Reconstruction_Segmentation(obj,send_to_Icy,~) 

            d = obj.AI_Powered_2D_SMLM_Reconstruction_vicinity_half_width; % = 5;
            s = obj.AI_Powered_2D_SMLM_Reconstruction_extraction_scale; % = 2; %non-super res pixels            
            K = obj.AI_Powered_2D_SMLM_Reconstruction_extraction_scale_ratio; % = 3.5; % ditto
            t = obj.AI_Powered_2D_SMLM_Reconstruction_extraction_threshold; %; % std factor?
            method = obj.AI_Powered_2D_SMLM_Reconstruction_extraction_method; % = 'Linear_TopHat';

            [sX,sY,~,~,sT] = size(obj.imgdata);
            XYF = [];
            icyvol = [];

            % maybe used for harvesting
            ycoor = repmat(1:sX,[sY 1]);
            xcoor = repmat((1:sY)',[1 sX]);
                        
            r = ceil(s/2+1);

                   n_frames = 10;
                    
                   frame_data = cell(n_frames,1);         
                   %parfor k=1:n_frames % only first 10 frames
                   for k=1:n_frames % only first 10 frames

                        frame = single(squeeze(obj.imgdata(:,:,1,1,k)));    
                    % the block BELOW should be a copy of the corresponding one from "Analysis" proc for consistency
                                                            
                        switch method
                            case 'Linear_TopHat'
                                z = linear_tophat(frame,s,K);                
                            case 'DOG'
                                g1=gsderiv(frame,s,0);
                                g2=gsderiv(frame,s*K,0);
                                z=(g1-g2);
                            case 'Primitive_Linear_Tophat'
                                z = conv2(frame,ones(s)/(s*s),'same');
                                z = (frame-z);
                                    z(1:r,:)=0;
                                    z(sX-r:sX,:)=0;
                                    z(:,1:r)=0;
                                    z(:,sY-r:sY)=0;                      
                        end
                        %       
                        z = z.*(z>0);
                        sample = z(z~=0);
                        thresh = t*std(sample(:)); 
                        z(z<thresh)=0;
                        %
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% harvesting
                            z = z > imdilate(z, [1 1 1; 1 0 1; 1 1 1]); 
                            z = z > 0;                                    
                            [x,y]=find(z==1);
                            f = k*ones(length(x),1);
                            frame_data{k} = [x y f];
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% harvesting
%                         L = bwlabel(z);
%                         xw = xcoor.*z;
%                         yw = ycoor.*z;
%                         ML = max(L(:));
%                         XYF_k = zeros(ML,3);
%                         for m=1:ML
%                             pixs_m = z(L==m);
%                             denom = sum(pixs_m(:));
%                             xw_m = xw(L==m);
%                             yw_m = yw(L==m);
%                             x = sum(xw_m(:))/denom;
%                             y = sum(yw_m(:))/denom;
%                             XYF_k(m,:) = round([x y k]); % pity;        
%                         end
%                         frame_data{k} = XYF_k; 
                   end 
                   XYF = [];
                   %parfor k=1:n_frames 
                   for k=1:n_frames 
                       XYF_k = frame_data{k};
                       XYF = [XYF; XYF_k];
                   end
                    %
                    indi = XYF(:,1)-d>=1 & XYF(:,2)-d>=1 & ... 
                           XYF(:,1)+d<=sX & XYF(:,2)+d<=sY;
                    XYF = XYF(indi,:);
                    % the block ABOVE should be a copy of the corresponding one from "Analysis" proc for consistency
                    %  
                    icyvol = zeros(sX,sY,2,1,n_frames,'single');
                    icyvol(:,:,1,1,1:n_frames) = single(obj.imgdata(:,:,1,1,1:n_frames));
                    for k=1:size(XYF,1)
                        cur_frame = XYF(k,3);
                        if cur_frame <= n_frames
                            icyvol(XYF(k,1),XYF(k,2),2,1,cur_frame) = 255;
                        end
                    end
                                                
                if send_to_Icy                
                    try
                        notification = [obj.current_filename ' - Segmentation: AI_Powered_2D_SMLM_Reconstruction'];
                        if isempty(obj.h_Icy_segmentation_adjustment)
                            obj.h_Icy_segmentation_adjustment = icy_imshow(icyvol,notification);                    
                        else
                            icy_imshow(obj.h_Icy_segmentation_adjustment,icyvol,notification);                    
                        end
                    catch
                        errordlg('problem with Icy, - might be not running');
                    end 
                end
end
%-------------------------------------------------------------------------%
function [datas, captions, table_names, fig] = analyze_AI_Powered_2D_SMLM_Reconstruction(obj,~,~) 

     t_start = tic;

disp('analyze_AI_Powered_2D_SMLM_Reconstruction - extraction started!');
     
     datas = [];
     captions = [];
     table_names = 'default';
     fig = [];
     
     [sX,sY,~,~,sT] = size(obj.imgdata);
     
     path_to_network = obj.AI_Powered_2D_SMLM_Reconstruction_network; % path to file
     upscale_fac = obj.AI_Powered_2D_SMLM_Reconstruction_upscale_factor;
     d = obj.AI_Powered_2D_SMLM_Reconstruction_vicinity_half_width;
     pix_size = obj.AI_Powered_2D_SMLM_Reconstruction_pixel_size; % nm/pixel
     R = obj.AI_Powered_2D_SMLM_Reconstruction_max_distance_to_spurious_pixl/pix_size*upscale_fac; % expressed in super res pixs
     block_size = obj.AI_Powered_2D_SMLM_Reconstruction_time_dependent_block_size; % 200 frames 
     visualisation_method = obj.AI_Powered_2D_SMLM_Reconstruction_image_formation_method; % = 'ASH'
     min_sigma = obj.AI_Powered_2D_SMLM_Reconstruction_min_sigma;
     max_sigma = obj.AI_Powered_2D_SMLM_Reconstruction_max_sigma;
      
     try
        load(path_to_network); % "net" object is there
     catch % try some defaults
        path_to_network = [pwd filesep 'TrainedNetworks' filesep 'trained_net_lambda_340_NA_1_49_nmppix_106_d_5_N_400000_sigma.mat'];
        if ~isfile(path_to_network)
            path_to_network = [pwd filesep 'AI_Powered_2D_SMLM_Reconstruction_default_trained_net.mat'];
        end
        load(path_to_network);
     end
     %
     if ~exist('net','var'), return, end
     
            s = obj.AI_Powered_2D_SMLM_Reconstruction_extraction_scale; % = 2; %non-super res pixels            
            K = obj.AI_Powered_2D_SMLM_Reconstruction_extraction_scale_ratio; % = 3.5; % ditto
            t = obj.AI_Powered_2D_SMLM_Reconstruction_extraction_threshold; %; % std factor
            method = obj.AI_Powered_2D_SMLM_Reconstruction_extraction_method; % = 'Linear_TopHat';
            
            r = ceil(s/2+1);
                        
            XYF = [];

            % maybe used for harvesting
            ycoor = repmat(1:sX,[sY 1]);
            xcoor = repmat((1:sY)',[1 sX]);
            
                temp_img = obj.imgdata; % ?!!!!
                % obj.imgdata = []; % remove from memory if image size is very large
                   
                   % original - total intensity image   
                   total_intensity = zeros(sX,sY,'single');
                   parfor k=1:sT
                   %for k=1:sT
                       total_intensity = total_intensity + single(squeeze(temp_img(:,:,1,1,k)));
                   end
                   total_intensity_UPS = imresize(total_intensity,upscale_fac,'box');            
                                                         
                   n_frames = sT;                
                   frame_data = cell(n_frames,1);         
                   parfor k=1:n_frames 
                   %for k=1:n_frames

                        frame = single(squeeze(temp_img(:,:,1,1,k)));
                                                                        
                        % the block BELOW should be a copy of the corresponding one from "Analysis" proc for consistency
                        switch method
                            case 'Linear_TopHat'
                                z = linear_tophat(frame,s,K);                
                            case 'DOG'
                                g1=gsderiv(frame,s,0);
                                g2=gsderiv(frame,s*K,0);
                                z=(g1-g2);
                            case 'Primitive_Linear_Tophat'
                                z = conv2(frame,ones(s)/(s*s),'same');
                                z = (frame-z);
                                    z(1:r,:)=0;
                                    z(sX-r:sX,:)=0;
                                    z(:,1:r)=0;
                                    z(:,sY-r:sY)=0;                      
                        end
                        %       
                        z = z.*(z>0);
                        sample = z(z~=0);
                        thresh = t*std(sample(:)); 
                        z(z<thresh)=0;
                        %
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% harvesting
                            z = z > imdilate(z, [1 1 1; 1 0 1; 1 1 1]); 
                            z = z > 0;                                    
                            [x,y]=find(z==1);
                            f = k*ones(length(x),1);
                            frame_data{k} = [x y f];
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% harvesting                        
%                         L = bwlabel(z);
%                         xw = xcoor.*z;
%                         yw = ycoor.*z;
%                         ML = max(L(:));
%                         XYF_k = zeros(ML,3);
%                         for m=1:ML
%                             pixs_m = z(L==m);
%                             denom = sum(pixs_m(:));
%                             xw_m = xw(L==m);
%                             yw_m = yw(L==m);
%                             x = sum(xw_m(:))/denom;
%                             y = sum(yw_m(:))/denom;
%                             XYF_k(m,:) = round([x y k]); % pity;        
%                         end
%                         frame_data{k} = XYF_k; 
                   end
                   XYF = [];
                   parfor k=1:n_frames 
                   %for k=1:n_frames 
                       XYF_k = frame_data{k};
                       XYF = [XYF; XYF_k];
                   end
                    %
                     indi = XYF(:,1)-d>=1 & XYF(:,2)-d>=1 & ... 
                            XYF(:,1)+d<=sX & XYF(:,2)+d<=sY;

                    XYF = XYF(indi,:);
                    % the block ABOVE should be a copy of the corresponding one from "Analysis" proc for consistency

disp(['extracted ' num2str(size(XYF,1)) ' emitters, time = ' num2str(toc(t_start)/60)]);
                
     VIC = zeros(2*d+1,2*d+1,size(XYF,1));
     x = XYF(:,1);
     y = XYF(:,2);
     f = XYF(:,3);
     x1 = x-d;
     x2 = x+d;
     y1 = y-d;
     y2 = y+d;
     for k=1:size(XYF,1)
        VIC(:,:,k) = single(squeeze(temp_img(x1(k):x2(k),y1(k):y2(k),1,1,f(k))));
     end
     
     clear('temp_img');
     
disp(['vicinities extracted, time = ' num2str(toc(t_start)/60)]);
     
     %%%%%%%%%%%%%%%%%% normalize image data
     VIC_norm = zeros(size(VIC));
     parfor k=1:size(VIC,3)
     %for k=1:size(VIC,3)         
        z = VIC(:,:,k);
        z = map(z,0,1);
        z = (z - mean(z(:)))/std(z(:));
        VIC_norm(:,:,k) = z;
        %disp([' normalizing vicinity images.. ' num2str(k)]);
     end
     %%%%%%%%%%%%%%%%%% normalize image data
disp(['vicinities normalized, time = ' num2str(toc(t_start)/60)]);     

     VIC_norm = reshape(VIC_norm,2*d+1,2*d+1,1,size(VIC, 3));
     
     prediction = predict(net,VIC_norm);
     if 2 == size(prediction,2)         
        dx_dy = prediction;
        XY = XYF(:,1:2) + dx_dy; 
        XY_nm = (XY - 0.5)*pix_size;
        F = XYF(:,3);
     else % use sigma to filter
        dx_dy = prediction(:,1:2);
        sigma = prediction(:,3);
        indi = sigma>min_sigma & sigma<max_sigma;
        XY = XYF(:,1:2) + dx_dy;
        F = XYF(:,3);                
        XY = XY(indi,:);
        XY_nm = (XY - 0.5)*pix_size;        
        F = F(indi);
        sigma = sigma(indi);
     end

        % and this is real thing in super res pixels
        % ATTENTION - XY_nm IS THE RESULT - will be used later        
        XY_UPS = round(XY_nm/pix_size*upscale_fac); % stands for UPSCALED
        F_UPS = F;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [XY_UPS,F_UPS] = AI_Powered_2D_SMLM_Reconstruction_merge_localisations_in_frames(XY_UPS,F_UPS,block_size); % block-wise "unique"
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['merged localisations in frames, time = ' num2str(toc(t_start)/60) ' #localisations = ' num2str(size(XY_UPS,1))]);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % R = 50/pix_size*upscale_fac; % R = 80 nm 
        % NB what it means - if single localisation doesn't have any localsations within R, it is removed
        % all distances are in super res pixels for this proc
        [XY_UPS,F_UPS] = AI_Powered_2D_SMLM_Reconstruction_remove_spurious_localisations(upscale_fac*sX,upscale_fac*sY,XY_UPS,F_UPS,R);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['spurious localisations removed, time = ' num2str(toc(t_start)/60) ' #localisations = ' num2str(size(XY_UPS,1))]);             

        % XY_UPS is expressed in super-res pixels
        DX_DY_DRIFT = AI_Powered_2D_SMLM_Reconstruction_calculate_XY_drift_increments(upscale_fac*sX,upscale_fac*sY,XY_UPS,F_UPS,sT,upscale_fac); 

        % loss of precision..
        XY_UPS = AI_Powered_2D_SMLM_Reconstruction_implement_drift_correction(XY_UPS,F_UPS,DX_DY_DRIFT); %will contain fractional part
        XY_UPS = round(XY_UPS);
disp('..drift corrected..');                     
                       
disp(['total execution time = ' num2str(toc(t_start)/60) ' min, #localisations = ' num2str(size(XY_UPS,1))]);

            % constructing output super res image            
                switch visualisation_method
                    case 'ASH'
                        scene_AI = AI_Powered_2D_SMLM_Reconstruction_ASH_2d(upscale_fac*sX,upscale_fac*sY,XY_UPS,max(3,round(upscale_fac/2)));
                    case 'Smoothed Histogram'
                        scene_AI = zeros(upscale_fac*sX,upscale_fac*sY,'single');
                        for k=1:size(XY_UPS,1)
                              x = round(XY_UPS(k,1)); %no need to divide by pix_size!
                              y = round(XY_UPS(k,2)); 
                            if x>=1 && x<=sX*upscale_fac && y>=1 && y<=sY*upscale_fac
                                scene_AI(x,y) = scene_AI(x,y) + 1;
                            end
                        end
                        s = 1;
                        scene_AI = gsderiv(scene_AI,s,0);                                                                        
                end
                        
            non_superres_features = zeros(size(total_intensity));
            for k=1:size(XYF,1)
                non_superres_features(XYF(k,1),XYF(k,2)) = non_superres_features(XYF(k,1),XYF(k,2))+1;
            end              
            non_superres_features = imresize(non_superres_features,upscale_fac,'box');

            fig = zeros(upscale_fac*sX,upscale_fac*sY,3,1,1,'single');
            fig(:,:,1,1,1) = total_intensity_UPS;
            fig(:,:,2,1,1) = scene_AI;
            fig(:,:,3,1,1) = non_superres_features;
        
            if ~exist('sigma','var')
                datas = [XY_UPS/upscale_fac*pix_size F_UPS];
                datas = num2cell(datas);
                captions = {'X [nm]','Y [nm]','frame'}; 
            else
                sigma = sigma(F_UPS);
                datas = [XY_UPS/upscale_fac*pix_size F_UPS sigma];
                datas = num2cell(datas);
                captions = {'X [nm]','Y [nm]','frame','sigma'};
            end
end
%-------------------------------------------------------------------------%
        
    end % methods
            
end