function [datas, captions, table_names, fig] = analyze_SIFNE(obj,~,~) 

    sgm = obj.do_SIFNE_Segmentation(false);
    AllFragments = squeeze(sgm(:,:,1,1,1)); 
    ROI_Mask = squeeze(sgm(:,:,2,1,1));
    ref_img = squeeze(sgm(:,:,3,1,1));
    clear('sgm');
    
    datas = [];
    captions = []; 
    table_names = [];
    fig = [];
    if 0 == sum(AllFragments,'All'), return, end
    
    save_dir = obj.DefaultDirectory;
    
    fname = obj.current_filename;
    fname = strsplit(fname,'.');
    fname = fname{1};
    
    R = obj.SIFNE_LFT_OFT_Radius_of_Filter;

            % Remove Junctions & Single Points button
            NOptsMargin = R;

            AllFragments(1:NOptsMargin,:) = 0;
            AllFragments((end-NOptsMargin):end,:) = 0;
            AllFragments(:,1:NOptsMargin) = 0;
            AllFragments(:,(end-NOptsMargin):end) = 0;

            [x,y] = find(AllFragments==1);
            Allpts = [x y];
            % remove margin information above

            % remove crossing points, a N-by-N region located at the base point will be removed
            R_Junc      = max(2,round(obj.SIFNE_SGM_Junction_Size/2)); 
            Size_Junc   = obj.SIFNE_SGM_Junction_Size;
            N1 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)-1,Allpts(:,2)-1));
            N2 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)-1,Allpts(:,2)));
            N3 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)-1,Allpts(:,2)+1));
            N4 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1),Allpts(:,2)-1));
            N5 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1),Allpts(:,2)+1));
            N6 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)+1,Allpts(:,2)-1));
            N7 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)+1,Allpts(:,2)));
            N8 = AllFragments(sub2ind(size(AllFragments),Allpts(:,1)+1,Allpts(:,2)+1));

            N = N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8;

            CrPts = Allpts(find(N>=3),1:2);

            h = waitbar(0,'Removing crossing points');
            for i = 1:size(CrPts,1)
                waitbar(i/size(CrPts,1),h);
                AllFragments(CrPts(i,1)-R_Junc:CrPts(i,1)+R_Junc,CrPts(i,2)-R_Junc:CrPts(i,2)+R_Junc) = 0;
            end

            close(h);
            [x y] = find(AllFragments==1);
            Allpts = [x y];

            N1 = AllFragments(sub2ind(size(AllFragments),x-1,y-1));
            N2 = AllFragments(sub2ind(size(AllFragments),x-1,y));
            N3 = AllFragments(sub2ind(size(AllFragments),x-1,y+1));
            N4 = AllFragments(sub2ind(size(AllFragments),x,y-1));
            N5 = AllFragments(sub2ind(size(AllFragments),x,y+1));
            N6 = AllFragments(sub2ind(size(AllFragments),x+1,y-1));
            N7 = AllFragments(sub2ind(size(AllFragments),x+1,y));
            N8 = AllFragments(sub2ind(size(AllFragments),x+1,y+1));

            N = N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8;

            SinglePts = Allpts(find(N==0),1:2);
            % remove single points
            h = waitbar(0,'Removing single points');
            for i = 1:size(SinglePts,1)
                waitbar(i/size(SinglePts,1),h);
                AllFragments(SinglePts(i,1),SinglePts(i,2)) = 0;
            end
            close(h);
            [x y] = find(AllFragments==1);
            Allpts = [x y];
            [L num] = bwlabel(AllFragments,8);

            RawCrPts = CrPts;
%             figure('name','Individual Filamentous Fragments');
%             OriginImg = imadjust(im2uint8(obj.imgdata));
%             imshow(mat2gray(OriginImg));hold on;plot(y-R,x-R,'r.');axis off;
            %
            % Remove Junctions & Single Points button - end
            
            % remove short line segments
            MIN_FragmentLength = obj.SIFNE_SGM_Minimum_Filaments_Number_of_Pixels;
            AllFragments = bwareaopen(AllFragments, MIN_FragmentLength);
            [L num] = bwlabel(AllFragments,8);
            [x y] = find(AllFragments==1);
            Allpts = [x y];

%             figure('name','Filtered Skeleton (Short Filaments Removed)');
%             imshow(mat2gray(OriginImg));hold on;
%             plot(y-R,x-R,'r.');axis off;            
            % remove short line segments
            
            % optional - iterative extraction of fragments 
            %
            % to do
            %
            % optional - iterative extraction of fragments 
            
            %Registration of tips direction
LL = L;
SkeTermi = bwmorph(AllFragments,'endpoints');
[x,y] = find(SkeTermi==1);
all_tips = [x y];

all_tips(:,3) = L(sub2ind(size(L),all_tips(:,1),all_tips(:,2)));

% may register tip orientation using parallel computing below
RegR = 2*R; % to register the direction, we need to consider only a local region around a tip
%
tempInfo = zeros(size(all_tips,1),3);

n_cores = feature('numCores');
MultiCore = n_cores>=4;
tic;
if ~MultiCore
    h = waitbar(0,'Registering Tips ...');
    for i = 1:size(all_tips,1)
        waitbar(i/size(all_tips,1),h);
        tempInfo(i,:) = TipReg(LL, all_tips(i,3), all_tips(i,1:2), RegR);
    end
    close(h);
else
    %delete(gcp);    
    %parpool ('local',feature('numCores'));
    parfor i = 1:size(all_tips,1)
        tempInfo(i,:) = TipReg(LL, all_tips(i,3), all_tips(i,1:2), RegR);
        i
    end
    %delete(gcp);
end
toc;
all_tips = [all_tips, tempInfo];
% may register tip orientation using parallel computing above

% figure('name','All Tips Detected (Orientations of 500 Tips Have Been Shown)');
% imshow(mat2gray(AllFragments));hold on; axis off;plot(all_tips(:,2),all_tips(:,1),'r+');
% for i = 1:100
%     idx = ceil(rand * size(all_tips,1));
%     text(all_tips(idx,2),all_tips(idx,1),num2str(all_tips(idx,4)),'color','g');
% end

            %Registration of tips direction
            
% fileID = fopen('UserSettings\GroupingSettings.txt','w');
% fprintf(fileID,['Radius of Searching Fan (pixels):              ',num2str(FanR),'\r\n']);
% fprintf(fileID,['Angle of Searching Fan (degrees):              ',num2str(EditFanAngle),'\r\n']);
% fprintf(fileID,['Criterion1 (Orientation Difference) (degrees): ',num2str(C1),'\r\n']);
% fprintf(fileID,['Criterion1 Weight:                             ',num2str(C1weight),'\r\n']);
% fprintf(fileID,['Criterion2 (Distance) (pixels):                ',num2str(C2),'\r\n']);
% fprintf(fileID,['Criterion3 (Gap Orientation) (degrees):        ',num2str(C3),'\r\n']);
% fprintf(fileID,['Criterion3 Weight:                             ',num2str(C3weight),'\r\n']);
% fprintf(fileID,['Pixel Size (um):                               ',num2str(EditPixelSize),'\r\n']);
% fprintf(fileID,['Short Fragments to Remove (pixels):            ',num2str(ShortFilamentEdit),'\r\n']);
% fprintf(fileID,['Maximum Curvature (radian/um):                 ',num2str(MaxCur),'\r\n']);
                                                        
% tips pairing - start

% auto set condition  
PixelSize = obj.microns_per_pixel; % str2num(get(handles.EditPixelSize,'String'));
pixelsize   = obj.microns_per_pixel; % :)
MaxCur      = obj.SIFNE_Max_Curvature;
FanRadius   = 1/MaxCur/pixelsize;
FanAngle    = 360/2/pi;     % C1          
                            % C3 = FanAngle/2
                            
MIN_Angle_Diff = obj.SIFNE_Orientation_Difference;      % str2num(get(handles.C1,'String'));
MIN_Dist = FanRadius;                                   % obj.SIFNE_Search_Radius; %str2num(get(handles.FanR,'String'));
                                                        % set(handles.ShortFilamentEdit,'string',MIN_Dist);
FanR = FanRadius;                                                        
                                                        
MIN_GapAngle_Diff = obj.SIFNE_Orientation_Difference;   % str2num(get(handles.C3,'String'));

MIN_Info = [MIN_Angle_Diff,MIN_Dist,MIN_GapAngle_Diff];

label_list = all_tips(:,3); % will be updated; help to find another tip of non-first filament during searching
L_GlobalIndex = zeros(size(L));
all_tips(:,7) = (1:size(all_tips,1))';
L_GlobalIndex(sub2ind(size(L_GlobalIndex),all_tips(:,1),all_tips(:,2))) = all_tips(:,7);
% go through all tips and register all searched information

%MaxCur = str2num(get(handles.MaxCurEdit,'String'));

new_partner_list = zeros(size(all_tips,1),100);         % reserve the memory; a maximum of 100 partner tips allowed
C1weight = obj.SIFNE_Orientation_Difference_Weight;     % str2num(get(handles.C1weightEdit,'String'));
C3weight = obj.SIFNE_Gap_Orientation_Weight;            % str2num(get(handles.C3weightEdit,'String'));

tic;
if ~MultiCore==1
    h = waitbar(0,'Tip Searching in Progress ...');
    for i = 1:size(all_tips,1)
        waitbar(i/size(all_tips,1),h);
        new_partner_list(i,:) = local_search_YA(L_GlobalIndex, ...            % image of skeleton and its tips are labeled with global index
            label_list, ...               % this list will be frequently updated (if a tip has been used, it will be removed from this list)
            FanAngle, ...                 % size of the fan-shape searching region
            all_tips, ...                 % information of all tips
            FanR,  ...                    % this radius defines the range for searching
            all_tips(i,4), ...            % orientation of that current tip
            all_tips(i,1:2), ...          % coordinate of current tip
            MIN_Info, ...                 % Minimum conditions that should satisfy
            C1weight, ...                 % weight of croterion 1
            C3weight);                    % weight of croterion 3
        
    end
    close(h);
else
    % delete(gcp);
    % parpool ('local',feature('numCores'));
    % parellel computing is used and it may take a few minutes for large data set
    parfor i = 1:size(all_tips,1)
        new_partner_list(i,:) = local_search_YA(L_GlobalIndex, ...            % image of skeleton and its tips are labeled with global index
            label_list, ...               % this list will be frequently updated (if a tip has been used, it will be removed from this list)
            FanAngle, ...                 % size of the fan-shape searching region
            all_tips, ...                 % information of all tips
            FanR,  ...                    % this radius defines the range for searching
            all_tips(i,4), ...            % orientation of that current tip
            all_tips(i,1:2), ...          % coordinate of current tip
            MIN_Info, ...                 % Minimum conditions that should satisfy
            C1weight, ...                 % weight of croterion 1
            C3weight);                    % weight of croterion 3
        
    end
    % delete(gcp);
end
toc;

for i = 1:size(new_partner_list,2)
    if sum(new_partner_list(:,i))==0
        new_partner_list(:,i:end) = [];
        break;
    end
end
%save data\new_partner_list.mat new_partner_list;

all_tips(:,7) = 1:size(all_tips,1); % add global index to all tips
all_tips(:,8) = ones(size(all_tips,1),1); % this reserves to indicates number of lives

Overlap = strcmp(obj.SIFNE_Filament_Overlap,'None (For Intricate Network)'); %get(handles.OverlapList,'Value');

h = waitbar(0,'Add Number of Lives for Each Filament ...');
% add the number of lives & partner tips to all tips
all_tips(:,9:end) = []; % clear previous record
for i = 1:size(all_tips,1)
    waitbar(i/size(all_tips,1),h);
    temp_partner = new_partner_list(i,:);
    temp_partner(find(temp_partner==0)) = [];
    if ~isempty(temp_partner)
        all_tips(i,8) = length(temp_partner);
        all_tips(i,9:(8 + length(temp_partner))) = temp_partner;
    end
end
close(h);
% add the number of lives & partner tips to all tips
[L num] = bwlabel(AllFragments,8);
L_GlobalIndex = zeros(size(L));
L_GlobalIndex(sub2ind(size(L),all_tips(:,1),all_tips(:,2))) = 1:length(all_tips(:,1));
%
disp('Tip Search Done ! Please Proceed to GROUPING !');

% assign number of lives to each fragments according to the max number of lives of its two tips
for i = 1:num
    pair_index = find(all_tips(:,3)==i);
    if isempty(pair_index)
        continue;
    end
    all_tips(pair_index(1),8) = max(all_tips(pair_index,8));
    all_tips(pair_index(2),8) = max(all_tips(pair_index,8));
end

all_tips = biDirPairing(all_tips);

if Overlap == 1
    all_tips(:,8) = ones(size(all_tips,1),1); % if overlap is not allowed, the number of lives should be '1'
end

% ****** structure of all_tips so far ******   #: Number
%           colume1            colume2           colume3         colume4         colume5            colume6        colume7          colume8          colume9               colume10                colume11          ...
% tip1     # of Row           # of Col           Labeled #     Orientation     Row # Center      Col # Center   global index      # of Lives     index of partner1     index of partner2       index of partner3     ...
% tip2     # of Row           # of Col           Labeled #     Orientation     Row # Center      Col # Center   global index      # of Lives     index of partner1     index of partner2       index of partner3     ...
% tip3     # of Row           # of Col           Labeled #     Orientation     Row # Center      Col # Center   global index      # of Lives     index of partner1     index of partner2       index of partner3     ...
% tip4     # of Row           # of Col           Labeled #     Orientation     Row # Center      Col # Center   global index      # of Lives     index of partner1     index of partner2       index of partner3     ...
% ......
% ****** structure of all_tips so far ******   #: Number
% assign number of lives to each fragments according to the max number of lives of its two tips

% tips pairing - ends
            
% Grouping - start
NofLives = all_tips(:,8);
NofLives(find(NofLives==0)) = 1;
all_tips(:,8) = NofLives;
% assign number of lives to each fragments according to the max number of lives of its two tips

new_partner_list = all_tips(:,9:end); % update the list of partners

all_filament = [];
all_connects = [];
[x y] = find(L~=0);
Lpts = [x, y, L(sub2ind(size(L),x,y))];
h = waitbar(0,'Grouping in Progress');
tic;
for i = 1:num
    waitbar(i/num,h);
    L_list = [];
    label_list = all_tips(:,3); % labels (corresponding to L) of all fragments the global index of which corresponds to all_tips
    xx = Lpts(find(Lpts(:,3)==i),1);
    yy = Lpts(find(Lpts(:,3)==i),2);
    if ~isempty(xx) % if this fragment still exists
        tips_index = find(all_tips(:,3)==i);% find global index
        newSinglefilament = []; % used to store a single filament
        newSinglefilament = [newSinglefilament; [xx yy]]; % add up to the single filament
        L_list = [L_list  i];
        tip_dir1 = tips_index(1); % change to the first tip
        tip_dir2 = tips_index(2);
        connect_tips = all_tips(tip_dir1,1:2); % the design of this connect_tips is for quick filling between gaps later
        if NofLives(tip_dir1)<2
            Lpts(find(Lpts(:,3)==i),3)=0;
            label_list(tips_index(1)) = 0;
            label_list(tips_index(2)) = 0;
            new_partner_list(find(new_partner_list==tip_dir1)) = 0;
            new_partner_list(find(new_partner_list==tip_dir2)) = 0;
        else
            NofLives(tips_index(1)) = NofLives(tips_index(1)) - 1; % reduce its number of lives by 1 if it is used
            NofLives(tips_index(2)) = NofLives(tips_index(2)) - 1; % reduce its number of lives by 1 if it is used
        end
        while 1+1==2
            can_index = new_partner_list(tip_dir1,:); % get global indices of possible partner tips
            can_index(find(can_index==0)) = []; % remove zeros
            if isempty(can_index) || ~isempty(find(L_list==all_tips(can_index(1),3)))  % if the current tip has no possible partner tip
                break;
            else
                optimal_index = can_index(1);
                xx = Lpts(find(Lpts(:,3)==label_list(optimal_index)),1);
                yy = Lpts(find(Lpts(:,3)==label_list(optimal_index)),2);
                newSinglefilament = [newSinglefilament; [xx yy]]; % add new fragment to the current single filament
                L_list = [L_list  label_list(optimal_index)];
                connect_tips = [connect_tips;all_tips(optimal_index,1:2)];
                ll = NofLives(optimal_index);
                if NofLives(optimal_index)<2
                    Lpts(find(Lpts(:,3)==label_list(optimal_index)),3)=0;
                    new_partner_list(find(new_partner_list==tip_dir1)) = 0;
                    label_list(optimal_index) = 0; % remove this label but leave the other one
                    new_partner_list(find(new_partner_list==optimal_index)) = 0;
                else
                    NofLives(find(label_list==all_tips(optimal_index,3))) = NofLives(find(label_list==all_tips(optimal_index,3))) - 1; % reduce its number of lives by 1 if it is used
                end
                new_partner_list(tip_dir1,  find(new_partner_list(tip_dir1,:)==optimal_index)  ) = 0; % remove the tip from the current list of partner tips
                current_label = all_tips(optimal_index,3); % get the current label
                tip_dir1 = find(label_list==current_label); % find the other tip since only one left now
                if length(tip_dir1)==2
                    tip_dir1 = tip_dir1(find(tip_dir1~=optimal_index)); % find the other end of the fragment
                end
                if ll<2
                    new_partner_list(find(new_partner_list==tip_dir1)) = 0;
                    label_list(tip_dir1) = 0;
                end
                connect_tips = [connect_tips;all_tips(tip_dir1,1:2)];
            end
        end
        tip_dir2 = tips_index(2);
        connect_tips = [all_tips(tip_dir2,1:2);connect_tips];
        while 1+1==2
            can_index = new_partner_list(tip_dir2,:); % get global indices of possible partner tips
            can_index(find(can_index==0)) = []; % remove zeros
            if isempty(can_index) || ~isempty(find(L_list==all_tips(can_index(1),3)))% if the current tip has no possible partner tip
                break;
            else
                optimal_index = can_index(1);
                xx = Lpts(find(Lpts(:,3)==label_list(optimal_index)),1);
                yy = Lpts(find(Lpts(:,3)==label_list(optimal_index)),2);
                newSinglefilament = [newSinglefilament; [xx yy]]; % add new fragment to the current single filament
                L_list = [L_list  label_list(optimal_index)];
                connect_tips = [all_tips(optimal_index,1:2);connect_tips];
                ll = NofLives(optimal_index);
                if NofLives(optimal_index)<2
                    Lpts(find(Lpts(:,3)==label_list(optimal_index)),3)=0;
                    new_partner_list(find(new_partner_list==tip_dir2)) = 0;
                    label_list(optimal_index) = 0; % remove this label but leave the other one
                    new_partner_list(find(new_partner_list==optimal_index)) = 0;
                else
                    NofLives(find(label_list==all_tips(optimal_index,3))) = NofLives(find(label_list==all_tips(optimal_index,3))) - 1; % reduce its number of lives by 1 if it is used
                end
                new_partner_list(tip_dir2,  find(new_partner_list(tip_dir2,:)==optimal_index)  ) = 0; % remove the tip from the current list of partner tips
                current_label = all_tips(optimal_index,3); % get the current label
                tip_dir2 = find(label_list==current_label); % find the other tip since only one left now
                if length(tip_dir2)==2
                    tip_dir2 = tip_dir2(find(tip_dir2~=optimal_index)); % find to other end of the fragment
                end
                if ll<2
                    new_partner_list(find(new_partner_list==tip_dir2)) = 0;
                    label_list(tip_dir2) = 0;
                end
                connect_tips = [all_tips(tip_dir2,1:2);connect_tips];
            end
        end
        if isempty(all_filament) && isempty(all_connects)
            all_filament(1:size(newSinglefilament,1),1:2,1) = newSinglefilament;
            all_connects(1:size(connect_tips,1),1:2,1) = connect_tips;
        else
            all_filament(1:size(newSinglefilament,1),1:2,size(all_filament,3)+1) = newSinglefilament;
            all_connects(1:size(connect_tips,1),1:2,size(all_connects,3)+1) = connect_tips;
        end
    end
end
toc;

disp('Grouping Done !');

close(h);

% filamentous fragment grouping ends

% Grouping - ends

% Sorting - start
maskI = zeros(size(L));
Fullength = 5*size(all_filament,1);
all_sorted_filament = zeros(Fullength,2,size(all_filament,3));

tic;
if ~MultiCore
    h = waitbar(0,'Analysis in Progress ...');
    for i = 1:size(all_filament,3)
        waitbar(i/size(all_filament,3),h);
        [all_sorted_filament(:,:,i)  AnalysisInfo(i,:)]= SortFilament(all_filament(:,:,i), maskI, all_connects(:,:,i), Fullength);
    end
    close(h);
else
    %delete(gcp);
    %parpool ('local',feature('numCores'));
    parfor i = 1:size(all_filament,3)
        [all_sorted_filament(:,:,i)  AnalysisInfo(i,:)]= SortFilament(all_filament(:,:,i), maskI, all_connects(:,:,i), Fullength);
        i
    end
    %delete(gcp);
end
toc;
for i = 1:size(all_sorted_filament,1)
    temp = all_sorted_filament(i,:,:);
    if(sum(temp(:)))==0
        all_sorted_filament(i:end,:,:) = [];
        break;
    end
end
%
disp('Analysis and Sorting Done !');

% Sorting - ends

% remove short filaments
ShortFilament = obj.SIFNE_Sorting_Minimum_Filament_Size;

RemoveIdx = [];

h = waitbar(0,'Removing Short Filaments...');
for i = 1:size(all_sorted_filament,3)
    waitbar(i/size(all_sorted_filament,3),h);
    if AnalysisInfo(i,2) <= ShortFilament
        RemoveIdx = [RemoveIdx i];
    end
end
close(h);
all_sorted_filament(:,:,RemoveIdx) = [];
AnalysisInfo(RemoveIdx,:) = [];
all_connects(:,:,RemoveIdx) = [];

% % remove ungrouped below
% RemoveUngrp = 1; % get(handles.RemoveUngrp,'Value');
% if RemoveUngrp==1
%     AllFragments = im2bw(AllFragments);
%     RemoveIdx = [];
%     h = waitbar(0,'Removing Ungrouped Fragments...');
%     for i = 1:size(all_sorted_filament,3)
%         waitbar(i/size(all_sorted_filament,3),h);
%         temp = all_sorted_filament(:,:,i);
%         temp(find(temp(:,1)==0),:) = [];
%         if size(temp,1)==sum(AllFragments(sub2ind(size(AllFragments),temp(:,1),temp(:,2))))
%             RemoveIdx = [RemoveIdx  i];
%         end
%     end
%     close(h);
%     all_sorted_filament(:,:,RemoveIdx) = [];
%     AnalysisInfo(RemoveIdx,:) = [];
%     all_connects(:,:,RemoveIdx) = [];
% end

% remove ungrouped above

% remove short filaments

% COMMENTED - PREV. ANALYSES - STARTS
% Analyses

% % all_sorted_filament
% % AllFragments
% % ROI_Mask
% % R
% %Analysis 1
%         warning off;
% %         close(figure(1));close(figure(2));
% %         load data\R.mat;
% %         load data\AllFragments.mat;
% %         load data\L.mat;
% %         load data\all_sorted_filament.mat all_sorted_filament;
% %         load data\AnalysisInfo.mat;
% %         load data\ROI_Mask;
%         % reconstruct network
%         allpts_IncludeDup = [];
%         Overlay_Map = zeros(size(L));
%         h = waitbar(0,'Retrieve All Filaments');
%         for i = 1:size(all_sorted_filament,3)
%             waitbar(i/size(all_sorted_filament,3),h);
%             temp1 = all_sorted_filament(:,:,i); % get all coordinates
%             temp1(find(temp1(:,1)==0),:) = []; % remove zeros
%             allpts_IncludeDup = [allpts_IncludeDup;temp1];
%             % check points used more than once
%             Overlay_Map(sub2ind(size(L),temp1(:,1),temp1(:,2))) = Overlay_Map(sub2ind(size(L),temp1(:,1),temp1(:,2))) + 1;
%         end
%         close(h);
%         
%         h22=figure(22);
%         imshow(mat2gray(AllFragments((R+1):(size(AllFragments,1)-R), (R+1):(size(AllFragments,2)-R))));hold on;axis off;
%         screen_size = get(0, 'ScreenSize'); 
%         set(h22, 'Position', [0 0 screen_size(3) screen_size(4) ] );
%         ColorList = rand(32,3);
%         temp = 1:size(all_sorted_filament,3);
%         for i = 1:32
%             idx = find(mod(temp,32)==(i-1));
%             temp1 = all_sorted_filament(:,1,idx);  temp1 = temp1(:);  temp1(find(temp1==0)) = [];
%             temp2 = all_sorted_filament(:,2,idx);  temp2 = temp2(:);  temp2(find(temp2==0)) = [];
%             plot(temp2-R,temp1-R,'.','color',ColorList(i,:),'MarkerSize',6);hold on;
%         end
%         %Overlap = get(handles.OverlapList,'Value');
%         %if Overlap==2
%         if Overlap~=1
%             [x y] = find(Overlay_Map>1);
%             plot(y-R,x-R,'r.','MarkerSize',6);hold on;
%         end
%         B = bwboundaries(ROI_Mask);
%         for k=1:numel(B)
%             B_k = B{k};
%             plot(B_k(:,2),B_k(:,1),'.','color',[1 1 1],'MarkerSize',6);hold on;
%         end          
%         set(gcf, 'Name', fname, 'numbertitle', 'off');
%         saveas(h22,[save_dir filesep fname '_' 'Extracted_Filaments.fig']);
%         close(h22);
% %Analysis 1
% 
% % all_sorted_filament
% % PixelSize
% % ROI_Mask
% % R
% %Analysis 2
%         warning off;
% %         load data\R.mat;
% %         load data\OriginImg.mat;
% %         load data\L.mat;
% %         load data\all_sorted_filament.mat all_sorted_filament;
% %         load data\ROI_Mask;
% 
%         PixelSize = obj.microns_per_pixel; % str2num(get(handles.EditPixelSize,'String'));
%                 
%         Overlay_Map = zeros(size(L));
%         h = waitbar(0,'Checking Junctions ...');
%         for i = 1:size(all_sorted_filament,3)
%             waitbar(i/size(all_sorted_filament,3),h);
%             temp1 = all_sorted_filament(:,:,i);
%             temp1(find(temp1(:,1)==0),:) = [];
%             Overlay_Map(sub2ind(size(L),temp1(:,1),temp1(:,2))) = Overlay_Map(sub2ind(size(L),temp1(:,1),temp1(:,2))) + 1;
%         end
%         close(h);
%         [x y] = find(Overlay_Map>1);
%         
%         CroMap = zeros(size(L));
%         CroMap(sub2ind(size(L),x,y)) = 1;
%         for i = 1:length(x)
%             temp = CroMap((x(i)-1):(x(i)+1), (y(i)-1):(y(i)+1));
%             if sum(temp(:))>1
%                 CroMap(x(i),y(i)) = 0;
%             end
%         end
%         [x y] = find(CroMap==1);
%         NewCrPts = [x y];
%         %save data\NewCrPts.mat NewCrPts;
%         
%         disF=bwdist(bwmorph(ROI_Mask,'remove'));
%         mask = ROI_Mask;
%         mask = single(mask);
%         mask(mask==0)=0;
%         disF = disF.*mask;
%         disF = disF*PixelSize;
%         
%         h22=figure(22);
%         imagesc(disF);colormap(jet);hold on;axis off;axis image;colorbar;
%         [xx yy] = find(Overlay_Map~=0);
%         plot(yy,xx,'k.');
%         plot(y,x,'g.','MarkerSize',15);title('Distribution of Junctions');
%         B = bwboundaries(ROI_Mask);
%         for k=1:numel(B)
%             B_k = B{k};
%             plot(B_k(:,2),B_k(:,1),'.','color',[1 1 1],'MarkerSize',6);hold on;
%         end        
%         
%         disF=bwdist(bwmorph(ROI_Mask,'remove'));
%         mask = ROI_Mask;
%         mask = single(mask);
%         mask(mask==0)=-1;
%         disF = disF.*mask;
%         d = disF(sub2ind(size(disF),x,y))*PixelSize;
%         set(gcf, 'Name', fname, 'numbertitle', 'off');
%         saveas(h22,[save_dir filesep fname '_' 'Distribution_of_Junctions.fig']);
%         close(h22);
% 
%         if ~isempty(d)
%             h22 = figure(22);
%             %         histogram(d);%axis square;
%             histfit(d,50,'kernel');xlim([0 inf]);
%             ylabel('Frequency');
%             xlabel('Distance to Cell Edge (\mum)');
%             title('Distribution of Junctions');
%             set(gcf, 'Name', fname, 'numbertitle', 'off');
%             saveas(h22,[save_dir filesep fname '_' 'Distribution_of_Junctions_Analysis.fig']);        
%             close(h22);
%         end
% %Analysis 2
% 
% %AnalysisInfo
% %ROI_Mask
% %all_sorted_filament
% %Analysis 3
% %         warning off;
% %         close(figure(1));close(figure(2));
% %         
% %         load data\all_sorted_filament.mat;
% %         load data\AnalysisInfo.mat;
% %         load data\L;
% %         load data\ROI_Mask;
%         % plot of information of all filaments
%         %PixelSize = str2num(get(handles.EditPixelSize,'String'));
%         h22 = figure;
%         h=rose(pi*AnalysisInfo(:,1)/180,30);axis off;
%         set(h,'linewidth',3)
%         axis square;
%         daspect([1 1 1]);
%         title('Histogram of Filaments Orientations');
%         set(gcf, 'Name', fname, 'numbertitle', 'off');
%         saveas(h22,[save_dir filesep fname '_' 'Histogram_of_Orientations.fig']);
%         close(h22);
%         
% %         mkdir result;
% %         saveas(figure(1),'result\Histogram_of_Orientations.fig');
%         
%         allc = AnalysisInfo(:,4:5);
%         disF=bwdist(bwmorph(ROI_Mask,'remove'));
%         mask = ROI_Mask;
%         mask = single(mask);
%         mask(mask==0)=-1;
%         disF = disF.*mask;
%         d = disF(sub2ind(size(L),ceil(allc(:,1)),ceil(allc(:,2))))*PixelSize;
%         
%         finalmap = zeros(181,ceil(max(d)));
%         
%         h = waitbar(0,'generating colormap ...');
%         for i = 1:length(d)
%             waitbar(i/size(all_sorted_filament,3),h);
%             if ceil(d(i))>0
%                 finalmap(ceil(AnalysisInfo(i,1)+91), ceil(d(i))) = finalmap(ceil(AnalysisInfo(i,1)+91), ceil(d(i))) + 1;
%             end
%         end
%         close(h);
%         
%         h22 = figure(22);
%         imagesc(finalmap);colormap(jet);colorbar;
%         set(gca,'ytick',[]);
%         xlabel('Distance to Cell Edge (\mum)');ylabel('-90 degrees to 90 degrees');
%         title('Distribution of Filament Orientations');
%         set(gcf, 'Name', fname, 'numbertitle', 'off');
%         saveas(h22,[save_dir filesep fname '_' 'Distribution_of_Orientations.fig']);
%         close(h22);
% %Analysis 3
% 
% %FanR
% %AllFragments
% %all_sorted_filament
% %PixelSize
% %ROI_Mask
% % Analysis 4
% %         warning off;
% %         close(figure(1));close(figure(2));
% %         load data\all_sorted_filament;
% %         load data\ROI_Mask;
% %         load data\AllFragments;
%         curR = round(FanR/2); % round(str2num(get(handles.FanR,'String'))/2);
%         
%         disF=bwdist(bwmorph(ROI_Mask,'remove'));
%         mask = ROI_Mask;
%         mask = single(mask);
%         mask(mask==0)=-1;
%         disF = disF.*mask;
%         
%         allcur = [];
%         alldSDF = [];
%         M = zeros(size(AllFragments));
%         h = waitbar(0,'Calculating Curvatures ...');
%         all_filament_curs = [];
%         all_mean_cur = zeros(1,size(all_sorted_filament,3));
%         for i = 1:size(all_sorted_filament,3)
%             waitbar(i/size(all_sorted_filament,3),h);
%             F = all_sorted_filament(:,:,i);
%             F(find(F(:,1)==0),:) = [];
%             Dlist = 0;
%             for j = 2:size(F,1)
%                 Dlist = [Dlist  Dlist(j-1)+pdist([F(j-1,:); F(j,:)])];
%             end
%             NoCurlist = unique([find(Dlist<=curR)  find(Dlist>=(Dlist(end)-curR))]);
%             if length(NoCurlist)~=size(F,1)
%                 for j = 1:size(F,1)
%                     if isempty(find(NoCurlist==j))
%                         x2 = F(j,1);
%                         y2 = F(j,2);
%                         D = Dlist - Dlist(j);
%                         D1 = abs(D - curR);
%                         idx1 = find(D1==min(D1));
%                         D2 = abs(D - (-curR));
%                         idx2 = find(D2==min(D2));
%                         x1 = F(idx1,1);
%                         y1 = F(idx1,2);
%                         x3 = F(idx2,1);
%                         y3 = F(idx2,2);
%                         alldSDF = [alldSDF  disF(sub2ind(size(disF),x2,y2))*PixelSize];
%                         cur = 2*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1)) ./sqrt(((x2-x1).^2+(y2-y1).^2)*((x3-x1).^2+(y3-y1).^2)*((x3-x2).^2+(y3-y2).^2));
%                         allcur = [allcur  cur/PixelSize];
%                         M(F(j,1),F(j,2)) = cur;
%                         all_filament_curs(j,1,i) = x2;
%                         all_filament_curs(j,2,i) = y2;
%                         all_filament_curs(j,3,i) = cur/PixelSize;
%                         all_filament_curs(j,4,i) = 1;% flag for calculation done for this point
%                     end
%                 end
%                 k = find(all_filament_curs(:,4,i)==1);
%                 all_mean_cur(i) = mean(all_filament_curs(k,3,i));
%             end
%         end
%         close(h);
%         M = M/PixelSize;
%         h22 = figure(22);
%         imagesc(M); colormap(jet);colorbar; title('Distribution of Curvatures');axis off;axis image;
%         set(gcf, 'Name', fname, 'numbertitle', 'off');
%         saveas(h22,[save_dir filesep fname '_' 'Distribution_of_Curvatures.fig']);
%         close(h22);
%           
%         All_Curs = [];
%         All_Curs_4plot = [];
%         disF=bwdist(bwmorph(ROI_Mask,'remove'));
%         mask = ROI_Mask;
%         mask = single(mask);
%         mask(mask==0)=-1;
%         disF = disF.*mask;
%         %
%         h = waitbar(0,'Calculating the Distribution of Curvatures ...');
%         for m = 1:size(all_filament_curs,3)
%             waitbar(m/size(all_filament_curs,3),h);
%             Cur = all_filament_curs(:,:,m);
%             F = all_sorted_filament;
%             F(find(F(:,1)==0),:) = [];
%             Cur(find(Cur(:,4)==0),:) = [];
%             if isempty(Cur)
%                 All_Curs_4plot = [All_Curs_4plot;  disF(sub2ind(size(disF),round((F(1,1)+F(end,1))/2),round((F(end,2)+F(1,2))/2)))   0];
%                 All_Curs = [All_Curs, 0];
%             else
%                 All_Curs_4plot = [All_Curs_4plot;  disF(sub2ind(size(disF),Cur(:,1),Cur(:,2)))   Cur(:,3)];
%                 All_Curs = [All_Curs, Cur(:,3)'];
%             end
%         end
%         close(h);
%         All_Curs_4plot(:,1) = round(All_Curs_4plot(:,1));
%         temp_mean_cur = [];
%         x1 = [];
%         y1 = [];
%         for i = min(All_Curs_4plot(:,1)): max(All_Curs_4plot(:,1))
%             k = find(All_Curs_4plot(:,1)==i);
%             temp_mean_cur = [temp_mean_cur  mean(All_Curs_4plot(k,2))];
%             if length(k)>5
%                 x1 = [x1  i];
%                 y1 = [y1  mean(All_Curs_4plot(k,2))];
%             end
%         end
% 
%         h22 = figure(22);
%         subplot(1,2,1);
%         histfit(all_mean_cur,50,'kernel');xlim([0 inf]);
%         xlabel('Curvature (unit: \mum^-^1)');ylabel('Frequency');
%         title('Distribution of Means of Filament Curvatures');
%         %        
%         
%         subplot(1,2,2);
%         plot(x1*0.02, y1,  'r');hold on;axis([0 inf 0 inf]);
%         xlabel('Distance to Cell Edge (unit: \mum)');ylabel('Mean Curvature (unit: \mum^-^1)');
%         title('Distribution of Curvatures');
%         set(gcf, 'Name', fname, 'numbertitle', 'off');
%         saveas(h22,[save_dir filesep fname '_' 'Distribution_of_Curvatures_Analysis.fig']);
%         close(h22);
%         %save data\all_filament_curs.mat all_filament_curs;
%         %mkdir result;
%         %saveas(figure(1),'result\Distribution_of_Curvatures.fig');
%         %saveas(figure(2),'result\Distribution_of_Curvatures_Analysis.fig');
% 
% % Analysis 4

% Analysis 5 %xls saving...
%         warning off;
%         close(figure(1));close(figure(2));
%         load data\AnalysisInfo;
%         load data\all_sorted_filament;
%         load data\NewCrPts;
%         load data\Size_Junc;
%         load data\ROI_Mask;
%         load data\L.mat;
%         load data\all_tips.mat;
%         load data\all_connects.mat;
%         load data\all_connects_shortlist.mat;
        
%         % generate information of all filaments
%         % PixelSize = str2num(get(handles.EditPixelSize,'String'));
%         titles = {'Filament ID','1st X pos','1st Y pos','last X pos','last Y pos','Orientation','Total Length','End-to-End Distance','Centroid X','Centroid Y'};
%         titles = [titles,repmat({'X','Y'},[1 size(all_sorted_filament,1)]) ];
%         InfoExcel = zeros(size(all_sorted_filament,3)*2, size(all_sorted_filament,1)*2+9);
%         
%         h = waitbar(0,'Generating information of composite filaments ...');
%         for i = 1:size(all_sorted_filament,3)
%             waitbar(i/size(all_sorted_filament,3),h);
%             curr_filament = all_sorted_filament(:,:,i);
%             curr_filament(find(curr_filament(:,1)==0),:) = [];
%             InfoExcel((i*2-1), 1:9) = [curr_filament(1,:), curr_filament(end,:), AnalysisInfo(i,:)];
%             temp = zeros(1,size(all_sorted_filament,1)*2);
%             temp((1:size(all_sorted_filament,1))*2-1) = all_sorted_filament(:,1,i);
%             temp((1:size(all_sorted_filament,1))*2) = all_sorted_filament(:,2,i);
%             InfoExcel((i*2-1), 10:end) = temp;
%             ks = dsearchn(curr_filament,NewCrPts);
%             for j = 1:length(ks)
%                 if pdist([curr_filament(ks(j),:); NewCrPts(j,:)])<=Size_Junc
%                     InfoExcel((i*2), ((9+ks(j)*2-1):(9+ks(j)*2))) = NewCrPts(j,:);
%                 end
%             end
%         end
%         FragID = 1:size(all_sorted_filament,3);
%         FragID = [FragID; FragID];
%         FragID = FragID(:);
%         InfoExcel = [FragID, InfoExcel];
%         %save data\InfoExcel.mat InfoExcel;
%         close(h);
%         InfoExcel = [titles;num2cell(InfoExcel)];
%         
%         
%         % ------ generate fragment information ------
%         FragmentInfo = zeros(num,size(L,1));
% 
%         if ~MultiCore
%             h = waitbar(0,'Fragment Information ...');
%             for i = 1:num
%                 waitbar(i/num,h);
%                 FragmentInfo(i,:) = GenFragmentInfo(L,i,all_tips,FragmentInfo(i,:));
%             end
%             close(h);
%         else
%             %delete(gcp);
%             %parpool ('local',feature('numCores'));            
%             % parellel computing is used and it may take a few minutes for large data set
% %             load data\L.mat;
% %             L = L;
% %             load data\all_tips.mat;
% %             all_tips = all_tips;
%             parfor i = 1:num
%                 FragmentInfo(i,:) = GenFragmentInfo(L,i,all_tips,FragmentInfo(i,:));
%             end
%         end
%         
%         
%         sumlist = sum(FragmentInfo,1);
%         FragmentInfo(:,(find(sumlist==0))) = [];
%         %save data\FragmentInfo.mat FragmentInfo;
%         MaxFragL = size(FragmentInfo,2)-6;
%         disp('Fragment Information Generated !');
%         titles = {'Fragment ID','# of Pixels','Beginning X','Beginning Y','Ending X','Ending Y'};
%         titles = [titles,repmat({'X','Y'},[1 (size(FragmentInfo,2)-6)/2]) ];
%         FragmentInfo = [titles;num2cell(FragmentInfo)];
%                        
%         % ------ generate linkage information before removing short------
%         
%         LinkInfo = zeros(size(all_connects,1),2+MaxFragL*2,size(all_connects,3));
%         if ~MultiCore
%             h = waitbar(0,'Linkage1 Information ...');
%             for i = 1:size(all_connects,3)
%                 waitbar(i/size(all_connects,3),h);
%                 LinkInfo(:,:,i) = GenLinkageInfo(i,L,all_connects(:,:,i),LinkInfo(:,:,i));
%             end
%             close(h);
%         else
% %             load data\L.mat;
% %             L = L;
%             parfor i = 1:size(all_connects,3)
%                 LinkInfo(:,:,i) = GenLinkageInfo(i,L,all_connects(:,:,i),LinkInfo(:,:,i));
%             end
%         end
%         
%         disp('Linkage1 Information Generated !');
%                 
%         LinkageInfo1 = [];
%         h = waitbar(0,'Integrating Linkage1 Information ...');
%         tic;
%         for i = 1:size(LinkInfo,3)
%             waitbar(i/size(LinkInfo,3),h);
%             tempLinkage = LinkInfo(:,:,i);
%             tempLinkage(find(tempLinkage(:,1)==0),:) = [];
%             LinkageInfo1 = [LinkageInfo1; tempLinkage];
%         end
%         close(h);
%         disp('Linkage1 Integration Done !');
%         sumlist = sum(LinkageInfo1,1);
%         LinkageInfo1(find(LinkageInfo1(:,1)==0),:) = [];
%         LinkageInfo1(:,(find(sumlist==0))) = [];
%         % save data\LinkageInfo1.mat LinkageInfo1;
%         titles = {'Composite Filament ID','Fragment ID'};
%         titles = [titles,repmat({'X','Y'},[1 (size(LinkageInfo1,2)-2)/2]) ];
%         LinkageInfo1 = [titles;num2cell(LinkageInfo1)];
%         
%         
%         % ------ generate linkage information after removing short------
%         
%         LinkInfo = zeros(size(all_connects_shortlist,1),2+MaxFragL*2,size(all_connects_shortlist,3));
%         if ~MultiCore
%             h = waitbar(0,'Linkage2 Information ...');
%             for i = 1:size(all_connects_shortlist,3)
%                 waitbar(i/size(all_connects_shortlist,3),h);
%                 LinkInfo(:,:,i) = GenLinkageInfo(i,L,all_connects_shortlist(:,:,i),LinkInfo(:,:,i));
%             end
%             close(h);
%         else
%             parfor i = 1:size(all_connects_shortlist,3)
%                 LinkInfo(:,:,i) = GenLinkageInfo(i,L,all_connects_shortlist(:,:,i),LinkInfo(:,:,i));
%             end
%         end
%         
%         disp('Linkage2 Information Generated !');
%         
%         
%         LinkageInfo2 = [];
%         h = waitbar(0,'Integrating Linkage2 Information ...');
%         tic;
%         for i = 1:size(LinkInfo,3)
%             waitbar(i/size(LinkInfo,3),h);
%             tempLinkage = LinkInfo(:,:,i);
%             tempLinkage(find(tempLinkage(:,1)==0),:) = [];
%             LinkageInfo2 = [LinkageInfo2; tempLinkage];
%         end
%         close(h);
%         disp('Linkage2 Integration Done !');
%         sumlist = sum(LinkageInfo2,1);
%         LinkageInfo2(find(LinkageInfo2(:,1)==0),:) = [];
%         LinkageInfo2(:,(find(sumlist==0))) = [];
%         % save data\LinkageInfo2.mat LinkageInfo2;
%         titles = {'Composite Filament ID','Fragment ID'};
%         titles = [titles,repmat({'X','Y'},[1 (size(LinkageInfo2,2)-2)/2]) ];
%         LinkageInfo2 = [titles;num2cell(LinkageInfo2)];
%         
%         % mkdir result;
%         % write inforamtion to excel
%         %cd result;
%         
%         xlswrite([save_dir filesep 'IntegratedInfo.xlsx'],InfoExcel,1,'A1');
%         xlswrite([save_dir filesep 'IntegratedInfo.xlsx'],FragmentInfo,2,'A1');
%         xlswrite([save_dir filesep 'IntegratedInfo.xlsx'],LinkageInfo1,3,'A1');
%         xlswrite([save_dir filesep 'IntegratedInfo.xlsx'],LinkageInfo2,4,'A1');
%         %curDir = pwd;
%         %filaname = [curDir,'\IntegratedInfo.xlsx'];
%         e = actxserver('Excel.Application');
%         ewb = e.Workbooks.Open([save_dir filesep 'IntegratedInfo.xlsx']);
%         ewb.Worksheets.Item(1).Name = 'Ultimate Filaments';
%         ewb.Worksheets.Item(2).Name = 'Fragment Info';
%         ewb.Worksheets.Item(3).Name = 'Linkage Info1';
%         ewb.Worksheets.Item(4).Name = 'Linkage Info2';
%         ewb.Save
%         ewb.Close(false);
%         e.Quit;
%         %cd ..;
%         %msgbox('Excel Generated !');
%         disp('Excel Generated !');
        
% Analysis 5


% all_sorted_filament
% AllFragments
% ROI_Mask
% R
% all_sorted_filament
% PixelSize
% ROI_Mask
% R
%AnalysisInfo
%ROI_Mask
%all_sorted_filament
%FanR
%AllFragments
%all_sorted_filament
%PixelSize
% %ROI_Mask
% 
% RawImg = obj.imgdata;
% 
% save([save_dir filesep fname '_raw_analysis_data'],...
% 'all_sorted_filament', ...
% 'AllFragments', ...
% 'ROI_Mask', ...
% 'RawImg', ...
% 'R', ...
% 'PixelSize', ...
% 'AnalysisInfo', ...
% 'FanR');

% Analyses

            % if ~isempty(hw), delete(hw), drawnow; end
            
% COMMENTED - PREV. ANALYSES - ENDS
                        
% refactored analyses            
RawImg = obj.imgdata;
[H,W] = size(RawImg);
intensity_img = zeros(H+R+R,W+R+R);
intensity_img(R+1:H+R,R+1:W+R) = RawImg;

%Overlap = strcmp(obj.SIFNE_Filament_Overlap,'None (For Intricate Network)'); %get(handles.OverlapList,'Value');
Overlap = 1; % 'None (For Intricate Network)'
%Overlap = 2; % 'Allowed'

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% local orientational jitter

        % reconstruct network
        allpts_IncludeDup = [];
        flab = zeros(size(L));
        joints = zeros(size(L));
        h = waitbar(0,'Scanning All Filaments');
        for i = 1:size(all_sorted_filament,3)
            waitbar(i/size(all_sorted_filament,3),h);
            temp1 = all_sorted_filament(:,:,i); % get all coordinates
            temp1(find(temp1(:,1)==0),:) = []; % remove zeros
            % check points used more than once
            flab(sub2ind(size(L),temp1(:,1),temp1(:,2))) = i;
            joints(sub2ind(size(L),temp1(:,1),temp1(:,2))) = joints(sub2ind(size(L),temp1(:,1),temp1(:,2))) + 1;
        end
        close(h);
        
        joints = imdilate(joints>1,strel('disk',2));        
        flab = flab.*(~joints);
                
        tic
        %A = get_blobs_adjacency(flab,[],true); % never       
        A = get_blobs_adjacency(flab,ROI_Mask,true);
        %icy_imshow(A);
        toc/60

loj = zeros(size(A,1),1); %local orientational jitter
labs = (1:size(A,1));
parfor k=1:size(A,1)
    tta_k = AnalysisInfo(k,1);
    adj_k = A(k,:).*labs;
    adj_k = adj_k(adj_k~=0);
    num = 0;
    denom = 0;
    for m=1:numel(adj_k)
        tta_m = AnalysisInfo(m,1);
        l_m = AnalysisInfo(m,2); % use length as weight        
        diff_angle = abs(tta_k - tta_m);
        if diff_angle > 90 
            diff_angle = 180 - diff_angle;
        end
        num = num + diff_angle*l_m;
        denom = denom + l_m;        
    end
    loj(k) = num/denom;
end

%         figure;
%         histogram(loj,'Normalization','probability');
%         xlabel('local orientational jitter [deg]');
%         ylabel('probability');
%         grid(gca,'on');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% local orientational jitter

% %Analysis 1
        % reconstruct network
        allpts_IncludeDup = [];
        Overlay_Map = zeros(size(L));
        h = waitbar(0,'Retrieve All Filaments');
        arrlength = zeros(1,size(all_sorted_filament,3));
        for i = 1:size(all_sorted_filament,3)
            waitbar(i/size(all_sorted_filament,3),h);
            temp1 = all_sorted_filament(:,:,i); % get all coordinates
            temp1(find(temp1(:,1)==0),:) = []; % remove zeros
            allpts_IncludeDup = [allpts_IncludeDup;temp1];
            % check points used more than once
            Overlay_Map(sub2ind(size(L),temp1(:,1),temp1(:,2))) = Overlay_Map(sub2ind(size(L),temp1(:,1),temp1(:,2))) + 1;
            arrlength(i) = size(temp1,1);
        end
        close(h);
        
%         figure;
%         histogram(arrlength*PixelSize,'Normalization','probability');
%         xlabel('filament length [\mum]');
%         ylabel('probability');
%         axis([0 14 0 0.4]);
%         grid(gca,'on');
                
        % intensity on original image (relative to average intensity of filaments)
        % std of this intensity
        %
        sample = ref_img(0~=Overlay_Map);
        mean_filament_intensity = mean(sample(:));
        %
        h = waitbar(0,'Scanning Filaments');
        arr_mean_filament_intensity = zeros(1,size(all_sorted_filament,3));
        arr_std_filament_intensity = zeros(1,size(all_sorted_filament,3));
        for i = 1:size(all_sorted_filament,3)
            waitbar(i/size(all_sorted_filament,3),h);
            temp1 = all_sorted_filament(:,:,i); % get all coordinates
            temp1(find(temp1(:,1)==0),:) = []; % remove zeros
            sample = ref_img(sub2ind(size(L),temp1(:,1),temp1(:,2)));
            arr_mean_filament_intensity(i) = mean(sample(:)); %/mean_filament_intensity;
            arr_std_filament_intensity(i) = std(sample(:)); % /mean_filament_intensity;
        end
        close(h);
        %
%         figure;
%         plot(arr_mean_filament_intensity,arr_std_filament_intensity,'k.','markersize',12);
%         xlabel('relative filament intensity: mean');
%         ylabel('relative filament intensity: std');
%         axis([0 5 0 3.5]);
%         grid(gca,'on');        
                        
%         h99=figure(99);
%         imshow(mat2gray(AllFragments((R+1):(size(AllFragments,1)-R), (R+1):(size(AllFragments,2)-R))));hold on;axis off;
%         screen_size = get(0, 'ScreenSize'); 
%         set(h99, 'Position', [0 0 screen_size(3) screen_size(4) ] );
%         ColorList = rand(32,3);
%         temp = 1:size(all_sorted_filament,3);
%         for i = 1:32
%             idx = find(mod(temp,32)==(i-1));
%             temp1 = all_sorted_filament(:,1,idx);  temp1 = temp1(:);  temp1(find(temp1==0)) = [];
%             temp2 = all_sorted_filament(:,2,idx);  temp2 = temp2(:);  temp2(find(temp2==0)) = [];
%             plot(temp2-R,temp1-R,'.','color',ColorList(i,:),'MarkerSize',6);hold on;
%         end
%         %
%         if Overlap~=1
%             [x y] = find(Overlay_Map>1);
%             plot(y-R,x-R,'r.','MarkerSize',6);hold on;
%         end
%         B = bwboundaries(ROI_Mask);
%         for k=1:numel(B)
%             B_k = B{k};
%             plot(B_k(:,2),B_k(:,1),'.','color',[1 1 1],'MarkerSize',6);hold on;
%         end          
%         set(gcf, 'Name', fname, 'numbertitle', 'off');
% %Analysis 1
% 
% %Analysis 2
                
        Overlay_Map = zeros(size(L));
        h = waitbar(0,'Checking Junctions ...');
        for i = 1:size(all_sorted_filament,3)
            waitbar(i/size(all_sorted_filament,3),h);
            temp1 = all_sorted_filament(:,:,i);
            temp1(find(temp1(:,1)==0),:) = [];
            Overlay_Map(sub2ind(size(L),temp1(:,1),temp1(:,2))) = Overlay_Map(sub2ind(size(L),temp1(:,1),temp1(:,2))) + 1;
        end
        close(h);
        [x y] = find(Overlay_Map>1);
        
        CroMap = zeros(size(L));
        CroMap(sub2ind(size(L),x,y)) = 1;
        for i = 1:length(x)
            temp = CroMap((x(i)-1):(x(i)+1), (y(i)-1):(y(i)+1));
            if sum(temp(:))>1
                CroMap(x(i),y(i)) = 0;
            end
        end
        [x y] = find(CroMap==1);
        NewCrPts = [x y];
        
%         disF=bwdist(bwmorph(ROI_Mask,'remove'));
%         mask = ROI_Mask;
%         mask = single(mask);
%         mask(mask==0)=0;
%         disF = disF.*mask;
%         disF = disF*PixelSize;
%         
%         h33=figure(33);
%         imagesc(disF);colormap(jet);hold on;axis off;axis image;colorbar;
%         [xx yy] = find(Overlay_Map~=0);
%         plot(yy,xx,'k.');
%         plot(y,x,'d','MarkerFaceColor','m','MarkerEdgeColor','k','MarkerSize',8);                
%         title('Distribution of Junctions');
%         B = bwboundaries(ROI_Mask);
%         for k=1:numel(B)
%             B_k = B{k};
%             plot(B_k(:,2),B_k(:,1),'.','color',[1 1 1],'MarkerSize',6);hold on;
%         end        
%         
%         disF=bwdist(bwmorph(ROI_Mask,'remove'));
%         mask = ROI_Mask;
%         mask = single(mask);
%         mask(mask==0)=-1;
%         disF = disF.*mask;
%         d = disF(sub2ind(size(disF),x,y))*PixelSize;
%         set(gcf, 'Name', fname, 'numbertitle', 'off');
% 
%         if ~isempty(d)
%             h44 = figure(44);
%             %         histogram(d);%axis square;
%             histfit(d,50,'kernel');xlim([0 inf]);
%             ylabel('Frequency');
%             xlabel('Distance to Cell Edge (\mum)');
%             title('Distribution of Junctions');
%             set(gcf, 'Name', fname, 'numbertitle', 'off');
%         end
% %Analysis 2
% 
% %Analysis 3
%         h55 = figure;
%         h=rose(pi*AnalysisInfo(:,1)/180,30);axis off;
%         set(h,'linewidth',3)
%         axis square;
%         daspect([1 1 1]);
%         title('Histogram of Filaments Orientations');
%         set(gcf, 'Name', fname, 'numbertitle', 'off');
                
%         allc = AnalysisInfo(:,4:5);
%         disF=bwdist(bwmorph(ROI_Mask,'remove'));
%         mask = ROI_Mask;
%         mask = single(mask);
%         mask(mask==0)=-1;
%         disF = disF.*mask;
%         d = disF(sub2ind(size(L),ceil(allc(:,1)),ceil(allc(:,2))))*PixelSize;
%         
%         finalmap = zeros(181,ceil(max(d)));
%         
%         h = waitbar(0,'generating colormap ...');
%         for i = 1:length(d)
%             waitbar(i/size(all_sorted_filament,3),h);
%             if ceil(d(i))>0
%                 finalmap(ceil(AnalysisInfo(i,1)+91), ceil(d(i))) = finalmap(ceil(AnalysisInfo(i,1)+91), ceil(d(i))) + 1;
%             end
%         end
%         close(h);
        
%         h2 = figure(2);
%         shape_factor = AnalysisInfo(:,2)./AnalysisInfo(:,3);
%         histogram(shape_factor,'Normalization','probability');
%         xlabel('tortuosity');
%         ylabel('probability');
%         axis([1 1.8 0 0.4]);
%         grid(gca,'on');        
%         
%                         
%         h66 = figure(66);
%         imagesc(finalmap);colormap(jet);colorbar;
%         set(gca,'ytick',[]);
%         xlabel('Distance to Cell Edge (\mum)');ylabel('-90 degrees to 90 degrees');
%         title('Distribution of Filament Orientations');
%         set(gcf, 'Name', fname, 'numbertitle', 'off');        
% %Analysis 3
% 
% %Analysis 4
        curR = round(FanR/2); % round(str2num(get(handles.FanR,'String'))/2);
        
        disF=bwdist(bwmorph(ROI_Mask,'remove'));
        mask = ROI_Mask;
        mask = single(mask);
        mask(mask==0)=-1;
        disF = disF.*mask;
        
        allcur = [];
        alldSDF = [];
        M = zeros(size(AllFragments));
        h = waitbar(0,'Calculating Curvatures ...');
        all_filament_curs = [];
        all_mean_cur = zeros(1,size(all_sorted_filament,3));
        for i = 1:size(all_sorted_filament,3)
            waitbar(i/size(all_sorted_filament,3),h);
            F = all_sorted_filament(:,:,i);
            F(find(F(:,1)==0),:) = [];
            Dlist = 0;
            for j = 2:size(F,1)
                Dlist = [Dlist  Dlist(j-1)+pdist([F(j-1,:); F(j,:)])];
            end
            NoCurlist = unique([find(Dlist<=curR)  find(Dlist>=(Dlist(end)-curR))]);
            if length(NoCurlist)~=size(F,1)
                for j = 1:size(F,1)
                    if isempty(find(NoCurlist==j))
                        x2 = F(j,1);
                        y2 = F(j,2);
                        D = Dlist - Dlist(j);
                        D1 = abs(D - curR);
                        idx1 = find(D1==min(D1));
                        D2 = abs(D - (-curR));
                        idx2 = find(D2==min(D2));
                        x1 = F(idx1,1);
                        y1 = F(idx1,2);
                        x3 = F(idx2,1);
                        y3 = F(idx2,2);
                        alldSDF = [alldSDF  disF(sub2ind(size(disF),x2,y2))*PixelSize];
                        cur = 2*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1)) ./sqrt(((x2-x1).^2+(y2-y1).^2)*((x3-x1).^2+(y3-y1).^2)*((x3-x2).^2+(y3-y2).^2));
                        allcur = [allcur  cur/PixelSize];
                        M(F(j,1),F(j,2)) = cur;
                        all_filament_curs(j,1,i) = x2;
                        all_filament_curs(j,2,i) = y2;
                        all_filament_curs(j,3,i) = cur/PixelSize;
                        all_filament_curs(j,4,i) = 1;% flag for calculation done for this point
                    end
                end
                k = find(all_filament_curs(:,4,i)==1);
                all_mean_cur(i) = mean(all_filament_curs(k,3,i));
            end
        end
        close(h);
%         M = M/PixelSize;
%         h88 = figure(88);
%         imagesc(M); colormap(jet);colorbar; title('Distribution of Curvatures');axis off;axis image;
%         set(gcf, 'Name', fname, 'numbertitle', 'off');
         
        All_Curs = [];
        All_Curs_4plot = [];
        disF=bwdist(bwmorph(ROI_Mask,'remove'));
        mask = ROI_Mask;
        mask = single(mask);
        mask(mask==0)=-1;
        disF = disF.*mask;
        %
        h = waitbar(0,'Calculating the Distribution of Curvatures ...');
        for m = 1:size(all_filament_curs,3)
            waitbar(m/size(all_filament_curs,3),h);
            Cur = all_filament_curs(:,:,m);
            F = all_sorted_filament;
            F(find(F(:,1)==0),:) = [];
            Cur(find(Cur(:,4)==0),:) = [];
            if isempty(Cur)
                All_Curs_4plot = [All_Curs_4plot;  disF(sub2ind(size(disF),round((F(1,1)+F(end,1))/2),round((F(end,2)+F(1,2))/2)))   0];
                All_Curs = [All_Curs, 0];
            else
                All_Curs_4plot = [All_Curs_4plot;  disF(sub2ind(size(disF),Cur(:,1),Cur(:,2)))   Cur(:,3)];
                All_Curs = [All_Curs, Cur(:,3)'];
            end
        end
        close(h);
        All_Curs_4plot(:,1) = round(All_Curs_4plot(:,1));
        temp_mean_cur = [];
        x1 = [];
        y1 = [];
        for i = min(All_Curs_4plot(:,1)): max(All_Curs_4plot(:,1))
            k = find(All_Curs_4plot(:,1)==i);
            temp_mean_cur = [temp_mean_cur  mean(All_Curs_4plot(k,2))];
            if length(k)>5
                x1 = [x1  i];
                y1 = [y1  mean(All_Curs_4plot(k,2))];
            end
        end

%         h77 = figure(77);
%         subplot(1,2,1);
%         histfit(all_mean_cur,50,'kernel');xlim([0 inf]);
%         xlabel('Curvature (unit: \mum^-^1)');ylabel('Frequency');
%         title('Distribution of Means of Filament Curvatures');
%         %             
%         subplot(1,2,2);
%         plot(x1*0.02, y1,  'r');hold on;axis([0 inf 0 inf]);
%         xlabel('Distance to Cell Edge (unit: \mum)');ylabel('Mean Curvature (unit: \mum^-^1)');
%         title('Distribution of Curvatures');
%         set(gcf, 'Name', fname, 'numbertitle', 'off');
        
% %Analysis 4

% THIS STRUCTURE IS TO BE EXTENDED..        
% orient = regionprops(maskI,'Orientation');
% ctrs = regionprops(maskI,'Centroid');
% ctrs = ctrs(1).Centroid;
% orient = orient(1).Orientation;
% TotalLength = d;
% EndToEndDist = pdist([tips(1,1),tips(1,2);tips(2,1),tips(2,2)]);
% 
% AnalysisInfo = [orient  TotalLength  EndToEndDist  ctrs(2)  ctrs(1)];        
% size(A)
% size(AnalysisInfo)
       
% tortuosity, 
% intensity_mean, 
% intensity_std, 
% local orientational jitter
% PATCH #

% PATCH #
patches_labels = bwlabel(ROI_Mask);

patches_area = zeros(1,max(patches_labels(:)));
for k=1:length(patches_area)
    patches_area(k) = sum(patches_labels==k,'all');
end

patches_inds = zeros(size(AnalysisInfo,1),1);
patches_areas = zeros(size(AnalysisInfo,1),1);
for k=1:size(AnalysisInfo,1)
     xc = round(AnalysisInfo(k,4));
     yc = round(AnalysisInfo(k,5));          
     % caution - rounding
     vic = patches_labels(xc-1:xc+1,yc-1:yc+1);
     ind = max(vic(:));     
     if 0==ind % should not happen
         ind=1;
     end
     patches_inds(k)=ind;
     patches_areas(k) = patches_area(ind); % pix
end

tortuosity = AnalysisInfo(:,2)./AnalysisInfo(:,3);

% filaments_data = zeros(size(AnalysisInfo,1),10);
% filaments_data(:,1:3)=AnalysisInfo(:,1:3);
% filaments_data(:,4)=AnalysisInfo(:,4)-R; % one needs original not padded image coordinates
% filaments_data(:,5)=AnalysisInfo(:,5)-R;
% filaments_data(:,6)=tortuosity;
% filaments_data(:,7)=arr_mean_filament_intensity;
% filaments_data(:,8)=arr_std_filament_intensity;
% filaments_data(:,9)=loj;
% filaments_data(:,10)=patches_inds;
% filaments_data(:,11)=patches_areas;

xls_filaments_data = [];
for k=1:size(AnalysisInfo,1)
    rec = {fname, ...
            patches_inds(k), ...
            AnalysisInfo(k,4)-R, ... % xc
            AnalysisInfo(k,5)-R, ... % yc
            all_mean_cur(k), ... % curvature
            AnalysisInfo(k,2)*PixelSize, ... % length           
            AnalysisInfo(k,3)*PixelSize, ... % ends distance
            tortuosity(k), ...
            arr_mean_filament_intensity(k), ...
            arr_std_filament_intensity(k), ...
            loj(k), ... % local orientational jitter
            patches_areas(k)}; % patch area 
            xls_filaments_data = [xls_filaments_data; rec];      
end
filaments_caption = {'filename','patch','xc [pix]','yc [pix]','curvature','length [um]','ends_dist [um]','tortuosity', ...
            'intensity mean','intensity std','local orientational jitter','patch area'};
xls_filaments_data = [filaments_caption; xls_filaments_data];
%xlswrite([save_dir filesep fname '_filaments_data'],xls_filaments_data);
xlwrite([save_dir filesep fname '_filaments_data.xls'],xls_filaments_data);

%%%%%%%%%%%%%%% junctions
% junctions analysis
% NewCrPts

if size(NewCrPts,1) >= 2

junctions = cell(max(patches_labels(:)),1);
for k=1:size(NewCrPts,1)
    xc=NewCrPts(k,1);
    yc=NewCrPts(k,2);
    lab = patches_labels(xc,yc);
    junctions{lab} = [junctions{lab}; [xc yc]];
end
%
max_length_pixels = 5/PixelSize; % 5 microns max distance between neighboring junctions  
junctions_quantification = cell(max(patches_labels(:)),1);
for k=1:max(patches_labels(:))         
    try
        junctions_quantification{k} = quantify_points_neighbours_and_density(junctions{k},max_length_pixels);
    catch
        disp(['no junctions in the patch ' num2str(k)]);
    end
end
%
xls_junctions_data = [];
for k=1:max(patches_labels(:))
    try
    d_k = junctions_quantification{k};
    for m=1:size(d_k,1)
        rec = {fname k d_k(m,1) d_k(m,2)/PixelSize^2};
        xls_junctions_data = [xls_junctions_data; rec];
    end
    catch
    end
end
%
if ~isempty(xls_junctions_data) 
    junctions_caption = {'filename','patch','nnghb','density [um^-2]'};
    xls_junctions_data = [junctions_caption; xls_junctions_data];
    %xlswrite([save_dir filesep fname '_junctions_data'],xls_junctions_data);
    xlwrite([save_dir filesep fname '_junctions_data.xls'],xls_junctions_data);
end
end
%%%%%%%%%%%%%%% junctions

%
fig = zeros(size(ROI_Mask,1),size(ROI_Mask,2),2,1,1);
z = patches_labels;
z(AllFragments>0) = 50;
z(flab>0) = z(flab>0) + 50;
if size(NewCrPts,1) ~= 0 
    x=NewCrPts(:,1);
    y=NewCrPts(:,2);
    z(sub2ind(size(L),x,y)) = 200;
end
fig(:,:,1,1,1) = intensity_img;
fig(:,:,2,1,1) = z;
fig = fig(R+1:H+R,R+1:W+R,:); % to keep same-size input and output images
bfsave(uint16(fig),[save_dir filesep fname '_image.ome.tif'],'Compression','LZW','BigTiff', true,'dimensionOrder','XYCZT');

% refactored analyses
                       
if exist('xls_junctions_data','var') && ~isempty(xls_junctions_data)
           datas = {xls_filaments_data xls_junctions_data};
           captions = {filaments_caption junctions_caption};
           table_names = {'filaments' 'junctions'};
else
           datas = {xls_filaments_data};
           captions = {filaments_caption};
           table_names = {'filaments'};    
end    

disp('analyze_SIFNE');          
end