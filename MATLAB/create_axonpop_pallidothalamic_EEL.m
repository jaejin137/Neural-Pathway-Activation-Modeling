%% Program that uses the XYZ coordinates of conturs to create axon streamlines
% Originally axon_populate_3D.m by Julia Slopsema 1/17/2019
% Adapted by Emily Lecy 11/2022
% Requires functions_populate_axons folder and functions_GP_emily to be added to path
% Use for pallidothalamic GPi output pathway (Ansa Lenticularis (AL) and Lenticular
% Fasciculus (LF) combination path).

%% ------------------------------------------------------------------------------------------------------------------------------
 % Step 0: Initialize variables -------------------------------------------------------------------------------------------------
clearvars; clc; close all;

user = 'jaejin'; 
subject_ID = 'Barb';
side = 'R';
subj= 'Barb';
nAxon = 10; %Number of axons to create (264 for AL and 536 for LF)
nSeedPoints=floor(nAxon*20); %Creates more seed points than axons to prevent strange streamlines - can adjust multiplier if needed
% geom_name = {'LF'}; %Name of contour path which drives output
geom_name = {'GPe_STN'}; %Name of contour path which drives output
axon_dia = 2; %Axon diameter 
rand_end = 0; %DO NOT CHANGE
nuclei_true=1; %DO NOT CHANGE

% cd(['C:\Users\' user '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data_Processing\STL Files'])
% [stl.GPe.Points,stl.GPe.ConnectivityList] = stlread1([subj '_' side '_GPe_DTIspace.stl']);
% [stl.GPi.Points,stl.GPi.ConnectivityList] = stlread1([subj '_' side '_GPi_DTIspace.stl']);
cd('C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\Data_Processing\STL')
[stl.GPe.Points,stl.GPe.ConnectivityList] = stlread1('GPe_R.stl');
[stl.GPi.Points,stl.GPi.ConnectivityList] = stlread1('GPi_R.stl');
% [stl.GPe.ConnectivityList, stl.GPe.Points] = stlread('GPe_R.stl');
% [stl.GPi.ConnectivityList, stl.GPi.Points] = stlread('GPi_R.stl');
stl.GPeTri = zeros(length(stl.GPe.ConnectivityList),9);
stl.GPiTri = zeros(length(stl.GPi.ConnectivityList),9);

% inFolder = ['C:\Users\' user '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data_Processing\'];
% outFolder = ['C:\Users\' user '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data_Processing\Axons\'];
inFolder = 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\Data_Processing\';
outFolder = 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\Data_Processing\Axons\';
    if ~exist(outFolder,'dir')
    mkdir(outFolder)
    end
cell_body = 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\Data_Processing\Axons\Cell_body';
    if ~exist(cell_body,'dir')
    mkdir(cell_body)
    end

    
%for ii=1:length(geom_name)
ii = 1;
% cd(['C:\Users\' user '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data_Processing\Contours\'])  
cd('C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\Data_Processing\Contours')

in_foldername = [inFolder 'Contours\' side '_' geom_name{ii} '\' ];
    
%% ------------------------------------------------------------------------------------------------------------------------------    
 % Step 1: import points, make contours, and plot them --------------------------------------------------------------------------
 % preallocate memory
    for i = 1:length(stl.GPiTri)
        stl.GPiTri(i,:) = [stl.GPi.Points(stl.GPi.ConnectivityList(i,1),:),...
        stl.GPi.Points(stl.GPi.ConnectivityList(i,2),:),...
        stl.GPi.Points(stl.GPi.ConnectivityList(i,3),:)];
%         stl.GPiTri(i,:) = [stl.GPi.Points.Points(stl.GPi.Points.ConnectivityList(i,1),:),...
%         stl.GPi.Points.Points(stl.GPi.Points.ConnectivityList(i,2),:),...
%         stl.GPi.Points.Points(stl.GPi.Points.ConnectivityList(i,3),:)];
    end
    
    bb.GPi.x = [min(stl.GPi.Points(:,1)),max(stl.GPi.Points(:,1))];
    bb.GPi.y = [min(stl.GPi.Points(:,2)),max(stl.GPi.Points(:,2))];
    bb.GPi.z = [min(stl.GPi.Points(:,3)),max(stl.GPi.Points(:,3))];
%     bb.GPi.x = [min(stl.GPi.Points.Points(:,1)),max(stl.GPi.Points.Points(:,1))];
%     bb.GPi.y = [min(stl.GPi.Points.Points(:,2)),max(stl.GPi.Points.Points(:,2))];
%     bb.GPi.z = [min(stl.GPi.Points.Points(:,3)),max(stl.GPi.Points.Points(:,3))];
    

    h = waitbar(0,'Please wait...');
    clear GPi
    num.GPi = nSeedPoints;
    sz=[5 (nSeedPoints)]; GPi=zeros(sz);
%     for i = 1:num.GPi
%             waitbar(i/num.GPi,h,['Checking GPi cell: ',num2str(i)]);
%             temp(1) = rand(1)*(bb.GPi.x(2)-bb.GPi.x(1)) + bb.GPi.x(1);
%             temp(2) = rand(1)*(bb.GPi.y(2)-bb.GPi.y(1)) + bb.GPi.y(1);
%             temp(3) = rand(1)*(bb.GPi.z(2)-bb.GPi.z(1)) + bb.GPi.z(1);
% 
%             % Check if cell is within the GPi
% %             while jordancurve(stl.GPiTri,temp)==0
%             while true
%                 temp(1) = rand(1)*(bb.GPi.x(2)-bb.GPi.x(1)) + bb.GPi.x(1);
%                 temp(2) = rand(1)*(bb.GPi.y(2)-bb.GPi.y(1)) + bb.GPi.y(1);
%                 temp(3) = rand(1)*(bb.GPi.z(2)-bb.GPi.z(1)) + bb.GPi.z(1);
%             end
%           GPi(1:3,i)=temp(1:3);
%     end
    GPi_temp = zeros(3, num.GPi);
    pidx_max = num.GPi;
    for i = 1:pidx_max
%             waitbar(i/pidx_max,h,['Checking GPi cell: ',num2str(i)]);
            temp = [rand(1)*(bb.GPi.x(2)-bb.GPi.x(1)) + bb.GPi.x(1),...
                rand(1)*(bb.GPi.y(2)-bb.GPi.y(1)) + bb.GPi.y(1),...
                rand(1)*(bb.GPi.z(2)-bb.GPi.z(1)) + bb.GPi.z(1)];

            % Check if cell is within the GPi
            while jordancurve(stl.GPiTri,temp)==0
%             while is_inside(stl.GPiTri,temp)
                temp(1) = rand(1)*(bb.GPi.x(2)-bb.GPi.x(1)) + bb.GPi.x(1);
                temp(2) = rand(1)*(bb.GPi.y(2)-bb.GPi.y(1)) + bb.GPi.y(1);
                temp(3) = rand(1)*(bb.GPi.z(2)-bb.GPi.z(1)) + bb.GPi.z(1);
            end
          GPi_temp(:,i) = temp(1:3);
    end
    GPi(1:3,:) = GPi_temp;
 close(h)   
    
 %% ------------------------------------------------------------------------------------------------------------------------------
  % STEP 2: assign axon type to each GPicell abnd to DV location. To do this, make a random matrix of numnbers -------------------
  %and use this random assay to assign axon types and DV projections
  %At the end of this, you should have a GPi grid where rows 1-3 are the
  %[xyz] points. Row 4 defines if the point will project out ventrally or dorsally.

 DV_choice=reshape(randperm(nSeedPoints),1,nSeedPoints);

% Add dorsal ventral specification to all axons
for w=1:length(DV_choice)
    ww=DV_choice(1,w);
    s = (-1)^ww; %if s=1, number is even (assign D, else if, assign odd)
        if s==1
            GPi(4,w) = 1;
        else
             GPi(4,w) = 2;
        end
end

%% ------------------------------------------------------------------------------------------------------------------------------ 
 % Step 3: import points, make contours, and plot them --------------------------------------------------------------------------
    % Goal 1: import points, make contours, and plot them
    nContour = 1; % number of contours
    nContour_pts = 100; % number of points per contour
    
    % preallocate mem
    contour = zeros(nContour_pts+1,3); % nContour_pts+1 b/c one is added to close the contour
    contour_pts = zeros(nSeedPoints,3);
    
   % h1=figure; hold on;
        
        %edits made by annie 10/19/2020 to account for contours with less than or more than 100 pts
        cont = load([in_foldername '/1.txt']); % load contour points
        if length(cont) >=100
%             contour(1:nContour_pts,:)=cont(1:100,:); % load contour points
            contour(1:nContour_pts,:)=cont(1:100,1:3); % load contour points
        elseif length(cont) == 99
%             contour(1:nContour_pts,:)=[cont(1:99,:);cont(99,:)];
            contour(1:nContour_pts,:)=[cont(1:99,1:3);cont(99,1:3)];

        end
        
        normal_of_contour_plane = affine_fit(contour(1:nContour_pts,:));
        contour(nContour_pts+1,:)=contour(1,:,1);                              % close the contours
        plot3(contour(:,1),contour(:,2),contour(:,3)); hold on            % plot the contour

        % Use axon_1_populate_contour_xy for countours that run along the z-axis
        % and axon_1_populate_contour_xz for contours that run along the x- or
        % y-axis. Either can be used on non-vertical or non-horizonatal running
        % contours.
        CONTOUR_DIRECTION_THRESHOLD = 0.7;
        if (abs(dot(normal_of_contour_plane,[0,0,1])) > CONTOUR_DIRECTION_THRESHOLD)
            contour_pts(:,:,1)=axon_1_populate_contour_xy(contour(:,:,1),nSeedPoints); % populate contour with random points
        elseif (abs(dot(normal_of_contour_plane,[0,1,0])) > CONTOUR_DIRECTION_THRESHOLD)
            contour_pts(:,:,1)=axon_1_populate_contour_xz(contour(:,:,1),nSeedPoints); % populate contour with random points
        else
            contour_pts(:,:,1)=axon_1_populate_contour_yz(contour(:,:,1),nSeedPoints); % populate contour with random points
        end      
        %plot3(contour_pts(:,1,1),contour_pts(:,2,1),contour_pts(:,3,1),'r.');  % plot contours
%     close all;
  
%% ------------------------------------------------------------------------------------------------------------------------------
 % STEP 4: define is contour points are in the ventral or dorsal aspect of the contour ------------------------------------------
    bb.contour.x = [min(contour(:,1)),max(contour(:,1))];
    bb.contour.y = [min(contour(:,2)),max(contour(:,2))];
    bb.contour.z = [min(contour(:,3)),max(contour(:,3))];
  
    DVdistance=bb.contour.z(:,2)-bb.contour.z(:,1); %finds center of DV regions
    DVgoal=DVdistance/2; DVgoal=bb.contour.z(:,1)+DVgoal;

    %define is z point is on the dorsal or ventral side, it then creates
    %row 4 in the contour_pts variable which assigns a point as dorsal (2)
    %or ventral (1). Row 5 is just there so matrices match-IGNORE IT
contour_pts=contour_pts';

for pp=1:nSeedPoints
    Z=contour_pts(3,pp);
   if Z>DVgoal
       contour_pts(4,pp)=1;
   elseif Z<DVgoal 
       contour_pts(4,pp)=2;
   end
end 
contour_pts(5,:)=ones;


%combine both dorsal and ventral side, the first n (dorsalnum) are dorsal
%points, and the next n (ventralnum) is the second half of th

    final_pts = zeros(5,nSeedPoints,2);
    final_pts(:,:,1)=GPi;
    final_pts(:,:,2)=contour_pts;
%     close all;

%% ------------------------------------------------------------------------------------------------------------------------------
 % Step 5: Connect the random points across contours ----------------------------------------------------------------------------
 % flip for TRIAL
[fiber_pts,final_connection_pts1,final_connection_pts2,final_connection_pts3] = axon_4_connect_curved_Pallidothalamic (final_pts,nAxon,nSeedPoints);
axonpop = axon_4_axonpop_struct_fiberD(fiber_pts, 2);

%% ------------------------------------------------------------------------------------------------------------------------------
 % STEP 6: Make figures to save  ----------------------------------------------- ------------------------------------------------
grey=[.5 .5 .5];
h1=figure;
% nuc = stlread2(['C:\Users\', user, '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data_Processing\STL Files\' subj '_' side '_GPi_DTIspace.stl']);
nuc = stlread2('C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\Data_Processing\STL\GPi_R.stl');

for i=1:nAxon
   plot3(fiber_pts{i}(:,1),fiber_pts{i}(:,2),fiber_pts{i}(:,3))
        for iContour=1:nContour
        plot3(contour(:,1,iContour),contour(:,2,iContour),contour(:,3,iContour)); hold on  
        hold on;
        end 
end  
    hold on;
    patch('Faces',nuc.faces,'Vertices',nuc.vertices,'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5],'FaceAlpha',0.2,'EdgeAlpha',0.2);
    hold on;

    cd(outFolder)
    savefig(h1,[subj '_' side '_GPiOutput'])
    

%% ------------------------------------------------------------------------------------------------------------------------------
 % Step 7: make cell body within axonpop matrix for use later--------------------------------------------------------------------
 % make cell body matrix for all cells
cellbod1=final_connection_pts1(1:3,:,1)';
cellbod2=final_connection_pts2(1:3,:,1)';
cellbod3=final_connection_pts3(1:3,:,1)';
cellbod_final=cat(1,cellbod1,cellbod2,cellbod3); 

 % add to axonpop struct
GPix=num2cell(cellbod_final(:,1));
    [axonpop(:).cellbod_GPix]= GPix{:};
GPiy=num2cell(cellbod_final(:,2));
    [axonpop(:).cellbod_GPiy]= GPiy{:};
GPiz=num2cell(cellbod_final(:,3));
    [axonpop(:).cellbod_GPiz]= GPiz{:};

%% ------------------------------------------------------------------------------------------------------------------------------
 % Step 8: clean and save workspace, figures, and variable final_connection_points ----------------------------------------------- -----------------------------------------------
 clear each; clear eachsad; clear i; clear ii; clear pp; clear r; clear rr; clear s; clear temp; clear sz; clear w; clear ww; clear Z; clear CONTOUR_DIRECTION_THRESHOLD;
 cd(outFolder)
 %savefig([(subj '_' side '_' geom_name '_fig')) 
 save([subj '_' side '_GPi_Output_Workspace'])
 save([subj '_' side '_GPi_Output_axonpop'], 'axonpop')
 cd(cell_body)
 save([subj '_' side '__GPi_Output_cellbod1'], 'final_connection_pts1')
 save([subj '_' side '__GPi_Output_cellbod2'], 'final_connection_pts2')
 save([subj '_' side '__GPi_Output_cellbod3'], 'final_connection_pts3')
 
%% ------------------------------------------------------------------------------------------------------------------------------
 % Step 9: delete intersecting axons---------------------------------------------------------------------------------------------
geom_name= {'PutGPeGPi'};

for i = 1:length(geom_name)
%     cd(['\Users\', user, '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data_Processing\STL Files']) 
    cd('C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\Data_Processing\STL') 

    % Temporarily disabled - JL
%     axons_that_intersect = []; 
%     surface_struct_DWI_space  = stlread2([subj '_' side '_encap_COMSOL.stl']);
%     axons_that_intersect = get_axons_that_intersect(axonpop,surface_struct_DWI_space,false); %if axon numbers don't make sense, change input from false to true
% 
%     axonpop(axons_that_intersect)=[];
%     outFolder = ['C:\Users\' user '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data Processing\Axons\Cleaned\'];
    outFolder = 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\Data_Processing\Axons\Cleaned\';

    mkdir(outFolder)
    cd(outFolder)
%     save(strcat(subj, '_', side, '_', geom_name), 'axonpop')
    save(strcat(subj, '_', side, '_', geom_name{i}), 'axonpop')

end
