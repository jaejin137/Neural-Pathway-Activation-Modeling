%% Program uses the XYZ coordinates of nuclei to create curved axon streamlines
% Emily Lecy 11/2022, functions in BOX
% Use for internal pallidal pathways (GPeGPiSTNSN, GPeSTNSN, PutGPe, PutGPeGPi)
% Takes a while becuase you need to populate a lot of cells (to connect
% successfully, just wait it out)

% If your axons are not evenly populating all nuclei, you can change the
% nSeedpoints to be largeer or smaller until you get your desired results.
% You can also change nAxons if needed, but then all axon ratios should be changes.  

%% ---------------------------------------------------------------------------------------------------
% Step:0 Initialize variables-------------------------------------------------------------------------
clearvars;
close all;
clc;

% user = 'elecy'; 
% subject_ID = 'UD1028_PD104';
% side = 'L';
% subj= 'UD1028';
user = 'jaejin'; 
subject_ID = 'Barb';
side = 'R';
subj= 'Barb';
% nAxon = 784; %# of axons (GPeGPi 750, GPeSTN 1250, PutGPe 784, PutGPeGPi 816)
nAxon = 1250; %# of axons (GPeGPi 750, GPeSTN 1250, PutGPe 784, PutGPeGPi 816)
nSeedPoints=floor(nAxon*25); %Creates more seed points than axons to prevent strange streamlines - can adjust multiplier if needed
geom_name = 'GPeSTN'; %Name (abbreviation) of the pathway to map
rand_end = 0; %DO NOT CHANGE
axon_dia = 2; %Axon diameter 
nuclei_true=1; % DO NOT CHANGE
nContour=2; %(number of nuclei connected; GPeGPi 3, GPeSTN 3, PutGPe 2, PutGPeGPi 4)

% cd(['C:\Users\' user '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data Processing\STL Files'])
% [stl.GPe.Points,stl.GPe.ConnectivityList] = stlread1([subj '_' side '_GPe_DTIspace.stl']);
% [stl.GPi.Points,stl.GPi.ConnectivityList] = stlread1([subj '_' side '_GPi_DTIspace.stl']);
% [stl.STN.Points,stl.STN.ConnectivityList] = stlread1([subj '_' side '_STN_DTIspace.stl']);
% [stl.SN.Points,stl.SN.ConnectivityList] = stlread1([subj '_' side '_SN_DTIspace.stl']);
% [stl.GPiRing.Points,stl.GPiRing.ConnectivityList] = stlread1([subj '_' side '_GPiRing.stl']);
% [stl.Put.Points,stl.Put.ConnectivityList] = stlread1([subj '_' side '_put_DTIspace.stl']);
cd('C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\Data_Processing\STL')
[stl.GPe.Points,stl.GPe.ConnectivityList] = stlread1('GPe_R.stl');
[stl.GPi.Points,stl.GPi.ConnectivityList] = stlread1('GPi_R.stl');
[stl.STN.Points,stl.STN.ConnectivityList] = stlread1('STN_R.stl');
[stl.GPiRing.Points,stl.GPiRing.ConnectivityList] = stlread1('R_GPiRing_new2.stl');

stl.GPeTri = zeros(length(stl.GPe.ConnectivityList),9);
stl.GPiTri = zeros(length(stl.GPi.ConnectivityList),9);
stl.STNTri = zeros(length(stl.STN.ConnectivityList),9);
%stl.SNTri = zeros(length(stl.SN.ConnectivityList),9);
stl.GPiRingTri = zeros(length(stl.GPiRing.ConnectivityList),9);
%stl.PutTri = zeros(length(stl.Put.ConnectivityList),9);

% inFolder = ['C:\Users\' user '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data Processing\'];
% outFolder = ['C:\Users\' user '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data Processing\Axons\'];
inFolder = 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\Data_Processing';
outFolder = 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\Data_Processing\Axons';
    if ~exist(outFolder,'dir')
    mkdir(outFolder)
    end

%% ----------------------------------------------------------------------------------------------------------------------------
% Step 1: import points, make contours, and plot them -------------------------------------------------------------------------
% preallocate memory

if strcmp (geom_name, 'GPeGPiSTNSN') ==1  %populate cells in GPe if applicable
    final_pts = zeros(3,nSeedPoints,(nContour));
        for i = 1:length(stl.GPeTri)
        stl.GPeTri(i,:) = [stl.GPe.Points(stl.GPe.ConnectivityList(i,1),:),...
        stl.GPe.Points(stl.GPe.ConnectivityList(i,2),:),...
        stl.GPe.Points(stl.GPe.ConnectivityList(i,3),:)];
        end
       for i = 1:length(stl.GPiTri)
        stl.GPiTri(i,:) = [stl.GPi.Points(stl.GPi.ConnectivityList(i,1),:),...
        stl.GPi.Points(stl.GPi.ConnectivityList(i,2),:),...
        stl.GPi.Points(stl.GPi.ConnectivityList(i,3),:)];
       end
       for i = 1:length(stl.STNTri)
        stl.STNTri(i,:) = [stl.STN.Points(stl.STN.ConnectivityList(i,1),:),...
        stl.STN.Points(stl.STN.ConnectivityList(i,2),:),...
        stl.STN.Points(stl.STN.ConnectivityList(i,3),:)];
       end
       for i = 1:length(stl.SNTri)
        stl.SNTri(i,:) = [stl.SN.Points(stl.SN.ConnectivityList(i,1),:),...
        stl.SN.Points(stl.SN.ConnectivityList(i,2),:),...
        stl.SN.Points(stl.SN.ConnectivityList(i,3),:)];
       end
end

if strcmp (geom_name, 'GPeSTNSN') ==1   %populate cells in GPe if applicable
    final_pts = zeros(3,nSeedPoints,(nContour));
    for i = 1:length(stl.GPeTri)
        stl.GPeTri(i,:) = [stl.GPe.Points(stl.GPe.ConnectivityList(i,1),:),...
        stl.GPe.Points(stl.GPe.ConnectivityList(i,2),:),...
        stl.GPe.Points(stl.GPe.ConnectivityList(i,3),:)];
    end
    for i = 1:length(stl.STNTri)
        stl.STNTri(i,:) = [stl.STN.Points(stl.STN.ConnectivityList(i,1),:),...
        stl.STN.Points(stl.STN.ConnectivityList(i,2),:),...
        stl.STN.Points(stl.STN.ConnectivityList(i,3),:)];
    end  
    for i = 1:length(stl.GPiRingTri)
        stl.GPiRingTri(i,:) = [stl.GPiRing.Points(stl.GPiRing.ConnectivityList(i,1),:),...
        stl.GPiRing.Points(stl.GPiRing.ConnectivityList(i,2),:),...
        stl.GPiRing.Points(stl.GPiRing.ConnectivityList(i,3),:)];
    end
end

if strcmp (geom_name, 'PutGPe') ==1   %populate cells in GPe if applicable
    final_pts = zeros(3,nSeedPoints,(nContour));
    for i = 1:length(stl.PutTri)
        stl.PutTri(i,:) = [stl.Put.Points(stl.Put.ConnectivityList(i,1),:),...
        stl.Put.Points(stl.Put.ConnectivityList(i,2),:),...
        stl.Put.Points(stl.Put.ConnectivityList(i,3),:)];
    end
    for i = 1:length(stl.GPeTri)
        stl.GPeTri(i,:) = [stl.GPe.Points(stl.GPe.ConnectivityList(i,1),:),...
        stl.GPe.Points(stl.GPe.ConnectivityList(i,2),:),...
        stl.GPe.Points(stl.GPe.ConnectivityList(i,3),:)];
    end 
end

if strcmp (geom_name, 'PutGPeGPi') ==1   %populate cells in GPe if applicable
    final_pts = zeros(3,nSeedPoints,(nContour));
    for i = 1:length(stl.PutTri)
        stl.PutTri(i,:) = [stl.Put.Points(stl.Put.ConnectivityList(i,1),:),...
        stl.Put.Points(stl.Put.ConnectivityList(i,2),:),...
        stl.Put.Points(stl.Put.ConnectivityList(i,3),:)];
    end
    for i = 1:length(stl.GPeTri)
        stl.GPeTri(i,:) = [stl.GPe.Points(stl.GPe.ConnectivityList(i,1),:),...
        stl.GPe.Points(stl.GPe.ConnectivityList(i,2),:),...
        stl.GPe.Points(stl.GPe.ConnectivityList(i,3),:)];
    end
     for i = 1:length(stl.GPiTri)
        stl.GPiTri(i,:) = [stl.GPi.Points(stl.GPi.ConnectivityList(i,1),:),...
        stl.GPi.Points(stl.GPi.ConnectivityList(i,2),:),...
        stl.GPi.Points(stl.GPi.ConnectivityList(i,3),:)];
     end
    for i = 1:length(stl.STNTri)
        stl.STNTri(i,:) = [stl.STN.Points(stl.STN.ConnectivityList(i,1),:),...
        stl.STN.Points(stl.STN.ConnectivityList(i,2),:),...
        stl.STN.Points(stl.STN.ConnectivityList(i,3),:)];
    end
    for i = 1:length(stl.SNTri)
        stl.SNTri(i,:) = [stl.SN.Points(stl.SN.ConnectivityList(i,1),:),...
        stl.SN.Points(stl.SN.ConnectivityList(i,2),:),...
        stl.SN.Points(stl.SN.ConnectivityList(i,3),:)];
    end  
end

% Added by JL
if strcmp (geom_name, 'GPeSTN') ==1   %populate cells in GPe if applicable
    final_pts = zeros(3,nSeedPoints,(nContour));
    for i = 1:length(stl.GPeTri)
        stl.GPeTri(i,:) = [stl.GPe.Points(stl.GPe.ConnectivityList(i,1),:),...
        stl.GPe.Points(stl.GPe.ConnectivityList(i,2),:),...
        stl.GPe.Points(stl.GPe.ConnectivityList(i,3),:)];
    end
    for i = 1:length(stl.STNTri)
        stl.STNTri(i,:) = [stl.STN.Points(stl.STN.ConnectivityList(i,1),:),...
        stl.STN.Points(stl.STN.ConnectivityList(i,2),:),...
        stl.STN.Points(stl.STN.ConnectivityList(i,3),:)];
    end  
    for i = 1:length(stl.GPiRingTri)
        stl.GPiRingTri(i,:) = [stl.GPiRing.Points(stl.GPiRing.ConnectivityList(i,1),:),...
        stl.GPiRing.Points(stl.GPiRing.ConnectivityList(i,2),:),...
        stl.GPiRing.Points(stl.GPiRing.ConnectivityList(i,3),:)];
    end
end
%% ----------------------------------------------------------------------------------------------------------------------------
% Step 2: randomly populate points on sides of each GPe, GPi, STN, and SN------------------------------------------------------
if strcmp (geom_name, 'GPeGPiSTNSN') ==1   %populate cells in GPe if applicable
[GPe,GPi,STN] = PopGPeGPi(stl,nSeedPoints);
end

if strcmp (geom_name, 'GPeSTNSN') ==1   %populate cells in GPe if applicable
[GPe,GPiRing,STN] = PopGPeSTN(stl,nSeedPoints);
end

if strcmp (geom_name, 'PutGPe') ==1   %populate cells in GPe if applicable
[Put,GPe] = PopPutGPe(stl,nSeedPoints);
end

if strcmp (geom_name, 'PutGPeGPi') ==1   %populate cells in GPe if applicable
[Put,GPe,GPi,STN] = PopPutGPeGPi(stl,nSeedPoints);
end

% Added by JL
if strcmp (geom_name, 'GPeSTN') ==1   %populate cells in GPe if applicable
[GPe,GPiRing,STN] = PopGPeSTN(stl,nSeedPoints);
end

%% ----------------------------------------------------------------------------------------------------------------------------
% Step 3: make points into final variable named final_pts (this will be input to connect_contour function----------------------
if strcmp (geom_name, 'GPeGPiSTNSN') ==1  
    final_pts(:,:,1)=GPe;
    final_pts(:,:,2)=GPi;
    final_pts(:,:,3)=STN;
end 

if strcmp (geom_name, 'GPeSTNSN') ==1  
    final_pts(:,:,1)=GPe;
    final_pts(:,:,2)=GPiRing;
    final_pts(:,:,3)=STN;
end

if strcmp (geom_name, 'PutGPe') ==1  
    final_pts(:,:,1)=Put;
    final_pts(:,:,2)=GPe;
end

if strcmp (geom_name, 'PutGPeGPi') ==1  
    final_pts(:,:,1)=Put;
    final_pts(:,:,2)=GPe;
    final_pts(:,:,3)=GPi;
    final_pts(:,:,4)=STN;
end

% Added by JL
if strcmp (geom_name, 'GPeSTN') ==1  
    final_pts(:,:,1)=GPe;
    final_pts(:,:,2)=GPiRing;
    final_pts(:,:,3)=STN;
end
%% ----------------------------------------------------------------------------------------------------------------------------
% Step 4: Connect the random points across contours and create axonpop strucutre-----------------------------------------------
% If this errors, just rerun the axon_4_connect function, this function is
% super picky and it is all about the randomness in your code, so by chance
% it sometimes fails- just rerun until it goes- there is no real error

if strcmp (geom_name, 'GPeGPiSTNSN') ==1  
    [fiber_pts,final_connection_pts1,Ax_choice] = axon_4_connect_curved_GPeGPiSTN (final_pts,nAxon,nSeedPoints);
    axonpop = axon_4_axonpop_struct_fiberD(fiber_pts, axon_dia); 
end

if strcmp (geom_name, 'GPeSTNSN') ==1  
    [fiber_pts,final_connection_pts1] = axon_4_connect_curved_GPeSTN (final_pts,nAxon,nSeedPoints);
    axonpop = axon_4_axonpop_struct_fiberD(fiber_pts, axon_dia); 
end

if strcmp (geom_name, 'PutGPe') ==1  
    [fiber_pts,final_connection_pts1] = axon_4_connect_curved_PutGPe (final_pts,nAxon,nSeedPoints);
    axonpop = axon_4_axonpop_struct_fiberD(fiber_pts, axon_dia); 
end

if strcmp (geom_name, 'PutGPeGPi') ==1  
    [fiber_pts,final_connection_pts1] = axon_4_connect_curved_PutGPeGPi(final_pts,nAxon,nSeedPoints);
    axonpop = axon_4_axonpop_struct_fiberD(fiber_pts, axon_dia); 
end

% Added by JL
if strcmp (geom_name, 'GPeSTN') ==1  
    [fiber_pts,final_connection_pts1] = axon_4_connect_curved_GPeSTN(final_pts,nAxon,nSeedPoints);
    axonpop = axon_4_axonpop_struct_fiberD(fiber_pts, axon_dia); 
end
%% ----------------------------------------------------------------------------------------------------------------------------
% Step 5: Create figures for sanity checks and save----------------------------------------------------------------------------
grey=[.7 .7 .7];
lightgrey=[.8 .8 .8];
moregrey=[.9 .9 .9];

if strcmp (geom_name, 'GPeGPiSTNSN') ==1  
    h1=figure; 
    g = trimesh(stl.GPe.ConnectivityList, stl.GPe.Points(:,1), stl.GPe.Points(:,2), stl.GPe.Points(:,3),'EdgeColor',moregrey);
    g.FaceAlpha = 0;
    hold on;
    
    g = trimesh(stl.GPi.ConnectivityList, stl.GPi.Points(:,1), stl.GPi.Points(:,2), stl.GPi.Points(:,3),'EdgeColor',moregrey);
    g.FaceAlpha = 0;
    hold on;
    
    g = trimesh(stl.STN.ConnectivityList, stl.STN.Points(:,1), stl.STN.Points(:,2), stl.STN.Points(:,3),'EdgeColor',moregrey);
    g.FaceAlpha = 0;
    hold on;

    for iAxon=1:nAxon
        plot3(fiber_pts{iAxon}(:,1),fiber_pts{iAxon}(:,2),fiber_pts{iAxon}(:,3),'LineWidth',1)
        hold on;
    end
    cd(outFolder)
    savefig(h1,[subj '_' side '_GPeGPiSTNSNAxon'])
end

if strcmp (geom_name, 'GPeSTNSN') ==1  
    h1=figure; 
    g = trimesh(stl.GPe.ConnectivityList, stl.GPe.Points(:,1), stl.GPe.Points(:,2), stl.GPe.Points(:,3),'EdgeColor',grey);
    g.FaceAlpha = 0;
    hold on;
    
    g = trimesh(stl.GPi.ConnectivityList, stl.GPi.Points(:,1), stl.GPi.Points(:,2), stl.GPi.Points(:,3),'EdgeColor',lightgrey);
    g.FaceAlpha = 0;
    hold on;
    
    g = trimesh(stl.STN.ConnectivityList, stl.STN.Points(:,1), stl.STN.Points(:,2), stl.STN.Points(:,3),'EdgeColor',moregrey);
    g.FaceAlpha = 0;
    hold on;

    for iAxon=1:nAxon
        plot3(fiber_pts{iAxon}(:,1),fiber_pts{iAxon}(:,2),fiber_pts{iAxon}(:,3))
        hold on;
    end
    cd(outFolder)
    savefig(h1,[subj '_' side '_GPeSTNSNAxon'])
end

if strcmp (geom_name, 'PutGPe') ==1  
h1=figure; 
g = trimesh(stl.Put.ConnectivityList, stl.Put.Points(:,1), stl.Put.Points(:,2), stl.Put.Points(:,3),'EdgeColor',lightgrey);
g.FaceAlpha = 0;
hold on;

g = trimesh(stl.GPe.ConnectivityList, stl.GPe.Points(:,1), stl.GPe.Points(:,2), stl.GPe.Points(:,3),'EdgeColor',grey);
g.FaceAlpha = 0;
hold on;

  for iAxon=1:nAxon
        plot3(fiber_pts{iAxon}(:,1),fiber_pts{iAxon}(:,2),fiber_pts{iAxon}(:,3))
        hold on;
  end
 cd(outFolder)
 savefig(h1,[subj '_' side '_PutGPeAxon'])
end

if strcmp (geom_name, 'PutGPeGPi') ==1  
h1=figure; 
g = trimesh(stl.Put.ConnectivityList, stl.Put.Points(:,1), stl.Put.Points(:,2), stl.Put.Points(:,3),'EdgeColor',lightgrey);
g.FaceAlpha = 0;
hold on;

g = trimesh(stl.GPe.ConnectivityList, stl.GPe.Points(:,1), stl.GPe.Points(:,2), stl.GPe.Points(:,3),'EdgeColor',grey);
g.FaceAlpha = 0;
hold on;

g = trimesh(stl.GPi.ConnectivityList, stl.GPi.Points(:,1), stl.GPi.Points(:,2), stl.GPi.Points(:,3),'EdgeColor',lightgrey);
g.FaceAlpha = 0;
hold on;

g = trimesh(stl.STN.ConnectivityList, stl.STN.Points(:,1), stl.STN.Points(:,2), stl.STN.Points(:,3),'EdgeColor',moregrey);
g.FaceAlpha = 0;
hold on;

  for iAxon=1:nAxon
        plot3(fiber_pts{iAxon}(:,1),fiber_pts{iAxon}(:,2),fiber_pts{iAxon}(:,3))
        hold on;
  end
 cd(outFolder)
 savefig(h1,[subj '_' side '_PutGPeGPiAxon'])
end

% Added by JL
if strcmp (geom_name, 'GPeSTN') ==1  
    h1=figure; 
%     g = trimesh(stl.GPe.ConnectivityList, stl.GPe.Points(:,1), stl.GPe.Points(:,2), stl.GPe.Points(:,3),'EdgeColor',grey);
    g = trisurf(stl.GPe.ConnectivityList, stl.GPe.Points(:,1), stl.GPe.Points(:,2), stl.GPe.Points(:,3),'FaceColor','b','EdgeColor',grey);
    g.FaceAlpha = 1;
    hold on;
    
%     g = trimesh(stl.GPi.ConnectivityList, stl.GPi.Points(:,1), stl.GPi.Points(:,2), stl.GPi.Points(:,3),'EdgeColor',lightgrey);
    g = trisurf(stl.GPi.ConnectivityList, stl.GPi.Points(:,1), stl.GPi.Points(:,2), stl.GPi.Points(:,3),'FaceColor','c','EdgeColor',grey);
    g.FaceAlpha = 1;
    hold on;
    
%     g = trimesh(stl.STN.ConnectivityList, stl.STN.Points(:,1), stl.STN.Points(:,2), stl.STN.Points(:,3),'EdgeColor',moregrey);
    g = trisurf(stl.STN.ConnectivityList, stl.STN.Points(:,1), stl.STN.Points(:,2), stl.STN.Points(:,3),'FaceColor','m','EdgeColor',grey);
    g.FaceAlpha = 1;
    hold on;

%     g = trimesh(stl.GPiRing.ConnectivityList, stl.GPiRing.Points(:,1), stl.GPiRing.Points(:,2), stl.GPiRing.Points(:,3),'EdgeColor',moregrey);
    g = trisurf(stl.GPiRing.ConnectivityList, stl.GPiRing.Points(:,1), stl.GPiRing.Points(:,2), stl.GPiRing.Points(:,3),'FaceColor','y','EdgeColor',grey);
    g.FaceAlpha = .3;
    hold on;

    if nAxon <= 64
        for iAxon=1:nAxon
            plot3(fiber_pts{iAxon}(:,1),fiber_pts{iAxon}(:,2),fiber_pts{iAxon}(:,3),LineWidth=1)
            hold on;
        end
    else
        for iAxon=1:nAxon
            plot3(fiber_pts{iAxon}(:,1),fiber_pts{iAxon}(:,2),fiber_pts{iAxon}(:,3))
            hold on;
        end
    end

    cd(outFolder)
    savefig(h1,[subj '_' side '_GPeSTNAxon'])
end

%% ------------------------------------------------------------------------------------------------------------------------------
% Step 6: clean and save workspace, figures, and variable final_connection_points ----------------------------------------------- -----------------------------------------------
clear each; clear eachsad; clear i; clear ii; clear pp; clear r; clear rr; clear s; clear temp; clear sz; clear w; clear ww; clear Z; clear CONTOUR_DIRECTION_THRESHOLD;
cd(outFolder)
if strcmp (geom_name, 'GPeGPiSTNSN') ==1  
    save([subj '_' side '_GPeGPiSTNSN_Output_Workspace'])
    save([subj '_' side '_GPeGPiSTNSN_Output_axonpop'], 'axonpop')
end

if strcmp (geom_name, 'GPeSTNSN') ==1  
    save([subj '_' side '_GPeSTNSN_Output_Workspace'])
    save([subj '_' side '_GPeSTNSN_Output_axonpop'], 'axonpop')
end


if strcmp (geom_name, 'PutGPe') ==1  
    save([subj '_' side '_PutGPe_Output_Workspace'])
    save([subj '_' side '_PutGPe_Output_axonpop'], 'axonpop')
end

if strcmp (geom_name, 'PutGPeGPi') ==1  
    save([subj '_' side '_PutGPeGPi_Output_Workspace'])
    save([subj '_' side '_PutGPeGPi_Output_axonpop'], 'axonpop')
end

 % Added by JL
if strcmp (geom_name, 'GPeSTN') ==1  
    save([subj '_' side '_GPeSTN_Output_Workspace'])
    save([subj '_' side '_GPeSTN_Output_axonpop'], 'axonpop')
end
 
 
%% ------------------------------------------------------------------------------------------------------------------------------
 % Step 7: make cell body within axonpop matrix for use later--------------------------------------------------------------------
 % make cell body matrix for all cells
 %GPe

cellbod_final=final_connection_pts1(1:3,:,1)';
 
if strcmp (geom_name, 'GPeGPiSTNSN') ==1
    % add to axonpop struct
    GPex=num2cell(cellbod_final(:,1));
        [axonpop(:).cellbod_GPex]= GPex{:};
    GPey=num2cell(cellbod_final(:,2));
        [axonpop(:).cellbod_GPey]= GPey{:};
    GPez=num2cell(cellbod_final(:,3));
        [axonpop(:).cellbod_GPez]= GPez{:};
    type=num2cell(Ax_choice);
        [axonpop(:).type]= type{:};
end

if strcmp (geom_name, 'GPeSTNSN') ==1
    % add to axonpop struct
    GPex=num2cell(cellbod_final(:,1));
        [axonpop(:).cellbod_GPex]= GPex{:};
    GPey=num2cell(cellbod_final(:,2));
        [axonpop(:).cellbod_GPey]= GPey{:};
    GPez=num2cell(cellbod_final(:,3));
        [axonpop(:).cellbod_GPez]= GPez{:};
end

if strcmp (geom_name, 'PutGPe') ==1 || strcmp (geom_name, 'PutGPeGPi') ==1
    % add to axonpop struct
    Putx=num2cell(cellbod_final(:,1));
        [axonpop(:).cellbod_Putx]= Putx{:};
    Puty=num2cell(cellbod_final(:,2));
        [axonpop(:).cellbod_Puty]= Puty{:};
    Putz=num2cell(cellbod_final(:,3));
        [axonpop(:).cellbod_Putz]= Putz{:};
end

% Added by JL
if strcmp (geom_name, 'GPeSTN') ==1
    % add to axonpop struct
    GPex=num2cell(cellbod_final(:,1));
        [axonpop(:).cellbod_GPex]= GPex{:};
    GPey=num2cell(cellbod_final(:,2));
        [axonpop(:).cellbod_GPey]= GPey{:};
    GPez=num2cell(cellbod_final(:,3));
        [axonpop(:).cellbod_GPez]= GPez{:};
end

%% ------------------------------------------------------------------------------------------------------------------------------
% Step 8: delete intersecting axons and save final cleaned axonpop struct with cell body location------------------------------
% Temporarily disabled by JL   
% for i = 1:length(nAxon)
%     cd(['\Users\', user, '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data Processing\STL Files']) 
%     axons_that_intersect = []; 
%     surface_struct_DWI_space  = stlread2([subj '_' side '_encap_COMSOL.stl']);
%     axons_that_intersect = get_axons_that_intersect(axonpop,surface_struct_DWI_space,false); %if axon numbers don't make sense, change input from false to true
%     axonpop(axons_that_intersect)=[];
% end

% outFolder = ['C:\Users\' user '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data Processing\Axons\Cleaned\'];
outFolder = 'C:\Users\jaejin\gdrive_jaejin\Work\PAM\Barb\Modeling_Data_rev1\Data_Processing\Axons\Cleaned\';
if ~exist(outFolder,'dir')
mkdir(outFolder)
end
cd(outFolder)
save(strcat(subj, '_', side, '_', geom_name, '_n', num2str(nAxon)), 'axonpop')
savefig(h1, strcat(subj, '_', side, '_', geom_name, '_n', num2str(nAxon)))
%save(strcat(subj, '_', side, '_', geom_name{1}), 'axonpop')
