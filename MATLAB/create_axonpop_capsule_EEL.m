% Program that uses the input varibles from the user and the XYZ
% coordinates of conturs to create axon geometries.
% Julia Slopsema 1/17/2019, adapted EEL 3/22 for GPDBS
% Use for capusle pathways in GP DBS (IC M1, IC SMA, IC PM).

%% ------------------------------------------------------------------------------------------------------------------------------
 %Step 0: initialize data and load inputs----------------------------------------------------------------------------------------
clearvars; close all; clc;
user = 'elecy'; 
nAxon=10; % number of axons to make
nSeedPoints=floor(nAxon*8); % number of seed points
% Use more seed points than axons to prevent 'wiggle' axons 
subject_ID = 'UD1018_PD090';
side = 'R';
subj= 'UD1018';
geom_name={'IC'};

inFolder = ['C:\Users\' user '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data Processing\'];
outFolder = ['C:\Users\' user '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data Processing\Axons\'];
lead_foldername=(['C:\Users\' user '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data Processing\STL Files\']);  

 
ii = 1;
cd(['C:\Users\' user '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data Processing\Contours\'])  

in_foldername = [inFolder 'Contours\' side '_' geom_name{ii}];
tract = geom_name{ii};
if ~exist(outFolder,'dir')
    mkdir(outFolder)
end
%% ------------------------------------------------------------------------------------------------------------------------------
    % Step 1: import points, make contours, and plot them
    nContour = numel(dir([in_foldername, '/*.txt'])); % number of contours
    nContour_pts = 100; % number of points per contour
    
    % preallocate mem
    contour = zeros(nContour_pts+1,3,nContour); % nContour_pts+1 b/c one is added to close the contour
    final_pts = zeros(nSeedPoints,3,nContour);
    
    
    kk=1;
    for iContour=1:nContour % use less contours if results are too wiggly
        
        %edits made by annie 10/19/2020 to account for contours with less than or more than 100 pts
        cont = load([in_foldername '/' num2str(iContour),'.txt']); % load contour points
        if length(cont) >=100
            contour(1:nContour_pts,:,kk)=cont(1:100,:); % load contour points
        elseif length(cont) == 99
            contour(1:nContour_pts,:,kk)=[cont(1:99,:);cont(99,:)];
        end
        
        normal_of_contour_plane = affine_fit(contour(1:nContour_pts,:,kk));
        contour(nContour_pts+1,:,kk)=contour(1,:,kk);                              % close the contours
        plot3(contour(:,1,kk),contour(:,2,kk),contour(:,3,kk)); hold on            % plot the contours
              
        
        % Use axon_1_populate_contour_xy for countours that run along the z-axis
        % and axon_1_populate_contour_xz for contours that run along the x- or
        % y-axis. Either can be used on non-vertical or non-horizonatal running
        % contours.
        CONTOUR_DIRECTION_THRESHOLD = 0.7;
        if (abs(dot(normal_of_contour_plane,[0,0,1])) > CONTOUR_DIRECTION_THRESHOLD)
            final_pts(:,:,kk)=axon_1_populate_contour_xy(contour(:,:,kk),nSeedPoints); % populate contour with random points
        elseif (abs(dot(normal_of_contour_plane,[0,1,0])) > CONTOUR_DIRECTION_THRESHOLD)
            final_pts(:,:,kk)=axon_1_populate_contour_xz(contour(:,:,kk),nSeedPoints); % populate contour with random points
        else
            final_pts(:,:,kk)=axon_1_populate_contour_yz(contour(:,:,kk),nSeedPoints); % populate contour with random points
        end
        plot3(final_pts(:,1,kk),final_pts(:,2,kk),final_pts(:,3,kk),'r.');  % plot contours
        kk=kk+1;
%         fprintf('contour #%02d...\n',iContour)
%         pause(0.5)
    end; hold off
    
%% ------------------------------------------------------------------------------------------------------------------------------
    % Step 2: Connect the random points across contours
    dist_between_contours =  sum(squeeze(mean(contour(:,:,1:end-1),1)-mean(contour(:,:,2:end),1)).^2,1).^(1/3);
   
    fiber_pts = axon_4_connect_contours(final_pts,nSeedPoints,nAxon);
            fiber_start = ones(1,nAxon); 
            fiber_end= zeros(1,length(fiber_pts));
            
    for i =1:nAxon
        fiber_pts{1,i} = flipud(fiber_pts{1,i}(fiber_start(i):end-fiber_end(i),:));
    end
    
    % plot axons
    h2= figure;
    for iAxon=1:nAxon
        plot3(fiber_pts{iAxon}(:,1),fiber_pts{iAxon}(:,2),fiber_pts{iAxon}(:,3))
        hold on;
    end
    cd(lead_foldername)
    nuc = stlread2(['C:\Users\', user, '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data Processing\STL Files\' subj '_' side '_GPi_DTIspace.stl']);
    patch('Faces',nuc.faces,'Vertices',nuc.vertices,'FaceColor',[1 0 0],'EdgeColor',[1 0 0],'FaceAlpha',0.2,'EdgeAlpha',0.2);
    hold on;
    for iContour=1:nContour
    plot3(final_pts(:,1,iContour),final_pts(:,2,iContour),final_pts(:,3,iContour),'r.');  % plot contours
    end 
    hold off; axis equal; axis tight; title('3D contour');
    xlabel('x (mm)');ylabel('y (mm)');zlabel('z (mm)'); view(3)
    savefig(h2,[outFolder '/' subj '_' side '_' tract '.fig'])

   
  %% ------------------------------------------------------------------------------------------------------------------------------
   % Step 3: Create axonpop structure for the axon population and save--------------------------------------------------------------
    axonpop = axon_4_axonpop_struct_fiberD(fiber_pts, 5.7);
    cd(outFolder)
    save([subj '_' side '_' tract '_axonpop.mat'],'axonpop')  % save workspace


 %% ------------------------------------------------------------------------------------------------------------------------------
    % Step 4: remove intersecting axons and save files----------------------------------------------------------------------------
for i = 1:length(geom_name)
    cd(['\Users\', user, '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data Processing\STL Files']) 
    axons_that_intersect = []; 
    surface_struct_DWI_space  = stlread2([subj '_' side '_encap_COMSOL.stl']);
    axons_that_intersect = get_axons_that_intersect(axonpop,surface_struct_DWI_space,false); %if axon numbers don't make sense, change input from false to true

    axonpop(axons_that_intersect)=[];
    outFolder = ['C:\Users\' user '\Box\DBS_Patient_Modeling\Globus_Pallidus\' subject_ID '\Data Processing\Axons\Cleaned\'];
    mkdir(outFolder)
    cd(outFolder)
    save(strcat(subj, '_', side, '_', geom_name{i}), 'axonpop')
end

