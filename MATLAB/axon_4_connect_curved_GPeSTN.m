function [fiber_pts,final_connection_pts1] = axon_4_connect_curved_GPeSTN (final_pts,nAxon,nSeedPoints)

%final pts, seed pts, axons
% Emily Lecy (Nov 2022)
%
% This program connects the Nr random points from a final_pts variable
%  and outputs a 3 x Nr x length of spline for fiber_pts_spline.
%
% Variables:
%   final_pts= 3d coordinated of all points in each contour
%   nContour= # of contours
%   nAxon= number of desired axons
%   fiber_pts= output of points to connect!
%   ALL ERRORS ARE JUST PREALLOCATION-OK.

Total_connect1=4; 
final_connection_pts1 = zeros(3,nAxon,(Total_connect1));
final_connection_pts = zeros(3,nAxon);

nuc1ptsleft = 1:nSeedPoints; %countdown for contour GPe
nuc2ptsleft = 1:nSeedPoints; %countdown for contour GPi
nuc3ptsleft = 1:nSeedPoints; %countdown for contour STN

%% Step 0: Choose random points for your axons to begin in GPe--------------------------------------------------------------------- 
%make sure these random points are more lateral in GPe, otherwise you get
%issues with connections
first=randperm(nSeedPoints); first=first(1:nAxon)';
  for w=1:length(final_connection_pts(:,:,1))
      ww=first(w,1);
      final_connection_pts(1:3,w,1)=final_pts(1:3,ww,1);
      nuc1ptsleft = nuc1ptsleft(nuc1ptsleft~=ww); %remove points used from nuc so they cannot be reused
  end 
    %find bottom 20% of y axis points, eliminate 
%   minY=min(final_pts(2,:,1));
%   Yfinal=minY*.10; Yfinal=minY+Yfinal;
% for ww=1:dubs
%     pt1=final_connection_pts(2,ww);
%     if pt1<Yfinal
%         bad(ww,:)=1;
%     else
%        bad(ww,:)=0;
%     end
% end
% kick=bad(:,1)==1;
% final_connection_pts(:,kick)= [];

%% Step 1: Distribute each axon based on connection points in smaller --> bigger spheres ----------------------------------------------
h = waitbar(0,'Please wait...');

for w=1:nAxon
        waitbar(w/nAxon,h,['Connecting Contour: ',num2str(w)]);
        pt1=final_connection_pts(1:3,w,1); [~ , ~]=size(nuc1ptsleft);
        for r=1:length(nuc1ptsleft) 
            rr=nuc1ptsleft(1,r);
            pt2=final_pts(1:3,rr,1);
            dist(r,:) = sqrt((pt1(1)-pt2(1))^2+(pt1(2)-pt2(2))^2+(pt1(3)-pt2(3))^2); %can be zero if you take your own point....ugh.....
        end 
        nucdist=(nuc1ptsleft)';
        distfin=cat(2,dist,nucdist);
        distfin=sortrows(distfin,1);
        select3=round(length(dist)*0.08); select3=round(select3); %find random point in 10% closest axons
        sphere3=distfin(1:select3,:); % sort distances        
        final_connection_pts1(1:3,w,1)=final_connection_pts(1:3,w);

        %another GPe point to connect to randomly that is within 10%
        %closest to it
            sphere3num=length(sphere3);
            sphere3_final=randi(numel(1:sphere3num)); 
            pick=sphere3(sphere3_final,2); % final_pts row of sphere 2 connections
            nuc1ptsleft = nuc1ptsleft(nuc1ptsleft~=pick); %remove first points used from nuc so they cannot be reused
            final_connection_pts1(1:3,w,2)=final_pts(1:3,pick,1);
            clear dist; clear nucdist; clear distfin; clear final_sphere; clear pick_sphere1; clear pick_sphere2; clear pick_sphere3;

        % connect GPe to nearest GPiRing
        pt1=final_connection_pts1(1:3,w,2);
        for r=1:length(nuc2ptsleft) 
            rr=nuc2ptsleft(1,r);
            pt2=final_pts(1:3,rr,2);
            dist(r,:) = sqrt((pt1(1)-pt2(1))^2+(pt1(2)-pt2(2))^2+(pt1(3)-pt2(3))^2); 
        end 
        nucdist=(nuc2ptsleft)';
        dist=cat(2,dist,nucdist);
        distfin=sortrows(dist,1); 
        rando=randi(numel(1:8));
        pick=distfin(rando,2);
        nuc2ptsleft = nuc2ptsleft(nuc2ptsleft~=pick); %remove first points used from nuc so they cannot be reused
        final_connection_pts1(1:3,w,3)=final_pts(1:3,pick,2);
        clear dist; clear nucdist; clear distfin; clear final_sphere; clear pick_sphere1; clear pick_sphere2; clear pick_sphere3;
   
        %pick closest STN
        pt1=final_connection_pts1(1:3,w,3);
        for r=1:length(nuc3ptsleft) 
            rr=nuc3ptsleft(1,r);
            pt2=final_pts(1:3,rr,3);
            dist(r,:) = sqrt((pt1(1)-pt2(1))^2+(pt1(2)-pt2(2))^2+(pt1(3)-pt2(3))^2); %can be zero if you take your own point....ugh.....
        end 
        nucdist=(nuc3ptsleft)';
        dist=cat(2,dist,nucdist);
        distfin=sortrows(dist,1); 
        pick=distfin(1,2);
        nuc3ptsleft = nuc3ptsleft(nuc3ptsleft~=pick); %remove first points used from nuc so they cannot be reused
        final_connection_pts1(1:3,w,4)=final_pts(1:3,pick,3);
        clear dist; clear nucdist; clear distfin; clear final_sphere; clear pick_sphere1; clear pick_sphere2; clear pick_sphere3;
end
close(h)
%% STEP 2: clean up grids so they  only have axon connection points in them (used to have zeroes)--------------------------------------------------------------------------------------------
FC1_1=final_connection_pts1(1:3,:,1); FC1_2=final_connection_pts1(1:3,:,2); FC1_3=final_connection_pts1(1:3,:,3); FC1_4=final_connection_pts1(1:3,:,4);
FC1_1=FC1_1(:,any(FC1_1)); FC1_2=FC1_2(:,any(FC1_2)); FC1_3=FC1_3(:,any(FC1_3)); FC1_4=FC1_4(:,any(FC1_4)); 
clear final_connection_pts1;
final_connection_pts1(:,:,1)=FC1_1; final_connection_pts1(:,:,2)=FC1_2; final_connection_pts1(:,:,3)=FC1_3; final_connection_pts1(:,:,4)=FC1_4; 
%% Step 4: interpolate-------------------------------------------------------------------------------------------------------------------
% initalize waitbar

[~, um, ~]=size(final_connection_pts1);

for iAxon=1:um
    rez = 100; % (pts/mm) approximatly 100 points per mm
    %nSplinePoints = dist_between_contours*rez;
    xyz = squeeze(final_connection_pts1(:,iAxon,:))';
    CS = cat(1,0,cumsum(sqrt(sum(diff(xyz,[],1).^2,2))));
    fiber_pts{iAxon} = interp1(CS, xyz, unique([CS(:)' linspace(0,CS(end),CS(end)*rez)]),'makima');
end
end %function end
