function fiber_pts = axon_4_connect_contours (contourpts,nSeedPoints,nAxon)

% axon_2_connect_contours()
%
% Matt Johnson (June 2007)
% edited by Allison Connolly (Oct 2011)
%
% This program connects the Nr random points generated by
% axon_1_populate_contour from one contour to the next and outputs a 3 x Nr
% x length of spline for fiber_pts_spline.
%
% Variables:
%   Nr=number of random points per contour
%   contourpts=[# of pts,[x,y,z] location of pts, # of contours]
%   p1=axon length multiplyer, you will have to manually change this
%   variable until your axons seem like the correct length

nContour=size(contourpts,3); % # of contours 
for iContour = 1:nContour-1 %for each contour

    seed_xvalues = squeeze(contourpts(:,1,iContour))';
    seed_yvalues = squeeze(contourpts(:,2,iContour))';
    seed_zvalues = squeeze(contourpts(:,3,iContour))';
    next_contour_xvalues = squeeze(contourpts(:,1,iContour+1))';
    next_contour_yvalues = squeeze(contourpts(:,2,iContour+1))';
    next_contour_zvalues = squeeze(contourpts(:,3,iContour+1))';

    % Create coordinate indices
    contour_xindices = next_contour_xvalues;
    contour_yindices = next_contour_yvalues;
    contour_zindices = next_contour_zvalues;

    % Define variables
    seed_pts=zeros(3,nSeedPoints);

    % Calculate distance between each seed point in one contour and 
    % all other points in the next contour
    dist = zeros(length(contour_xindices),nSeedPoints);
    for iSeedPoints = 1:nSeedPoints
        if iContour == 1 % in first contour
            seed_pts(:,iSeedPoints) = [seed_xvalues(iSeedPoints),seed_yvalues(iSeedPoints),seed_zvalues(iSeedPoints)]';
        else
            seed_pts(:,iSeedPoints) = min_connection_pts(:,iSeedPoints);
        end
        for ii = 1:length(contour_xindices)
            xd = contour_xindices(ii) - seed_pts(1,iSeedPoints);
            yd = contour_yindices(ii) - seed_pts(2,iSeedPoints);
            zd = contour_zindices(ii) - seed_pts(3,iSeedPoints);
            dist(ii,iSeedPoints) = sqrt(xd^2 + yd^2 + zd^2);
        end
    end
    
    c1ptsleft = 1:nSeedPoints; %countdown for contour 1
    c2ptsleft = 1:nSeedPoints; %countdown for contour 2
    for pt = 1:nSeedPoints
        
        clear mindist;
        ctr = 1;
        for iSeedPoints = c1ptsleft,
            mindist(ctr) = min(dist(c2ptsleft,iSeedPoints)); %for each c1 pt, the minimum distance to any c2 pt
            ctr = ctr + 1;
        end
        
        [~,maxind] = max(mindist); % the index of the c1 pt farthest from any c2 pt
        %c1_index = c1ptsleft(maxind); % choose the index out of c1 pts left that we will find a c2 match pt for
        c1_index = c1ptsleft(randperm(length(c1ptsleft),1));%(maxind); % choose the index out of c1 pts left that we will find a c2 match pt for

        [~,minind] = min(dist(c2ptsleft,c1_index)); % the index of the c2 pt closest to the c1 pt chosen above
        c2_index = c2ptsleft(minind); % index out of c2 pts left
        min_connection_pts(:,c1_index) = [contour_xindices(c2_index),...
            contour_yindices(c2_index),contour_zindices(c2_index)]'; % save the connection
        
        c1ptsleft = c1ptsleft(c1ptsleft~=c1_index); %remove the used pt from contour 1
        c2ptsleft = c2ptsleft(c2ptsleft~=c2_index); %remove the used pt from contour 2
    end
         
    if iContour == 1
        final_connection_pts(:,:,1) = [seed_xvalues ; seed_yvalues ; seed_zvalues];
    end
    X1 = squeeze(min_connection_pts(1,:));
    Y1 = squeeze(min_connection_pts(2,:));
    Z1 = squeeze(min_connection_pts(3,:));
    final_connection_pts(:,:,iContour+1) = [X1;Y1;Z1];
end 

% initalize waitbar
hWait=waitbar(0,'1','Name','Connecting points...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(hWait,'canceling',0)

for iAxon=1:nAxon
    
    % Check for waitbar Cancel button press
    if getappdata(hWait,'canceling');break;end
    
    % Report progress in the waitbar's message field
    waitbar(iAxon/nAxon,hWait,['Connecting points for axon ' num2str(iAxon) ' of ' num2str(nAxon)])
    
    rez = 100; % (pts/mm) approximatly 100 points per mm
    %nSplinePoints = dist_between_contours*rez;
    xyz = squeeze(final_connection_pts(:,iAxon,:))'; 
    CS = cat(1,0,cumsum(sqrt(sum(diff(xyz,[],1).^2,2))));
    fiber_pts{iAxon} = interp1(CS, xyz, unique([CS(:)' linspace(0,CS(end),CS(end)*rez)]),'pchip');
end

delete(hWait) % DELETE the waitbar; don't try to CLOSE it.

end %function end


