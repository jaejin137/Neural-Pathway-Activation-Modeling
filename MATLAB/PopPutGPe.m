function [Put,GPe] = PopPutGPe(stl,nSeedPoints)
    bb.GPe.x = [min(stl.GPe.Points(:,1)),max(stl.GPe.Points(:,1))];
    bb.GPe.y = [min(stl.GPe.Points(:,2)),max(stl.GPe.Points(:,2))];
    bb.GPe.z = [min(stl.GPe.Points(:,3)),max(stl.GPe.Points(:,3))];
    bb.Put.x = [min(stl.Put.Points(:,1)),max(stl.Put.Points(:,1))];
    bb.Put.y = [min(stl.Put.Points(:,2)),max(stl.Put.Points(:,2))];
    bb.Put.z = [min(stl.Put.Points(:,3)),max(stl.Put.Points(:,3))];
%% pop Put
  h = waitbar(0,'Please wait...');
    clear Put
    num.Put = nSeedPoints;
    for i = 1:num.Put
    
    % Update the waitbar
    waitbar(i/num.Put,h,['Checking Put cell: ',num2str(i)]);
    
    temp(1) = rand(1)*(bb.Put.x(2)-bb.Put.x(1)) + bb.Put.x(1);
    temp(2) = rand(1)*(bb.Put.y(2)-bb.Put.y(1)) + bb.Put.y(1);
    temp(3) = rand(1)*(bb.Put.z(2)-bb.Put.z(1)) + bb.Put.z(1);

    % Check if cell is within the STN
    while jordancurve(stl.PutTri,temp)==0

        temp(1) = rand(1)*(bb.Put.x(2)-bb.Put.x(1)) + bb.Put.x(1);
        temp(2) = rand(1)*(bb.Put.y(2)-bb.Put.y(1)) + bb.Put.y(1);
        temp(3) = rand(1)*(bb.Put.z(2)-bb.Put.z(1)) + bb.Put.z(1);

    end
    %Put(i).xyz = temp;
    Put(1:3,i)=temp(1:3);
    end
    close(h)   
%% pop GPe
  h = waitbar(0,'Please wait...');
    clear GPe
    num.GPe = nSeedPoints;
    for i = 1:num.GPe
    
    % Update the waitbar
    waitbar(i/num.GPe,h,['Checking GPe cell: ',num2str(i)]);
    
    temp(1) = rand(1)*(bb.GPe.x(2)-bb.GPe.x(1)) + bb.GPe.x(1);
    temp(2) = rand(1)*(bb.GPe.y(2)-bb.GPe.y(1)) + bb.GPe.y(1);
    temp(3) = rand(1)*(bb.GPe.z(2)-bb.GPe.z(1)) + bb.GPe.z(1);

    % Check if cell is within the STN
    while jordancurve(stl.GPeTri,temp)==0

        temp(1) = rand(1)*(bb.GPe.x(2)-bb.GPe.x(1)) + bb.GPe.x(1);
        temp(2) = rand(1)*(bb.GPe.y(2)-bb.GPe.y(1)) + bb.GPe.y(1);
        temp(3) = rand(1)*(bb.GPe.z(2)-bb.GPe.z(1)) + bb.GPe.z(1);

    end
    %GPe(i).xyz = temp;
    GPe(1:3,i)=temp(1:3);
    end
    close(h)
    
end % end function