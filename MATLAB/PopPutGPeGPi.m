function [Put,GPe,GPi,STN] = PopPutGPeGPi(stl,nSeedPoints)
    bb.GPe.x = [min(stl.GPe.Points(:,1)),max(stl.GPe.Points(:,1))];
    bb.GPe.y = [min(stl.GPe.Points(:,2)),max(stl.GPe.Points(:,2))];
    bb.GPe.z = [min(stl.GPe.Points(:,3)),max(stl.GPe.Points(:,3))];
    bb.Put.x = [min(stl.Put.Points(:,1)),max(stl.Put.Points(:,1))];
    bb.Put.y = [min(stl.Put.Points(:,2)),max(stl.Put.Points(:,2))];
    bb.Put.z = [min(stl.Put.Points(:,3)),max(stl.Put.Points(:,3))];
    bb.GPi.x = [min(stl.GPi.Points(:,1)),max(stl.GPi.Points(:,1))];
    bb.GPi.y = [min(stl.GPi.Points(:,2)),max(stl.GPi.Points(:,2))];
    bb.GPi.z = [min(stl.GPi.Points(:,3)),max(stl.GPi.Points(:,3))];
    bb.STN.x = [min(stl.STN.Points(:,1)),max(stl.STN.Points(:,1))];
    bb.STN.y = [min(stl.STN.Points(:,2)),max(stl.STN.Points(:,2))];
    bb.STN.z = [min(stl.STN.Points(:,3)),max(stl.STN.Points(:,3))];
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

%% pop GPi
  h = waitbar(0,'Please wait...');
    clear GPi
    num.GPi = nSeedPoints;
    for i = 1:num.GPi
    
    % Update the waitbar
    waitbar(i/num.GPi,h,['Checking GPi cell: ',num2str(i)]);
    
    temp(1) = rand(1)*(bb.GPi.x(2)-bb.GPi.x(1)) + bb.GPi.x(1);
    temp(2) = rand(1)*(bb.GPi.y(2)-bb.GPi.y(1)) + bb.GPi.y(1);
    temp(3) = rand(1)*(bb.GPi.z(2)-bb.GPi.z(1)) + bb.GPi.z(1);

    % Check if cell is within the STN
    while jordancurve(stl.GPiTri,temp)==0

        temp(1) = rand(1)*(bb.GPi.x(2)-bb.GPi.x(1)) + bb.GPi.x(1);
        temp(2) = rand(1)*(bb.GPi.y(2)-bb.GPi.y(1)) + bb.GPi.y(1);
        temp(3) = rand(1)*(bb.GPi.z(2)-bb.GPi.z(1)) + bb.GPi.z(1);

    end
    %GPi(i).xyz = temp;
    GPi(1:3,i)=temp(1:3);
    end
    close(h)   
%% find STN and SN seedpoints
    h = waitbar(0,'Please wait...');
    clear STN
    num.STN = nSeedPoints;

for i = 1:num.STN
    % Update the waitbar
    waitbar(i/num.STN,h,['Checking STN cell: ',num2str(i)]);
    
    temp(1) = rand(1)*(bb.STN.x(2)-bb.STN.x(1)) + bb.STN.x(1);
    temp(2) = rand(1)*(bb.STN.y(2)-bb.STN.y(1)) + bb.STN.y(1);
    temp(3) = rand(1)*(bb.STN.z(2)-bb.STN.z(1)) + bb.STN.z(1);

    % Check if cell is within the STN
    while jordancurve(stl.STNTri,temp)==0

        temp(1) = rand(1)*(bb.STN.x(2)-bb.STN.x(1)) + bb.STN.x(1);
        temp(2) = rand(1)*(bb.STN.y(2)-bb.STN.y(1)) + bb.STN.y(1);
        temp(3) = rand(1)*(bb.STN.z(2)-bb.STN.z(1)) + bb.STN.z(1);

    end
    
    %STN(i).xyz = temp;
    STN(1:3,i)=temp(1:3);

end
    close(h)

% h = waitbar(0,'Please wait...');
%     clear SN
%     num.SN = nSeedPoints;for i = 1:num.SN
%     
%     % Update the waitbar
%     waitbar(i/num.SN,h,['Checking SN cell: ',num2str(i)]);
%     
%     temp(1) = rand(1)*(bb.SN.x(2)-bb.SN.x(1)) + bb.SN.x(1);
%     temp(2) = rand(1)*(bb.SN.y(2)-bb.SN.y(1)) + bb.SN.y(1);
%     temp(3) = rand(1)*(bb.SN.z(2)-bb.SN.z(1)) + bb.SN.z(1);
% 
%     % Check if cell is within the STN
%     while jordancurve(stl.SNTri,temp)==0
% 
%         temp(1) = rand(1)*(bb.SN.x(2)-bb.SN.x(1)) + bb.SN.x(1);
%         temp(2) = rand(1)*(bb.SN.y(2)-bb.SN.y(1)) + bb.SN.y(1);
%         temp(3) = rand(1)*(bb.SN.z(2)-bb.SN.z(1)) + bb.SN.z(1);
% 
%     end
%     %SN(i).xyz = temp;
%     SN(1:3,i)=temp(1:3);
%     end
%     close(h)
%     
end % end function