%function [STN,STNMCPath,MC] = PopSTNMC(stl,nSeedPoints)
function [STN,MC] = PopSTNMC(stl,nSeedPoints)
    bb.STN.x = [min(stl.STN.Points(:,1)),max(stl.STN.Points(:,1))];
    bb.STN.y = [min(stl.STN.Points(:,2)),max(stl.STN.Points(:,2))];
    bb.STN.z = [min(stl.STN.Points(:,3)),max(stl.STN.Points(:,3))];
%    bb.STNMCPath.x = [min(stl.STNMCPath.Points(:,1)),max(stl.STNMCPath.Points(:,1))];
%    bb.STNMCPath.y = [min(stl.STNMCPath.Points(:,2)),max(stl.STNMCPath.Points(:,2))];
%    bb.STNMCPath.z = [min(stl.STNMCPath.Points(:,3)),max(stl.STNMCPath.Points(:,3))];
    bb.MC.x = [min(stl.MC.Points(:,1)),max(stl.MC.Points(:,1))];
    bb.MC.y = [min(stl.MC.Points(:,2)),max(stl.MC.Points(:,2))];
    bb.MC.z = [min(stl.MC.Points(:,3)),max(stl.MC.Points(:,3))];

    %% pop STN
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

%    %% pop STNMCPath
%    h = waitbar(0,'Please wait...');
%    clear STNMCPath
%    num.STNMCPath = nSeedPoints;
%    for i = 1:num.STNMCPath
%        
%        % Update the waitbar
%        waitbar(i/num.STNMCPath,h,['Checking STNMCPath cell: ',num2str(i)]);
%        
%        temp(1) = rand(1)*(bb.STNMCPath.x(2)-bb.STNMCPath.x(1)) + bb.STNMCPath.x(1);
%        temp(2) = rand(1)*(bb.STNMCPath.y(2)-bb.STNMCPath.y(1)) + bb.STNMCPath.y(1);
%        temp(3) = rand(1)*(bb.STNMCPath.z(2)-bb.STNMCPath.z(1)) + bb.STNMCPath.z(1);
%    
%        % Check if cell is within the STNMCPath
%        while jordancurve(stl.STNMCPathTri,temp)==0
%    
%            temp(1) = rand(1)*(bb.STNMCPath.x(2)-bb.STNMCPath.x(1)) + bb.STNMCPath.x(1);
%            temp(2) = rand(1)*(bb.STNMCPath.y(2)-bb.STNMCPath.y(1)) + bb.STNMCPath.y(1);
%            temp(3) = rand(1)*(bb.STNMCPath.z(2)-bb.STNMCPath.z(1)) + bb.STNMCPath.z(1);
%    
%        end
%        %STNMCPath(i).xyz = temp;
%        STNMCPath(1:3,i)=temp(1:3);
%    end
%    close(h)

    %% pop MC
    h = waitbar(0,'Please wait...');
    clear MC
    num.MC = nSeedPoints;
    for i = 1:num.MC
    
        % Update the waitbar
        waitbar(i/num.MC,h,['Checking MC cell: ',num2str(i)]);
        
        temp(1) = rand(1)*(bb.MC.x(2)-bb.MC.x(1)) + bb.MC.x(1);
        temp(2) = rand(1)*(bb.MC.y(2)-bb.MC.y(1)) + bb.MC.y(1);
        temp(3) = rand(1)*(bb.MC.z(2)-bb.MC.z(1)) + bb.MC.z(1);
    
        % Check if cell is within the STN
        while jordancurve(stl.MCTri,temp)==0
    
            temp(1) = rand(1)*(bb.MC.x(2)-bb.MC.x(1)) + bb.MC.x(1);
            temp(2) = rand(1)*(bb.MC.y(2)-bb.MC.y(1)) + bb.MC.y(1);
            temp(3) = rand(1)*(bb.MC.z(2)-bb.MC.z(1)) + bb.MC.z(1);

        end
        %MC(i).xyz = temp;
        MC(1:3,i)=temp(1:3);
    end
    close(h)
    
%   h = waitbar(0,'Please wait...');
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

end % end function
