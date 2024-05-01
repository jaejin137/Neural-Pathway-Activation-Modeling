function [GPe,GPiRing,STN] = PopGPeSTN(stl,nSeedPoints)
    bb.GPe.x = [min(stl.GPe.Points(:,1)),max(stl.GPe.Points(:,1))];
    bb.GPe.y = [min(stl.GPe.Points(:,2)),max(stl.GPe.Points(:,2))];
    bb.GPe.z = [min(stl.GPe.Points(:,3)),max(stl.GPe.Points(:,3))];
    bb.GPiRing.x = [min(stl.GPiRing.Points(:,1)),max(stl.GPiRing.Points(:,1))];
    bb.GPiRing.y = [min(stl.GPiRing.Points(:,2)),max(stl.GPiRing.Points(:,2))];
    bb.GPiRing.z = [min(stl.GPiRing.Points(:,3)),max(stl.GPiRing.Points(:,3))];
    bb.STN.x = [min(stl.STN.Points(:,1)),max(stl.STN.Points(:,1))];
    bb.STN.y = [min(stl.STN.Points(:,2)),max(stl.STN.Points(:,2))];
    bb.STN.z = [min(stl.STN.Points(:,3)),max(stl.STN.Points(:,3))];
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
    
    %% pop GPiRing
    h = waitbar(0,'Please wait...');
    clear GPiRing
    num.GPiRing = nSeedPoints;
    for i = 1:num.GPiRing
        
        % Update the waitbar
        waitbar(i/num.GPiRing,h,['Checking GPiRing cell: ',num2str(i)]);
        
        temp(1) = rand(1)*(bb.GPiRing.x(2)-bb.GPiRing.x(1)) + bb.GPiRing.x(1);
        temp(2) = rand(1)*(bb.GPiRing.y(2)-bb.GPiRing.y(1)) + bb.GPiRing.y(1);
        temp(3) = rand(1)*(bb.GPiRing.z(2)-bb.GPiRing.z(1)) + bb.GPiRing.z(1);
    
        % Check if cell is within the GPiRing
        while jordancurve(stl.GPiRingTri,temp)==0
    
            temp(1) = rand(1)*(bb.GPiRing.x(2)-bb.GPiRing.x(1)) + bb.GPiRing.x(1);
            temp(2) = rand(1)*(bb.GPiRing.y(2)-bb.GPiRing.y(1)) + bb.GPiRing.y(1);
            temp(3) = rand(1)*(bb.GPiRing.z(2)-bb.GPiRing.z(1)) + bb.GPiRing.z(1);
    
        end
        %GPiRing(i).xyz = temp;
        GPiRing(1:3,i)=temp(1:3);
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