function [efPxTrace,fPxTrace]=mnl_MakeSeperateMasksForEachImage(Trace,Scale,dim)
% Function to make maps for each trace 
% Inputs
% Trace - The Sturctured Matrix with the points
% Scale - The scale
% Dim - dimensions of the Image

%% Step 1 - Convert the Numbers back to Pixels
szT=size(Trace,2);
PxTrace=struct('Points',[],'Diameter',[],'OriginalTrace',[],'LengthVx',[],'LengthReal',[]);

for i=1:szT
    PxTrace(i).Points=Trace(i).RealPoints./Scale;
    PxTrace(i).Diameter=Trace(i).Diameter;
    PxTrace(i).OriginalTrace=Trace(i).OriginalTrace;
    PxTrace(i).LengthVx=Trace(i).LengthVx;
    PxTrace(i).LengthReal=Trace(i).LengthReal;
    PxTrace(i).SetId=Trace(i).SetId;
    PxTrace(i).TypeId=Trace(i).TypeId;
end
%Flip Y and Z because they are negative
fPxTrace=struct('Points',[]);
for i=1:szT
    XPointList=PxTrace(i).Points(:,1);
    YPointList=PxTrace(i).Points(:,2)*-1;
    ZPointList=PxTrace(i).Points(:,3)*-1;
    fPxTrace(i).Points=[round(XPointList) round(YPointList) round(ZPointList)];
    clear XPointList YPointList ZPointList
    fPxTrace(i).Diameter=PxTrace(i).Diameter;
    fPxTrace(i).OriginalTrace=PxTrace(i).OriginalTrace;
    fPxTrace(i).LengthVx=PxTrace(i).LengthVx;
    fPxTrace(i).LengthReal=PxTrace(i).LengthReal;
    fPxTrace(i).SetId=PxTrace(i).SetId;
    fPxTrace(i).TypeId=PxTrace(i).TypeId;
end
%% Step 2 - Make Mask
xlimit=dim(1);
ylimit=dim(2);
zlimit=dim(4);
efPxTrace=struct('Points',[],'Diameter',[],'OriginalTrace',[],'LengthVx',[],'LengthReal',[],'AllVoxels',[]);
for i=1:szT
    efPxTrace(i).Points=fPxTrace(i).Points;
    efPxTrace(i).Diameter=fPxTrace(i).Diameter;
    efPxTrace(i).OriginalTrace=fPxTrace(i).OriginalTrace;
    efPxTrace(i).LengthVx=fPxTrace(i).LengthVx;
    efPxTrace(i).LengthReal=fPxTrace(i).LengthReal;
    efPxTrace(i).SetId=fPxTrace(i).SetId;
    efPxTrace(i).TypeId=fPxTrace(i).TypeId;
    nPoints=size(fPxTrace(i).Points,1);
    %Alternate Grid Strat
    n=1;
    for j=1:nPoints
        rXY=Trace(i).Diameter(j); %Extract the diameter of the axon
        r1=round(rXY/Scale(1)/2); %Convert to pixels
        zr=r1+(round(2*Scale(3)/2)); %Z variable in the future shold be variable on the NA of the lens
        r=[r1 r1 zr]; % radius spread in [x y z] pixels        
        %Central Positions
        Xpos=fPxTrace(i).Points(j,1)+1; %Add one to counter the fact that MATLAB starts from (1,1,1) not (0,0,0)
        Ypos=fPxTrace(i).Points(j,2)+1; %Add one to counter the fact that MATLAB starts from (1,1,1) not (0,0,0)
        Zpos=fPxTrace(i).Points(j,3)+1; %Add one to counter the fact that MATLAB starts from (1,1,1) not (0,0,0)
        %limits
        Xmin=Xpos-r(1);Xmax=Xpos+r(1);
        if Xmin<1
            Xmin=1;
        elseif Xmin>xlimit
            Xmin=xlimit;
        end
        if Xmax<1
            Xmax=1;
        elseif Xmax>xlimit
            Xmax=xlimit;
        end
        
        Ymin=Ypos-r(2);Ymax=Ypos+r(2);
        if Ymin<1
            Ymin=1;
        elseif Ymin>ylimit
            Ymin=ylimit;
        end
        if Ymax<1
            Ymax=1;
        elseif Ymax>ylimit
            Ymax=ylimit;
        end
        
        Zmin=Zpos-r(3);Zmax=Zpos+r(3);
        if Zmin<1
            Zmin=1;
        elseif Zmin>zlimit
            Zmin=zlimit;
        end
        if Zmax<1
            Zmax=1;
        elseif Zmax>zlimit
            Zmax=zlimit;
        end
        %Grid Strat
        [x,y,z]=ndgrid(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax);
        tVx=[x(:),y(:),z(:)];
        numVx=size(tVx,1);
        n2=n+numVx-1;
        tAllVoxels(n:n2,:)=tVx;
        clear tVx
        n=n2+1;
    end
    %Now remove the duplicates
    [AllVoxels]=unique(tAllVoxels,'rows');
    efPxTrace(i).AllVoxels=AllVoxels;
    clear AllVoxels tAllVoxels
end
end
