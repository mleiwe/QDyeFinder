function [fTrace]=mnl_DoubleCheckLongTraceConsistency_v3(Trace,lim,LengthThresh,Data,BkgMean,MaxValues,Scale,EuThresh)
%Function to assess the colour consistency of a trace and if it is too
%different find the point where to cut it
%Inputs
% Trace - the structure of the traces extracted from Neurolucida
% lim - the minimum safe distance calculated from mnl_EvaluateFragmentsColourConsistency4
% LengthThresh - The length at which the code checks to see if the trace is inaccurate
% Data - the Image
% BkgMean - the mean background values
% MaxValues - The maximum values of the image
% Scale - The scale (x,y, and z)
% EuThresh - The Euclidean Threshold
%
%Outputs
%fTrace -Updated Trace structure
szT=size(Trace,2);
fN=1;
dim=size(Data);
nXFP=dim(3);
for i=1:szT
    %% Find Traces above the threshold
    if Trace(i).LengthReal>=LengthThresh
        %% Now evaluate the colour changes incrementally at each point
        NumPoints=size(Trace(i).Points,1);
        PointList=Trace(i).Points;
        [~,DistanceBetweenPoints]=mnl_CalculateLength(PointList);
        Distances(:,1)=1:NumPoints-1;
        Distances(:,2)=2:NumPoints;
        Distances(:,3)=DistanceBetweenPoints;
        for j=1:NumPoints-1            
            Point1=j;
            for k=j+1:NumPoints
                Point2=k;
                ListOfDistances=Distances(Point1:Point2-1,3);
                FinalDist=sum(ListOfDistances);
                if FinalDist>=lim
                    EdPoint=Point2;
                    break
                end
                if k==NumPoints
                    EdPoint=Point2;
                end
            end
            Group1_Points=PointList(Point1:EdPoint,:);
            DiameterList=Trace(i).Diameter(Point1:EdPoint);
            [ChanMean,~,~,~,VoxelVals,~]=mnl_ExtractKeyVoxels(Group1_Points,DiameterList,Scale,Data); %NB Points are in um
            Group_Mean(j,:)=ChanMean;
            tBkgSub=ChanMean-BkgMean;
            idx=tBkgSub<=0;
            tBkgSub(idx)=0;
            Group_BkgSub(j,:)=tBkgSub;
            if sum(tBkgSub)>0
                Group_NormMean(j,:)=Group_BkgSub(j,:)./MaxValues;
                Group_VecNorm(j,:)=mnl_NormaliseVectors(Group_NormMean(j,:));
            else
                Group_NormMean(j,:)=zeros(1,nXFP);
                Group_VecNorm(j,:)=zeros(1,nXFP);
            end
            Group(j).Points=Group1_Points;
            Group(j).Diameters=DiameterList;
        end
        szGroups=size(Group_VecNorm,1);
        j=1;
        Gp1=Group_VecNorm(j,:);
        Gp1PointsDiam(:,1:3)=[Group(j).Points];
        Gp1PointsDiam(:,4)=Group(j).Diameters;
        while j<szGroups            
            Gp2=Group_VecNorm(j+1,:);
            Gp2PointsDiam(:,1:3)=[Group(j+1).Points];
            Gp2PointsDiam(:,4)=Group(j+1).Diameters;
            
            EuD=mnl_EuclideanDistance(Gp1,Gp2);
            if EuD<=EuThresh
                %Merge the Groups              
                Merge=[Gp1PointsDiam;Gp2PointsDiam];
                Merge=mnl_UniqueRows(Merge);
                MergePoints=Merge(:,1:3);
                MergeDiams=Merge(:,4);
                %Calculate the new colour
                [ChanMean,~,~,~,~,~]=mnl_ExtractKeyVoxels(MergePoints,MergeDiams,Scale,Data);
                Merge_Mean=ChanMean;
                tBkgSub=ChanMean-BkgMean;
                idx=tBkgSub<=0;
                tBkgSub(idx)=0;
                Merge_BkgSub=tBkgSub;
                Merge_NormMean=Merge_BkgSub./MaxValues;
                Merge_VecNorm=mnl_NormaliseVectors(Merge_NormMean);
                %Re-assign Group 1 to include Group 2
                Gp1=Merge_VecNorm;
                Gp1PointsDiam=[];
                Gp1PointsDiam(:,1:3)=MergePoints;
                Gp1PointsDiam(:,4)=MergeDiams;
                %Clear Group 2 and Merge
                Gp2PointsDiam=[];
                clear Merge MergePoints MergeDiams
            else
                %Filter out Group2 points from Group 1
                Gp1PointsDiam=mnl_RemoveGroup2Rows(Gp1PointsDiam,Gp2PointsDiam);
                %Calculate the new Group 1 colour
                [ChanMean,~,~,~,~,~]=mnl_ExtractKeyVoxels(Gp1PointsDiam(:,1:3),Gp1PointsDiam(:,4),Scale,Data);
                Gp1=ChanMean; %Calculate the new colour
                
                %Need to save it Group 1 as a trace
                fTrace(fN).Points=Gp1PointsDiam(:,1:3)./Scale;
                fTrace(fN).OriginalTrace=Trace(i).OriginalTrace;
                fTrace(fN).Diameter=Gp1PointsDiam(:,4);
                fTrace(fN).RealPoints=Gp1PointsDiam(:,1:3);
                
                szPointChecker=size(fTrace(fN).Points,1);
                if szPointChecker>1
                    %Calculate Voxel Length
                    [dist,~]=mnl_CalculateLength(fTrace(fN).Points);
                    fTrace(fN).LengthVx=dist;
                    %Calculate Real Length
                    [dist,~]=mnl_CalculateLength(fTrace(fN).RealPoints);
                    fTrace(fN).LengthReal=dist;
                else
                    fTrace(fN).LengthVx=NaN;
                    fTrace(fN).LengthReal=NaN;
                end
                fTrace(fN).SetId=Trace(i).SetId;
                fTrace(fN).TypeId=Trace(i).TypeId;
                %Now make group 2 into group 1
                Gp1=Gp2;
                Gp1PointsDiam=[];
                Gp1PointsDiam=Gp2PointsDiam(:,1:3);
                Gp1PointsDiam(:,4)=Gp2PointsDiam(:,4);
                %And clear excess variables
                fN=fN+1;
                clear Gp2 Gp2PointsDiam
            end
            j=j+1;
        end
        
        fTrace(fN).Points=Gp1PointsDiam(:,1:3)./Scale;
        fTrace(fN).OriginalTrace=Trace(i).OriginalTrace;
        fTrace(fN).Diameter=Gp1PointsDiam(:,4);
        fTrace(fN).RealPoints=Gp1PointsDiam(:,1:3);
        [dist,~]=mnl_CalculateLength(fTrace(fN).Points);
        fTrace(fN).LengthVx=dist;
        [dist,~]=mnl_CalculateLength(fTrace(fN).RealPoints);
        fTrace(fN).LengthReal=dist;
        fTrace(fN).SetId=Trace(i).SetId;
        fTrace(fN).TypeId=Trace(i).TypeId;
        fN=fN+1;
        clear Distances Gp1PointsDiam Gp2PointsDiam
    else
        fTrace(fN).Points=Trace(i).Points./Scale;
        fTrace(fN).OriginalTrace=Trace(i).OriginalTrace;
        fTrace(fN).Diameter=Trace(i).Diameter;
        fTrace(fN).RealPoints=Trace(i).Points;
        fTrace(fN).LengthVx=Trace(i).LengthVx;
        fTrace(fN).LengthReal=Trace(i).LengthReal;
        fTrace(fN).SetId=Trace(i).SetId;
        fTrace(fN).TypeId=Trace(i).TypeId;
        fN=fN+1;
    end
    clear Group_Mean Group_BkgSub Group_NormMean Group_VecNorm
    mnl_InsertProgressTrackerInLoops(i,szT)
end
end
function [dist,DistanceBetweenPoints]=mnl_CalculateLength(Points)
nPoints=size(Points,1);
for i=1:nPoints-1
    Pos1=Points(i,:);%The xyz position
    Pos2=Points(i+1,:);
    Diff=Pos1-Pos2;
    Sq=Diff.^2;
    SS=sum(Sq);
    DistanceBetweenPoints(i)=sqrt(SS);
end
dist=sum(DistanceBetweenPoints);
end
function [ChanMean,nVoxels,Length_vx,Length,VoxelVals,VoxelsLoc]=mnl_ExtractKeyVoxels(PointsList,DiameterList,Scale,Data)
%Function to extract the voxels of the point list
% Inputs
% - PointsList - Points Extracted from Neurolucida
% - DiameterList - The diameter around each point, from Neurolucida
% - Scale - The size of each voxel [x y z]
% - Data - The image used to extract the colours from
% Outputs
% - ChanMean - Mean Colour
% - nVoxels - number of Voxels
% - Length_vx - The length of the plot in voxels
% - Length - The length in microns
% - VoxelVals - a matrix containing all the colour values of each voxel
% - BW - The mask used to generate the image
%% Step 1 - Calculate the total length of the Fragment
[Length]=mnl_MeasureDistaceOfPoints(PointsList);   
%% Step 2 - Convert the Numbers back to Pixels
VoxList=PointsList./Scale;
[Length_vx]=mnl_MeasureDistaceOfPoints(VoxList);   
%Flip Y and Z because they are negative
XPointList=VoxList(:,1);
YPointList=VoxList(:,2)*-1;
ZPointList=VoxList(:,3)*-1;
Points=[round(XPointList) round(YPointList) round(ZPointList)]; %NB Not flipped yet
%% Step 3 - Find the Key Voxels
dim=size(Data);
xlimit=dim(2);%NB X and Y are flipped in the image
ylimit=dim(1);
nColours=dim(3);
zlimit=dim(4);
numPoints=size(PointsList,1);
n=1;
for j=1:numPoints
    rXY=DiameterList(j); %Extract the diameter of the axon
    r1=round(rXY/Scale(1)/2); %Convert to pixels for xy
    r2=round(rXY/Scale(3)/2); %Convert to pixels for z
    ZPSF=round(2*Scale(3)/2); %Z point spread function in um converted to pixels
    zr=r2+ZPSF; %Z variable in the future shold be variable on the NA of the lens
    r=[r1 r1 zr]; % radius spread in [x y z] pixels
    %Central Positions
    Xpos=Points(j,1);
    Ypos=Points(j,2);
    Zpos=Points(j,3);
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
    %Make a 3D grid
    [x,y,z]=ndgrid(Xmin:Xmax,Ymin:Ymax,Zmin:Zmax);
    tVx=[x(:),y(:),z(:)];
    numVx=size(tVx,1);
    n2=n+numVx-1;
    tAllVoxels(n:n2,:)=tVx;
    clear tVx
    n=n2+1;
end
[VoxelsLoc]=unique(tAllVoxels,'rows'); %
%% Step 4 - Obtain the colour values of that voxel
nPoints=size(VoxelsLoc,1); %The number of points
TempColours=zeros(nPoints,nColours);
for j=1:nPoints
    TempColours(j,:)=Data(VoxelsLoc(j,2),VoxelsLoc(j,1),:,VoxelsLoc(j,3)); %NB cData(y,x,:,z)
end
VoxelVals=TempColours;
nVoxels=nPoints;
%% Step 5 - Calculate average and SD of each channel
ChanMean=mean(VoxelVals,1,'omitnan');
end
function [dist]=mnl_MeasureDistaceOfPoints(PointsList)
nPoints=size(PointsList,1);
dist=0;
for i=1:nPoints-1
    Pos1=PointsList(i,:);
    Pos2=PointsList(i+1,:);
    [D]=mnl_EuclideanDistance(Pos1,Pos2);
    distPoints(i)=D;
    dist=dist+D;   
end
end
function [EuD]=mnl_EuclideanDistance(X,Y)
EuD=sqrt(sum((X-Y).^2));
end
function [Matrix]=mnl_EuclideanDistanceMatrix(Values)
szV=size(Values,1);
for i=1:szV
    Pos1=Values(i,:);
    for j=1:szV
        Pos2=Values(j,:);
        Matrix(i,j)=mnl_EuclideanDistance(Pos1,Pos2);
    end
end
end
function [NewMatrix]=mnl_UniqueRows(Matrix)
numRows=size(Matrix,1);
numCols=size(Matrix,2);
nM=1;
for i=1:numRows
    Row1=Matrix(i,:);
    if isnan(Row1)==0
        for j=i+1:numRows
            Row2=Matrix(j,:);
            if isequal(Row1,Row2)==1
                Matrix(j,1:numCols)=nan(1,numCols);
            end
        end
    end    
end
for i=1:numRows
    if isnan(Matrix(i,:))==0
        NewMatrix(nM,:)=Matrix(i,:);
        nM=nM+1;
    end
end
end
function [NewGroup]=mnl_RemoveGroup2Rows(Group1,Group2)
%Function to remove the group 2 points from the group 1 set
Gp1_NumRows=size(Group1,1);
Gp1_NumCols=size(Group1,2);
Gp2_NumRows=size(Group2,1);
%Search through each group 2 value to find if there are any matches in group1
n=1;
for i=1:Gp1_NumRows
    Row1=Group1(i,:);
    for j=1:Gp2_NumRows
        Row2=Group2(j,:);
        if isequal(Row1,Row2)==1
            Group1(i,1:Gp1_NumCols)=nan(1,Gp1_NumCols);
        end
    end
end
for i=1:Gp1_NumRows
    if isnan(Group1(i,:))==0
        NewGroup(n,:)=Group1(i,:);
        n=n+1;
    end
end
end