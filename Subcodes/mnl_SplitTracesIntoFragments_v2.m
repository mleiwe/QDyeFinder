function [FragmentValues,BkgMean,MaxValues]=mnl_SplitTracesIntoFragments_v2(Trace,Scale,Data,SampleBkgVxList)
% Function that splits every trace into fragments of all possible lengths
% based on the number of points
%
%Inputs
% Trace - List of the Traces
% Scale - The scale of the image (Data)
% Data - The image to get the raw colours from
% SampleBkgVxList - The sample list of voxels to represent the background
%
%Outputs
% FragmentValues - Calculate the mean color of each sub fragment and
% measure the distance to the ground truth. Also calculate the distance,
% number of voxels recruited, and the brightness magnitude
%
%Marcus Leiwe, Kyushu University - 24th December 2020
    
%% Extract the Raw Colour Values so we can Normalise
szT=size(Trace,2);
dim=size(Data);
n=1; %counter for collating all the voxels
FragmentValues=struct('All',[],'FragSize',[]);
for i=1:szT %Do this for each trace
    %% Extract the Mask for the Whole Fragment
    AllPointsList=Trace(i).Points; %Extract the Points Desired
    AllDiameterList=Trace(i).Diameter; %Extract the Diameters Desired
    [AllChanMean,nVoxels,Length_vx,Length,VoxelVals,VoxelLocs]=mnl_ExtractKeyVoxels(AllPointsList,AllDiameterList,Scale,Data);
    %Save it to the structure
    FragmentValues(i).All.Raw.Mean=AllChanMean;
    FragmentValues(i).All.Raw.VoxelVals=VoxelVals;
    FragmentValues(i).All.nVoxels=nVoxels;
    FragmentValues(i).All.Length_vx=Length_vx;
    FragmentValues(i).All.Length=Length;
    FragmentValues(i).All.AllVoxels=VoxelLocs;
    FragmentValues(i).All.Points=Trace(i).Points;
    %Collate all the values
    szV=size(VoxelVals,1);
%     AllVoxelValues(n:n+szV-1,:)=VoxelVals;
%     AllVx(n:n+szV-1,:)=VoxelLocs;
    n=n+szV;
    clear AllPointsList AllDiameterList VoxelLocs VoxelVals
end
%% Calculate the values
% For the Background
disp('Calculating the Background Mean and Maximum Trace Values...')
NumBkgVx=size(SampleBkgVxList,1);
Background=nan(NumBkgVx,dim(3));
for i=1:NumBkgVx
    Pos=SampleBkgVxList(i,:);
    temp=Data(Pos(2),Pos(1),:,Pos(3));
    Background(i,:)=temp;
end
%Calculate the background mean from the sample
mBkg=median(Background,'omitnan');
BkgMean=mBkg; %This value will be used a lot
%Then subtract this from the mean trace value
tMatrix=nan(szT,dim(3));
for i=1:szT
    tVal=FragmentValues(i).All.Raw.Mean-BkgMean;
    idx_subzero=find(tVal<0);
    if isempty(idx_subzero)==0
        tVal(idx_subzero)=0;
    end
    FragmentValues(i).All.BkgRemovedMean=tVal;
    tMatrix(i,:)=FragmentValues(i).All.BkgRemovedMean;
end
%Find the maximum average value
MaxValues=max(tMatrix,[],'omitnan');
disp('Processing All Traces...')
%Normalise to the max average value - also create a temporary Matrix of all the traces, (and their sub-fragments?) 
tMatrix=nan(szT,dim(3));
for i=1:szT
    NormValues=FragmentValues(i).All.BkgRemovedMean./MaxValues;
    FragmentValues(i).All.NormMean=NormValues;
    tMatrix(i,:)=NormValues;
    FragmentValues(i).All.NormMeanMagnitude=vectorNorm(NormValues);
end
%Then Vector Normalise
for i=1:szT
    FragmentValues(i).All.VecNormMean=mnl_SingleVectorNormalise(tMatrix(i,:));
end

%% Now measure the spread of colours for traces and fragments
disp('Fragmenting Traces...')
% Set a limit to stop over processing
NumTraces=size(Trace,2);
%Now begin measuring the spread
for i=1:NumTraces
    nPoints=size(Trace(i).Points,1); %Total Number of pointsi
    for j=1:nPoints
        FragmentValues(i).FragSize(j).PointNum=j;
        % Calculate the all the potential groups of points
        nGroups=floor(nPoints/j);%Total number of potential groups
        if nGroups==1
            st=1;ed=j;
            PointList=Trace(i).Points(st:ed,:); %Extract the Points Desired
            DiameterList=Trace(i).Diameter(st:ed); %Extract the Diameters Desired
            [ChanMean,nVoxels,Length_vx,Length,~,VoxelLocs]=mnl_ExtractKeyVoxels(PointList,DiameterList,Scale,Data);            
            %Raw Information
            FragmentValues(i).FragSize(j).nVoxels=nVoxels; %The number of voxels
            FragmentValues(i).FragSize(j).Length_vx=Length_vx; %The length of the fragment in voxels
            FragmentValues(i).FragSize(j).Length=Length; %The length of the fragment in microns
            FragmentValues(i).FragSize(j).NormalisedLength=Length/Trace(i).LengthReal;
            FragmentValues(i).FragSize(j).RawMean(1,:)=ChanMean; %The mean colour of the fragment
            FragmentValues(i).FragSize(j).VoxelLocs(1).AllVoxels=VoxelLocs;
            %Subtract the Mean Background
            BkgRemovedMean=ChanMean-BkgMean;
            FragmentValues(i).FragSize(j).BkgRemovedMean=BkgRemovedMean; % Subtracting the mean background value
            %Normalise to the maximum value
            NormValues=BkgRemovedMean./MaxValues;
            FragmentValues(i).FragSize(j).NormMean=NormValues; %The normalised mean colour of the fragment
            %Calculate the Vector magnitude
            FragmentValues(i).FragSize(j).NormMeanMagnitude=vectorNorm(NormValues);
            %Vector Normalise with its Ground Truth
            VecNormTruth=FragmentValues(i).All.VecNormMean; 
            [VecNormFrag]=mnl_SingleVectorNormalise(FragmentValues(i).FragSize(j).NormMean(1,:));
            FragmentValues(i).FragSize(j).VecNormMean(1,:)=VecNormFrag;
            FragmentValues(i).FragSize(j).VecNormTruth(1,:)=VecNormTruth;
            % Calculate Distance to the True Colour
            [EuD]=mnl_EuclideanDistance(VecNormTruth,VecNormFrag);
            FragmentValues(i).FragSize(j).EuDistance=EuD; %Euclidean Distance
            clear PointsList DiameterList
        elseif nGroups>1 %If there is a possibility to generate more than one group
            %Firstly subsample to save time
            Zpercent=95; %confidence interval
            N=nGroups; %population size
            e=1; %margin of error in percent
            TempSampleSize=mnl_DetermineSampleSize(Zpercent,N,e);
            TempSampleSize=round(TempSampleSize);
            if TempSampleSize<1
                TempSampleSize=1;
            end
            GroupIDs=nan(nGroups,2); %Preallocate a matrix that has the start (column 1) and end(column 2) positions of each group
            ed=0;
            for k=1:nGroups
                st=ed+1;
                ed=st+j-1;
                GroupIDs(k,:)=[st ed];
            end
            SubGroupIDs=randperm(nGroups,TempSampleSize);
            %Ok now do it for each group
            for k=1:TempSampleSize
                id=SubGroupIDs(k);
                st=GroupIDs(id,1);
                ed=GroupIDs(id,2);
                PointList=Trace(i).Points(st:ed,:); %Extract the Points Desired                
                DiameterList=Trace(i).Diameter(st:ed); %Extract the Diameters Desired
                [ChanMean,nVoxels,Length_vx,Length,~,VoxelLocs]=mnl_ExtractKeyVoxels(PointList,DiameterList,Scale,Data);
                FragmentValues(i).FragSize(j).nVoxels(k)=nVoxels; %The number of voxels
                FragmentValues(i).FragSize(j).Length_vx(k)=Length_vx; %The length of the fragment in voxels
                FragmentValues(i).FragSize(j).Length(k)=Length; %The length of the fragment in microns
                FragmentValues(i).FragSize(j).NormalisedLength(k)=Length/Trace(i).LengthReal;
                FragmentValues(i).FragSize(j).RawMean(k,:)=ChanMean; %The mean colour of the fragment
                FragmentValues(i).FragSize(j).VoxelLocs(k).AllVoxels=VoxelLocs;
                %Subtract the Mean Background
                BkgRemovedMean=ChanMean-BkgMean;
                FragmentValues(i).FragSize(j).BkgRemovedMean(k,:)=BkgRemovedMean; % Subtracting the mean background value
                %Normalise to the maximum value
                NormValues=BkgRemovedMean./MaxValues;
                FragmentValues(i).FragSize(j).NormMean(k,:)=NormValues; %The normalised mean colour of the fragment
                %Measure the brightness
                FragmentValues(i).FragSize(j).NormMeanMagnitude(k)=vectorNorm(NormValues);
                %Vector Normalise with its Ground Truth
                VecNormTruth=FragmentValues(i).All.VecNormMean;
                [VecNormFrag]=mnl_SingleVectorNormalise(FragmentValues(i).FragSize(j).NormMean(k,:));
                FragmentValues(i).FragSize(j).VecNormMean(k,:)=VecNormFrag;
                FragmentValues(i).FragSize(j).VecNormTruth(k,:)=VecNormTruth;
                % Calculate Distance to the True Colour
                [EuD]=mnl_EuclideanDistance(VecNormTruth,VecNormFrag);
                FragmentValues(i).FragSize(j).EuDistance(k)=EuD; %Euclidean Distance
                clear PointsList DiameterList                
            end
        end
    end
    fprintf('%s%d%s%d\n','Completed Trace ',i,' out of ',NumTraces)
end
end
%% Subfunctions
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
    Xpos=Points(j,1)+1;
    Ypos=Points(j,2)+1;
    Zpos=Points(j,3)+1;
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
function [NormVec]=mnl_SingleVectorNormalise(Vector)
Mag=vectorNorm(Vector);
NormVec=Vector./Mag;
end
function [EuD]=mnl_EuclideanDistance(X,Y)
EuD=sqrt(sum((X-Y).^2));
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