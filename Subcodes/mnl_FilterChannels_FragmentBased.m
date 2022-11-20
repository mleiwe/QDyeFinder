function [BadChannels,SampleBkg,PxTrace,BkgSigCumulativeDistribution]=mnl_FilterChannels_FragmentBased(Trace,Scale,Data)
% Function to automatically or manually filter out the bad channels. It
% will also save a MATLAB figure and png of the signal to noise
% distributions of each channel
%
% Created by Marcus Leiwe, Kyushu Univiersity, 2021
%
% Inputs
%  Trace - the Structured Trace extracted from "mnl_ReadXML_BranchesSeperately_WhileLoop"
%  Scale - the Scale of the image
%  Data - the chromatically corrected image data
%
% Outputs
%  Bad Channels - a list of the bad channels
%  SampleBkg - A chosen sample of background voxels
%  PxTrace - Structure for the Neurolucida traces
%  BkgSigCumulativeDistribution - Structure containing the cumulative distributions for the signal and noise per channel

%% Extract the Raw Colour Values Of the Traces
szT=size(Trace,2);
dim=size(Data);
n=1; %counter for collating all the voxels
PxTrace=struct('Points',[],'Diameter',[],'OriginalTrace',[],'SetId',[],'LengthVx',[],'LengthReal',[],'ChanMean',[]);
for i=1:szT %Do this for each trace
    %% Extract the Mask for the Whole Fragment
    AllPointsList=Trace(i).Points; %Extract the Points Desired
    %AllDiameterList=Trace(i).Diameter; %Extract the Diameters Desired
    AllDiameterList=Trace(i).Diameter./Trace(i).Diameter; %Extract the Diameters Desired %Filtered to widen the width
    [ChanMean,~,LengthVx,Length,VoxelVals,VoxelLocs]=mnl_ExtractKeyVoxels(AllPointsList,AllDiameterList,Scale,Data);
    PxTrace(i).Points=Trace(i).Points./Scale;
    PxTrace(i).Diameter=Trace(i).Diameter;
    PxTrace(i).OriginalTrace=Trace(i).OriginalTrace;
    PxTrace(i).SetId=Trace(i).SetId;
    PxTrace(i).LengthVx=LengthVx;
    PxTrace(i).LengthReal=Length;
    PxTrace(i).ChanMean=ChanMean;
    %Collate all the values
    szV=size(VoxelVals,1);
    AllVx(n:n+szV-1,:)=VoxelLocs;
    n=n+szV;
    clear AllPointsList AllDiameterList VoxelLocs VoxelVals
end
%% Calculate the Background Values
%Find the Voxels of Interest
AllVx=unique(AllVx,'rows'); %Remove Duplicates
nVx=size(AllVx,1); %Final Number of voxels
% Now the background locations
prompt='Do you want to manually select the background? (y/n)';
UserInput=input(prompt,'s');
ManBkg=strcmp(UserInput,'y');
if ManBkg==1
    %Make a MIP from the Data
    MIP=zeros(dim(1),dim(2));
    for i=1:dim(1)
        for j=1:dim(2)
            a(1:dim(3),1:dim(4))=Data(i,j,:,:);
            MIP(i,j)=max(a(:));
        end
    end    
    %Draw an ROI around the background of the MIP
    figure('Name','Please Draw an ROI of the background')
    imagesc(MIP)
    title('Grayscale MIP')
    [~,~,Bkg_BW,~,~]=roipoly; %Get the BW mask for the background ROI
    %Now assign the locations to an AllBkg list
    [r,c]=find(Bkg_BW==1);
    manBkg=size(r,1);
    nBkg=manBkg*dim(4);
    Chosen_Bkg=nan(nBkg,3);
    Pos=1;
    for i=1:manBkg
        tempBkgList=ones(dim(4),3);
        temp=ones(dim(4),1);
        tempBkgList(:,1)=temp.*r(i);
        tempBkgList(:,2)=temp.*c(i);
        tempBkgList(:,3)=1:1:dim(4);
        Chosen_Bkg(Pos:(Pos-1+dim(4)),:)=tempBkgList;
        Pos=Pos+dim(4);
    end
    %SampleN_Bkg=nBkg-1;
    SampleN_Bkg=size(Chosen_Bkg,1);
else
    disp('Finding the Background...')
    TotVx=dim(1)*dim(2)*dim(4);%Total Number of Voxels of Interest
    nBkgVx=TotVx-nVx; %Number of background voxels
    [SampleN_Bkg]=mnl_DetermineSampleSize(97.5,nBkgVx,1); %Backgroud Subsampling
    xyzdim=[dim(2),dim(1),dim(4)];
    [Chosen_Bkg]=mnl_SubSampleBkgFromMask(round(SampleN_Bkg),xyzdim,AllVx);
end

%% Now Compare the channels
% For the Signal - create a list of the mean fragment colours for each channel
Signal=nan(szT,dim(3));
for i=1:szT
    Signal(i,:)=PxTrace(i).ChanMean;
end
% For the Background - NB for the background it is quicker to assign it directly, but I used the loop because the memory won't allocate more than
% ~1 million values
disp('Now Calculating the Background Colour')
Background=nan((round(SampleN_Bkg)),dim(3));%Preallocate background for speed
SampleBkg=nan((round(SampleN_Bkg)),3);%Preallocate background for speed
for i=1:round(SampleN_Bkg)
    Pos=Chosen_Bkg(i,:);
    SampleBkg(i,:)=Chosen_Bkg(i,:);
    temp=Data(Pos(2),Pos(1),:,Pos(3));%Remember it is flipped in the data
    Background(i,:)=temp;
    mnl_InsertProgressTrackerInLoops(i,round(SampleN_Bkg))
end
% Plot the differences
figure('Name','Scatter Plots to look for correlation between channels')
c=1;
for i=1:dim(3)
    for j=1:dim(3)
        subplot(dim(3),dim(3),c)
        scatter(Background(:,j),Background(:,i),'.k')
        hold on
        scatter(Signal(:,j),Signal(:,i),'.r')
        c=c+1;
        %Add a legend to the final plot
        if i==dim(3) && j==dim(3)
            legnames={'Background','Signal'};
            legend(legnames,'Location','southeast')
        end
    end
end
% So now we can evaluate the quality of the signal
prompt='Do you want to filter the channels manually (y/n)';
UserInput=input(prompt,'s');
Manual=strcmp(UserInput,'y');
cGood=1;
cBad=1;
Sig2NoiseThresh=3;
BadChannels=[];
GoodChannels=[];
Sig2NoisePeak=zeros(1,dim(3));
BkgSigCumulativeDistribution=struct('BackgroundPercentiles',[],'SignalPercentiles',[]);
for i=1:dim(3) %for each channel
    [Sig2NoisePeak(i),NoisePercentiles,SignalPercentiles]=mnl_EvaluateCumulativeSignal2Noise5(Background(:,i),Signal(:,i),i); %The signal to noise peak anywhere between the 81st to 100th percentile
    BkgSigCumulativeDistribution(i).BackgroundPercentiles=NoisePercentiles;
    BkgSigCumulativeDistribution(i).SignalPercentiles=SignalPercentiles;
    if Manual==1
        txt=sprintf('%s%d%s','Do you want to remove channel ',i,'?, y/n');
        prompt=txt;
        a=input(prompt,'s');
        Bd=strcmp(a,'y');
        Gd=strcmp(a,'n');
        if Bd==1
            BadChannels(cBad)=i;
            cBad=cBad+1;
        elseif Gd==1
            GoodChannels(cGood)=i;
            cGood=cGood+1;
        else
            disp('Wrong Input Detected - one chance to correct it')
            txt=sprintf('%s%d%s','Do you want to remove channel ',i,'?, y/n');
            prompt=txt;
            a=input(prompt,'s');
            Bd=strcmp(a,'y');
            Gd=strcmp(a,'n');
            if Bd==1
                BadChannels(cBad)=i;
                cBad=cBad+1;
            elseif Gd==1
                GoodChannels(cGood)=i;
                cGood=cGood+1;
            end
        end
    else
        %Automated Detection
        if Sig2NoisePeak(i)>=Sig2NoiseThresh
            GoodChannels(cGood)=i;
            cGood=cGood+1;
        else
            BadChannels(cBad)=i;
            cBad=cBad+1;
        end
    end 
end
h=figure('Name','All the Signal to Noise Peaks');
bar(Sig2NoisePeak);
hold on
x=0:dim(3)+1;
y=ones(1,dim(3)+2)*Sig2NoiseThresh;
p=plot(x,y,'--r');
xlabel('Channel Number')
ylabel('Signal to Noise')
title('Peak Signal to Noise Ratio')
pleg=sprintf('%s%s','Signal to Noise Threshold - ',num2str(round(Sig2NoiseThresh,2)));
legend(p,pleg);
savefig(h,'Sig2Noise_AllPeaks');
saveas(h,'Sig2Noise_AllPeaks','png');
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
    if isnan(rXY)==1
        rXY=nanmean(DiameterList);
        if isnan(rXY)==1
            rXY=1;
        end
    end
    r1=round(rXY/Scale(1)/2); %Convert to pixels for xy
    r2=round(rXY/Scale(3)/2); %Convert to pixels for z
    ZPSF=round(2*Scale(3)/2); %Z point spread function in um converted to pixels
    zr=r2+ZPSF; %Z variable in the future shold be variable on the NA of the lens
    r=[r1 r1 zr]; % radius spread in [x y z] pixels
    %Central Positions
    Xpos=Points(j,1)+1; %Add one to compensate that MATLAB starts at 1 while the image measurements start at 0
    Ypos=Points(j,2)+1; %Add one to compensate that MATLAB starts at 1 while the image measurements start at 0
    Zpos=Points(j,3)+1; %Add one to compensate that MATLAB starts at 1 while the image measurements start at 0
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
clear TempColours
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