function [cTrace]=mnl_AdjustTraceLocations(Trace,Data,scale)
%This function is designed to check the location of the tracing as
%occaisonally the point locations from the neurolucida detections my not be
%accurate. So we search the nearby area to check that the point is along
%the brightest point along the y axis (rows of Data)
% Outputs
% cTrace - The corrected traces locations
%
% Inputs 
% Trace - The original neurolucida trace points
% Data - The 4d image
% scale - The scale of the image, used to convert the locations to voxel locations

%% Get the key information
nTraces=size(Trace,2);
dim=size(Data);
%Data=double(Data);
%Get the values to normalise to
nVal=nan(1,dim(3));
for i=1:dim(3)
    temp(:,:,:)=Data(:,:,i,:);
    nVal(i)=prctile(temp(:),99.95);
end
%% Now correct the traces
cTrace=Trace;
for i=1:nTraces
    %How many points
    nPoints=size(Trace(i).Points,1);
    for j=1:nPoints
        Loc=Trace(i).Points(j,:);
        Diam=abs(Trace(i).Diameter(j)); %Made sure that the value will be positive
        %Convert to voxels
        VxLoc=Loc./scale;
        VxDiam=[Diam Diam Diam]./scale;
        Regions=[1 VxDiam(2)*10 2];
        %Flip Y and Z axis because they will be negative
        VxLoc=VxLoc.*[1 -1 -1];
        %Get the small image area
        [tData]=nested_Mnl_ExtractImageLocation(VxLoc,Regions,Data,dim);
        %Convert to grayscale
        sztD=size(tData);
        nData=zeros(sztD(1),sztD(2),sztD(3),sztD(4));
        for ci=1:dim(3)
            temp=tData(:,:,ci,:);
            nData(:,:,ci,:)=temp/nVal(ci);
        end
        idx=nData>1;
        nData(idx)=1;
        gsData(:,:,:)=max(nData,[],3); %Find the maximum for each channel
        %Compress along Y axis (rows)
        for yi=1:sztD(1)
            XZ(:,:)=gsData(yi,:,:);
            Yvals(yi)=max(XZ(:));
        end
        clear XZ
        %Smooth the trace via median
        MedFilt=3;
        for yi=1:sztD(1)
            minY=yi-MedFilt;
            if minY<1
                minY=1;
            end
            maxY=yi+MedFilt;
            if maxY>sztD(1)
                maxY=sztD(1);
            end
            mYvals(yi)=median(Yvals(minY:maxY));
        end
        %Now find the centre of mass
        [CoM]=Nested_mnl_FindCoM(mYvals);

        %Find the new position
        InitialPos=sztD(1)/2;
        DifferenceVx=CoM-InitialPos;
        Difference=DifferenceVx*scale(1);
        Difference=Difference*-1;%Remember to flip it back to negative
        newY=Loc(2)+Difference;
        nLoc=[Loc(1) newY Loc(3)];
        cTrace(i).Points(j,:)=nLoc;

        clear gsData Yvals mYvals
        
    end
    mnl_InsertProgressTrackerInLoops(i,nTraces)
end
end

function [tData]=nested_Mnl_ExtractImageLocation(VxLoc,VxDiam,Data,dim)
%Function to extract the area surrounding the image
%Calculate the min and max values
Xmin=round(VxLoc(1)-VxDiam(1));
if Xmin<1
    Xmin=1;
end
Xmax=round(VxLoc(1)+VxDiam(1));
if Xmax>dim(2)
    Xmax=dim(2);
end
Ymin=round(VxLoc(2)-VxDiam(2));
if Ymin<1
    Ymin=1;
end
Ymax=round(VxLoc(2)+VxDiam(2));
if Ymax>dim(1)
    Ymax=dim(1);
end
Zmin=round(VxLoc(3)-VxDiam(3));
if Zmin<1
    Zmin=1;
end
Zmax=round(VxLoc(3)+VxDiam(3));
if Zmax>dim(4)
    Zmax=dim(4);
end
tData(:,:,:,:)=double(Data(Ymin:Ymax,Xmin:Xmax,:,Zmin:Zmax));
end
function [CoM]=Nested_mnl_FindCoM(Yvals)
TotalMass=sum(Yvals);
nPoints=size(Yvals,2);
Val=nan(nPoints,1);
for i=1:nPoints
    Val(i)=i*Yvals(i);
end
CoM=sum(Val)/TotalMass;
end