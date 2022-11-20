function mnl_ImagesOfClusters(efPxTrace,dim,Clusters,ChosenClusters,ScaleFactor)
%Specify the new dimensions
NewDim=[round(dim(1)*ScaleFactor) round(dim(2)*ScaleFactor)];
%White Background
BW=ones(NewDim(2),NewDim(1),3);
%% Draw in all the traces
%Get all the Points
nTrace=size(efPxTrace,2);
VxList=[];
disp('Drawing in all the traces...')
for i=1:nTrace
    TempVx=efPxTrace(i).AllVoxels(:,1:2);
    VxList=[VxList;TempVx];
    mnl_InsertProgressTrackerInLoops(i,nTrace);
end
VxList=unique(VxList,'rows');
%Now Scale Down
sVxList=round(VxList*ScaleFactor);
sVxList=unique(sVxList,'rows');
%Check that the image doesn't exceed the x and y dimensions
[XExceedList]=find(sVxList(:,1)>NewDim(1));
[YExceedList]=find(sVxList(:,2)>NewDim(2));
if isempty(XExceedList)==0
    sVxList(XExceedList,1)=NewDim(1);
end
if isempty(YExceedList)==0
    sVxList(YExceedList,2)=NewDim(2);
end
%Check if there are zero values
[XMinList]=find(sVxList(:,1)<=0);
[YMinList]=find(sVxList(:,2)<=0);
if isempty(XMinList)==0
    sVxList(XMinList,1)=1;
end
if isempty(YMinList)==0
    sVxList(YMinList,2)=1;
end
%Draw it in
[BW]=mnl_DrawInTrace_2D(BW,sVxList,[0.75 0.75 0.75]);
%% Now draw in the chosen traces
nClust=size(ChosenClusters,2);
Nz=length(num2str(nClust));
Clust_fmt=num2str(Nz,'%%0%dg');

for i=1:nClust
    ClustId=ChosenClusters(i);
    dn=sprintf('%s%d','Drawing in Cluster ',ClustId);
    tClustId=ClustId;
    disp(dn)
    %Does this cluster exist?
    if isempty(Clusters(ClustId).Centroid)~=1
        tClustId=sprintf(Clust_fmt,i);
        fn=sprintf('%s%s%s','Filtered Cluster ',tClustId,'_Scaled.tiff');
        %Find which traces to add        
        TraceList=Clusters(ClustId).Traces;
        NumTraces=size(TraceList,1);
        TraceVx=[];
        for j=1:NumTraces
            Temp=efPxTrace(TraceList(j)).AllVoxels(:,1:2);
            TraceVx=[TraceVx;Temp];
        end
        TraceVx=unique(TraceVx,'rows');
        %Now Scale Down
        sTraceVxList=round(TraceVx*ScaleFactor);
        sTraceVxList=unique(sTraceVxList,'rows');
        %Check that the image doesn't exceed the x and y dimensions
        [XExceedList]=find(sTraceVxList(:,1)>NewDim(1));
        [YExceedList]=find(sTraceVxList(:,2)>NewDim(2));
        if isempty(XExceedList)==0
            sVxList(XExceedList,1)=NewDim(1);
        end
        if isempty(YExceedList)==0
            sVxList(YExceedList,2)=NewDim(2);
        end
        %Check if there are zero values
        [XMinList]=find(sTraceVxList(:,1)<=0);
        [YMinList]=find(sTraceVxList(:,2)<=0);
        if isempty(XMinList)==0
            sTraceVxList(XMinList,1)=1;
        end
        if isempty(YMinList)==0
            sTraceVxList(YMinList,2)=1;
        end
        %Now draw it in
        tBW=mnl_DrawInTrace_2D(BW,sTraceVxList,[0 1 0]); %Green Trace
        imwrite(tBW,fn)
    end
end
end
function [BW]=mnl_DrawInTrace_2D(BW,VxList,Colour)
%Draw in your chosen voxels
% Inputs
%  BW - BW matrix of image
%  VxList - The list of voxels to draw in
%  Colour - The chosen colour (e.g. [1 0 0]=Red)
%
% Outputs
%  BW - Updated image
tVxList=VxList(:,1:2); %ignore the z levels
fVxList=unique(tVxList,'rows'); %flatted along Z
szL=size(fVxList,1);
for i=1:szL
    Pos=fVxList(i,:);
    BW(Pos(2),Pos(1),:)=Colour;
end
end
