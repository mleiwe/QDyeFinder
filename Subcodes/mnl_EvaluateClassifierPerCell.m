function [Kappa,F1score,Recall,Precision,TruePositiveRate,FalsePositiveRate]=mnl_EvaluateClassifierPerCell(efPxTrace,ClusterIDs,Clusters,dim,DisplayFig,EuThresh,NumTraceLim)
%NB now excludes cells which have fewer traces than the minimum trace
nTrace=size(efPxTrace,2);
EuThresh=num2str(round(EuThresh,2));
%% Create the white for all the dendrites
% AllTraces=[];
% for i=1:nTrace
%     tTrace=efPxTrace(i).AllVoxels(:,1:2);
%     AllTraces=[AllTraces;tTrace];
% end
% AllTraces=unique(AllTraces,'rows');
% [BW]=mnl_DrawInTrace_2D(BW,AllTraces,[1 1 1]);
%% So now I want to trace each manual one and see how accurate the trace is
% Find out the name of each set Id
IdList={};
SetCounter=1;
for i=1:nTrace
    if isempty(efPxTrace(i).SetId)~=1
        NewSet=1;
        tSetId=efPxTrace(i).SetId;
        %Is it a new Set?
        if isempty(IdList)==1
            NewSet=1;
        else
            szId=size(IdList,2);
            for j=1:szId
                if strcmp(IdList{j},tSetId)==1
                    NewSet=0;
                end
            end
        end
        if NewSet==1
            IdList{SetCounter}=tSetId;
            SetCounter=SetCounter+1;
        end
    end
end
%Create a structure with all the required traces                
SetIds=struct('SetName',[],'Traces',[]);
nSets=SetCounter-1;
for i=1:nSets
    SetName=IdList{i};
    TraceList=[];
    TraceCounter=1;
    for j=1:nTrace
        tSetId=efPxTrace(j).SetId;
        if isempty(tSetId)==0
            if strcmp(tSetId,SetName)==1
                TraceList(TraceCounter,1)=j;
                TraceCounter=TraceCounter+1;
            end
        end
    end
    %Now Add to the structure
    SetIds(i).SetName=SetName;
    SetIds(i).Traces=TraceList;
end

%% Now evaluate each manually drawn neuron (Set)
if strcmp(DisplayFig,'y')==1
    ScaleFactor=0.25;
    NewDim=[round(dim(1)*ScaleFactor) round(dim(2)*ScaleFactor)];
    % Draw all the traces in light grey
    VxList=[];
    disp('Drawing in all the traces...')
    for i=1:size(efPxTrace,2)
        tVxList=efPxTrace(i).AllVoxels(:,1:2);
        VxList=[VxList;tVxList];
        mnl_InsertProgressTrackerInLoops(i,nTrace);
    end
    %Scale Down
    [sVxList]=mnl_NestedScalePixels(VxList,ScaleFactor,NewDim);
    %Create the white image
    BW=ones(NewDim(2),NewDim(1),3);
    %Draw the all the traces as gray
    [BW]=mnl_DrawInTrace_2D(BW,sVxList,[0.75 0.75 0.75]);
    %Now draw for each cell
    for i=1:nSets
        nTrace=size(SetIds(i).Traces,1);
        VxList=[];
        for j=1:nTrace
            num=SetIds(i).Traces(j);
            tVxList=efPxTrace(num).AllVoxels(:,1:2);
            VxList=[VxList;tVxList];
        end
        [sVxList]=mnl_NestedScalePixels(VxList,ScaleFactor,NewDim);
        [tempBW]=mnl_DrawInTrace_2D(BW,sVxList,[0 0 0]); %Draw the correct traces as black
        clear VxList tVxList sVxList
        %Find the closest matching cluster
        ChosenTraces=SetIds(i).Traces;
        AssociatedClusters=ClusterIDs(ChosenTraces);
        BestCluster(i)=mode(AssociatedClusters);
        if BestCluster(i)==0 %If most of the traces are unassigned, then find the next best option
            NonZeros=find(AssociatedClusters~=0);
            fAssociatedClusters=AssociatedClusters(NonZeros);
            BestCluster(i)=mode(fAssociatedClusters);
        end
        [idx]=find(AssociatedClusters==BestCluster(i));
        NumClusterTraces=size(idx,1);
        ClusterAccuracy(i)=(NumClusterTraces/nTrace)*100;        
        %Now Plot all of the chosen cluster as red (=false positive, but will overwrite the correct ones to green)
        if ~isnan(BestCluster(i)) %If there is a best cluster
            Ctraces=Clusters(BestCluster(i)).Traces;
            SetTraces=SetIds(i).Traces(j);
            nCtraces=size(Clusters(BestCluster(i)).Traces,1);       
            VxList=[];
            Fp=0;%False positive counter
            for j=1:nCtraces
                num=Clusters(BestCluster(i)).Traces(j);
                tVxList=efPxTrace(num).AllVoxels(:,1:2);
                VxList=[VxList;tVxList];
                if strcmp(efPxTrace(num).SetId,SetIds(i).SetName)==0
                    Fp=Fp+1;
                end
            end
            [sVxList]=mnl_NestedScalePixels(VxList,ScaleFactor,NewDim);
            [tempBW]=mnl_DrawInTrace_2D(tempBW,sVxList,[1 0 0]);
            %Now Plot the correct ones as green
            VxList=[];
            %Find the correct traces
            tp=0;
            Fn=0;
            for j=1:nTrace
                SetTrace=SetIds(i).Traces(j);
                if ClusterIDs(SetTrace)==BestCluster(i)
                    tp=tp+1;
                    TpTrace(tp)=SetTrace;
                else
                    Fn=Fn+1;
                end
            end
            for j=1:tp
                num=TpTrace(j);
                tVxList=efPxTrace(num).AllVoxels(:,1:2);
                VxList=[VxList;tVxList];
            end
            [sVxList]=mnl_NestedScalePixels(VxList,ScaleFactor,NewDim);
            [tempBW]=mnl_DrawInTrace_2D(tempBW,sVxList,[0 1 0]);
        end
        fn=sprintf('%s%s%d%s%s%s',IdList{i},' - BestCluster',BestCluster(i),'EuThresh_',EuThresh,'.tiff');
        fn2=sprintf('%s%s%d%s',IdList{i},' - BestCluster',BestCluster(i),'EuThresh_',EuThresh);
        fn3=sprintf('%s%s%d',IdList{i},' - BestCluster',BestCluster(i));
        Htemp=figure;
        subplot(1,2,1)
        for j=1:nTrace
            num=SetIds(i).Traces(j);
            TracePoints=efPxTrace(num).Points;
            Diameters=efPxTrace(num).Diameter;
            mnl_PlotSingleTrace2D(TracePoints,Diameters,[0 0 1])           
        end
        axis equal
        axis ij
        xlim([0 dim(1)])
        ylim([0 dim(2)])
        title(IdList{i})
        grid on
        subplot(1,2,2)
        image(tempBW)
        szBW=size(tempBW);
        axis equal
        xlim([0 szBW(2)])
        ylim([0 szBW(1)])
        axis off
        
        %Calculate the F1 Score
        tRecall=tp/(tp+Fn);
        tPrecision=tp/(tp+Fp);
        tF1score=2*((tPrecision*tRecall)/(tPrecision+tRecall));
        tF1=num2str(round(tF1score,3));
        %Now add to title
        title(sprintf('%s%d%s%s','Cluster ',BestCluster(i),' F1 score - ',tF1))

        savefig(Htemp,sprintf('%s%s',fn2,'.fig'))
        mnl_ExportEPSdense(Htemp,fn3)
        imwrite(tempBW,fn)
        close(Htemp)
    end        
end
%% Create a matrix of scores (True Positive, True Negative, False,Positive, and False Negative (NB Only Manually scored ones)
ConfusionMatrices=[];%[TruePos FalseNeg FalsePos TrueNeg]
cK=1;
AllTrace=size(efPxTrace,2);
for i=1:nSets
    nTrace=size(SetIds(i).Traces,1);
    if nTrace>=NumTraceLim
        ChosenTraces=SetIds(i).Traces;
        AssociatedClusters=ClusterIDs(ChosenTraces);
        SetIds(i).BestCluster=mode(AssociatedClusters);
        if SetIds(i).BestCluster==0 %If most of the traces are unassigned, then find the next best option
            NonZeros=find(AssociatedClusters~=0);
            fAssociatedClusters=AssociatedClusters(NonZeros);
            SetIds(i).BestCluster=mode(fAssociatedClusters);
        end
        [idx]=find(AssociatedClusters==SetIds(i).BestCluster);
        ClustTraces=find(ClusterIDs==SetIds(i).BestCluster);
        ClustSize=size(ClustTraces,1);
        %True Positives
        Tp=size(idx,1);
        SetIds(i).TruePositive=Tp;
        %False Negatives
        Fn=nTrace-Tp;
        SetIds(i).FalseNegative=Fn;
        %False Positives
        Fp=0;
        for j=1:ClustSize
            ClusterSetId=efPxTrace(ClustTraces(j)).SetId;
            if strcmp(ClusterSetId,SetIds(i).SetName)==0
                Fp=Fp+1;
            end
        end
        SetIds(i).FalsePositives=Fp;
        %True Negatives
        Tn=AllTrace-Tp-Fn-Fp;
        SetIds(i).TrueNegatives=Tn;
        SetIds(i).ClusterAccuracy=(ClustSize/nTrace)*100;
        %Confusion Matrix
        ConfusionMatrices(i,:)=[Tp Fn Fp Tn];%[TruePos FalseNeg FalsePos TrueNeg]
    else
        SetIds(i).BestCluster=NaN;
        SetIds(i).TruePositive=NaN;
        SetIds(i).FalseNegative=NaN;
        %False Positives
        SetIds(i).FalsePositives=NaN;
        %True Negatives
        SetIds(i).TrueNegatives=NaN;
        SetIds(i).ClusterAccuracy=NaN;
        %Confusion Matrix
        ConfusionMatrices(i,:)=nan(1,4);%[TruePos FalseNeg FalsePos TrueNeg]
    end
end
%% Trace Distribution
PercentTraces=zeros(nSets,size(Clusters,2));
nClust=size(Clusters,2);
for i=1:nSets
    SelectedTraces=SetIds(i).Traces;
    TotTraces=size(SelectedTraces,1);
    SelectedClusterIDs=ClusterIDs(SelectedTraces);

    [ClusterCounts,ClusterNames]=groupcounts(SelectedClusterIDs);

    nDiffClust=size(ClusterNames,1);
    for j=1:nDiffClust
        Percent=(ClusterCounts(j)/TotTraces)*100;
        PercentTraces(i,ClusterNames(j))=Percent;
    end   
end
if strcmp(DisplayFig,'y')==1
    figure('Name','Distribution of Traces by Cluster')
    imagesc(PercentTraces,[0 100])
    axis equal
    colormap(magma)
    colorbar
    xlabel('Cluster Number')
    ylabel('Set Number')
    ylim([0.5 nSets+0.5])
    xlim([0.5 nClust+0.5])
end
%% Now Calculate Cohen's Kappa per Trace
Kappa=[];
for i=1:nSets
    if isnan(ConfusionMatrices(i,:))==0
        [Kappa(i)]=mnl_FindCohensKappa(ConfusionMatrices(i,1),ConfusionMatrices(i,3),ConfusionMatrices(i,2),ConfusionMatrices(i,4));
    else
        Kappa(i)=NaN;
    end
end
%% Calculate F1* Score
for i=1:nSets
    if ~isnan(ConfusionMatrices(i,:))
        Tp=ConfusionMatrices(i,1);
        Fn=ConfusionMatrices(i,2);
        Fp=ConfusionMatrices(i,3);
        Recall(i)=Tp/(Tp+Fn);
        Precision(i)=Tp/(Tp+Fp);
        F1score(i)=2*((Precision(i)*Recall(i))/(Precision(i)+Recall(i)));
    else
        Recall(i)=NaN;
        Precision(i)=NaN;
        F1score(i)=NaN;
    end
end
%View individual trace relationships between Precision and Recall
if strcmp(DisplayFig,'y')==1
    figure
    subplot(2,2,1)
    mnl_boxplot_NoChoices(Kappa',{'Score Per Cell'},'Kappa');
    ylim([0 1])
    xlim([0 2])
    ylabel('Kappa')
    title('Cohen Kappa')
    subplot(2,2,2)
    scatter(Recall,Precision,'.k')
    title('Recall vs Precision - Per Neuron')
    xlabel('Recall')
    ylabel('Precision')
    axis equal
    xlim([0 1])
    ylim([0 1])
    subplot(2,2,3)
    mnl_boxplot_NoChoices(F1score',{'Score Per Cell'},'F1* Score');
    title('F1 Scores Each Neuron')
    %ylim([0 1])
    xlim([0 2])
    ylabel('F1 Score')
end
%% Calculate the True Positive and False Positive Rates for ROC curve analysis
TruePositiveRate=Recall;
FalsePositiveRate=[];
for i=1:nSets
    Fp=ConfusionMatrices(i,3);
    Tn=ConfusionMatrices(i,4);
    FalsePositiveRate(i)=Fp/(Fp+Tn);
end
if strcmp(DisplayFig,'y')==1
    subplot(2,2,4)
    scatter(FalsePositiveRate,TruePositiveRate,'.k')
    axis equal
    xlim([0 1])
    ylim([0 1])
    title('ROC plot each neuron')
    ylabel('True Positive Rate')
    xlabel('False Positive Rate')
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
function [sVxList]=mnl_NestedScalePixels(VxList,ScaleFactor,NewDim)
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
end
function mnl_PlotSingleTrace2D(TracePoints,Diameters,Colour)
%Inputs
% TracePoints - The position of the points
% Diameters - The diameter of the trace
% Colour - the RGB values that you want to plot the trace in
%% Now make the volume for each trace
PointList=TracePoints(:,1:2);
szD=size(Diameters);
[nD,p]=max(szD);
%Work if the diameters are encoded in rows (p=1) or columns (p=2)
if p==1
    for j=1:nD
        Diameter(j,1)=Diameters(j,1);
    end
elseif p==2
    for j=1:nD
        Diameter(j,1)=Diameters(1,j);
    end
end
%Left Side
Xleft=PointList(:,1)-(Diameter/2);
%Right Side
Xright=PointList(:,1)+(Diameter/2);
X=[Xleft;flipud(Xright)];
%Sort Out Y
Y=[PointList(:,2);flipud(PointList(:,2))];
%Merge
patch(X,Y,Colour,'EdgeColor','none')
end