function [Thresh]=mnl_EvaluateColourDistribution(nDim,efPxTrace)
%% Option 1 - Evaluate the range of colours in the traces - cumulative plot
% Possion distribution at various copy numbers
Spreads=[0.5 1 2 4 6 8 10];
NumPoints=size(efPxTrace,2); %Do the same number as the number of somas
szSp=size(Spreads,2);
% Generate Poission Distributions For Each Copy Number
[ModelCells]=mnl_GeneratePossionsNchannels(NumPoints,Spreads,nDim);
for i=1:szSp
    ModelDist(:,i)=ModelCells(i).NormXFPvals(:,1);
    SpreadId{i}=sprintf('%s%s',num2str(ModelCells(i).CopyNumber),' Copies');
end
% Get the colours
szT=size(efPxTrace,2);
Trace_VecNormMatrix=nan(szT,nDim);
for i=1:szT
    Trace_VecNormMatrix(i,:)=efPxTrace(i).VecNormMean;
end
for i=1:nDim
    ChanID{i}=sprintf('%s%d','Ch ',i);
end
%Now the figure
figure('Name','Colour Distributions of the Traces compared to the model')
%Subplot 1 - The model only
subplot(1,3,1)
mnl_CumulativePlotMatrix(ModelDist)
legend(SpreadId)
xlabel('Vector Normalised Value')
title('Model Cumulative Distribution')
%Subplot 2 - The traces only
subplot(1,3,2)
mnl_CumulativePlotMatrix(Trace_VecNormMatrix)
legend(ChanID)
xlabel('Vector Normalised Value')
title('Trace Cumulative Distribution')
%Subplot 3 - Two distinct boxplots
subplot(2,3,3)
mnl_boxplot_NoChoices(ModelDist,SpreadId,'Vec Norm Value');
title('Per Copy Number')
subplot(2,3,6)
mnl_boxplot_NoChoices(Trace_VecNormMatrix,ChanID,'Vec Norm Value');
title('Per Channel')
%Another figure as a boxplot on the same plot space
AllVecNorm=[];
for i=1:nDim
    AllVecNorm=[AllVecNorm;Trace_VecNormMatrix(:,i)];
end
figure('Name','All Together')
sz1=size(AllVecNorm,1);
jsz=max([sz1 NumPoints]);
JointBoxplot=nan(jsz,szSp+1);
for i=1:szSp
    JointBoxplot(1:NumPoints,i)=ModelDist(:,i);
end
JointBoxplot(1:sz1,szSp+1)=AllVecNorm;
JointLegend=SpreadId;
JointLegend{szSp+1}='Model Data';
mnl_boxplot_NoChoices(JointBoxplot,JointLegend,'Vec Norm Value');
ylim([0 1])

%% Option 2 - Evaluate the discriminability of the traces labelled (if they exist)
EuDist=0:0.01:0.5;
szEu=size(EuDist,2);
%For model
for i=1:szSp
    Matrix=ModelCells(i).NormXFPvals;
    [Model_Discrim_Pc(i,:),Model_Discrim_PcStd(i,:),~]=mnl_EvaluateThePermutation(Matrix,EuDist);        
    clear Matrix
end
%For Traces
[Trace_Discrim_Pc,Trace_Discrim_PcStd,~]=mnl_EvaluateThePermutation(Trace_VecNormMatrix,EuDist);
%Now figure
figure('Name','Discriminability of Traces')
subplot(1,2,1)
cmap=colormap(jet(szSp+1));
for i=1:szSp
    yMean=Model_Discrim_Pc(i,:);
    yStd=Model_Discrim_PcStd(i,:);
    P2=patch([EuDist fliplr(EuDist)], [yMean+yStd fliplr(yMean-yStd)], cmap(i,:),'EdgeColor','none');
    P2.FaceAlpha=0.2;
    hold on
    p(i)=plot(EuDist,yMean,'LineWidth',2,'Color',cmap(i,:));
    legname{i}=sprintf('%s%s',num2str(ModelCells(i).CopyNumber),' Copies');
    clear yMean yStd
end
%Now plot traces
yMean=Trace_Discrim_Pc;
yStd=Trace_Discrim_PcStd;
P2=patch([EuDist fliplr(EuDist)], [yMean+yStd fliplr(yMean-yStd)], cmap(szSp+1,:),'EdgeColor','none');
P2.FaceAlpha=0.2;
p(szSp+1)=plot(EuDist,yMean,'LineWidth',3,'Color',cmap(szSp+1,:));
legname{szSp+1}='Traces';
ylabel('Percent Discriminable')
xlabel('Euclidean Distance Threshold')
legend(p,legname)
ylim([0 100])
xlim([0 0.5])
subplot(1,2,2)
for i=1:szSp
    yMean=Model_Discrim_Pc(i,:);
    yStd=Model_Discrim_PcStd(i,:);
    P2=patch([EuDist fliplr(EuDist)], [yMean+yStd fliplr(yMean-yStd)], cmap(i,:),'EdgeColor','none');
    P2.FaceAlpha=0.2;
    hold on
    p(i)=plot(EuDist,yMean,'LineWidth',2,'Color',cmap(i,:));
    legname{i}=sprintf('%s%s',num2str(round(ModelCells(i).CopyNumber,1)),' Copies');
    clear yMean yStd
end
%Now plot traces
yMean=Trace_Discrim_Pc;
yStd=Trace_Discrim_PcStd;
P2=patch([EuDist fliplr(EuDist)], [yMean+yStd fliplr(yMean-yStd)], cmap(szSp+1,:),'EdgeColor','none');
P2.FaceAlpha=0.2;
p(szSp+1)=plot(EuDist,yMean,'LineWidth',2,'Color',cmap(szSp+1,:));
legname{szSp+1}='Traces';
ylim([90 100])
xlim([0 0.5])
ylabel('Percent Discriminable')
xlabel('Euclidean Distance Threshold')
%% Find where the traces cut the chosen percentile
PcChosen=95;
idx_Trace=find(Trace_Discrim_Pc>=PcChosen,1,'last');
Thresh=EuDist(idx_Trace);
%Now for each model too
Model_Thresh=nan(szSp,2);
idx_Model=nan(1,szSp);
for i=1:szSp
    idx_Model(i)=find(Model_Discrim_Pc(i,:)>=PcChosen,1,'last');
    Model_Thresh(i,1)=Spreads(i);
    if isempty(idx_Model(i))==0
        Model_Thresh(i,2)=EuDist(idx_Model(i));
    else
        Model_Thresh(i,2)=NaN;
    end
end
%Is there sufficient spread?
if Thresh<0.2
    fprintf('%s\n','Warning! There is limited colour spread in these traces')
end
%Which copy number is it closest too?
Diffs=Model_Thresh(:,2)-Thresh;
Dists=sqrt(Diffs.^2);
[~,I]=min(Dists);
ClosestCopy=sprintf('%s%s',num2str(Model_Thresh(I,1)),' copies');
fprintf('%s%s\n','The closest copy number at the 95% discriminability level is... ',ClosestCopy);
end
function [Discrim_Pc,Discrim_PcStd,PcUnique]=mnl_EvaluateThePermutation(SomaVals,EuThresh)
%Now measure the distances
%[EuD_All,Eu_Matrix]=mnl_GroupColourEuclidean_ListAndMatrix(SomaVals);
[EuD_Matrix]=mnl_GroupEuclidean_Matrixv2(SomaVals,SomaVals);

%% For Each EuThresh
nThresh=size(EuThresh,2);
nSomas=size(SomaVals,1);
PcUnique=nan(1,nThresh);
Discrim_Pc=nan(1,nThresh);
Discrim_Pcstd=nan(1,nThresh);
for i=1:nThresh
    Thresh=EuThresh(i);
    idx=EuD_Matrix<=Thresh;
    nsame=sum(idx,1);
    %Calculate the Percent Unique
    UniqueIdx=nsame==1;
    UniqueCounter=sum(UniqueIdx);
    PcUnique(i)=(UniqueCounter/nSomas)*100;
    %Calculate the discriminability
    nPairs=nsame-1;
    PcPairs=(nPairs./(nSomas-1)).*100;
    Discrim_Percent=100-PcPairs;
    Discrim_Pc(i)=mean(Discrim_Percent);
    Discrim_PcStd(i)=std(Discrim_Percent);
    
%     %legacy code much slower
%     %Calculate the percent unique
%     UniqueCounter=0;
%     for j=1:nSomas
%         EuVals=EuD_Matrix(j,:);
%         [~,loc]=find(EuVals<=Thresh);
%         nPairs=find(loc~=j);
%         szNP=size(nPairs,2);
%         if isempty(nPairs)==1
%             UniqueCounter=UniqueCounter+1;
%             Discrim_Percent(j)=100;
%         else
%             Discrim_Percent(j)=((nSomas-1-szNP)/(nSomas-1))*100;
%         end
%     end
%     Discrim_Pc(i)=mean(Discrim_Percent);
%     Discrim_PcStd(i)=std(Discrim_Percent);
%     PcUnique(i)=(UniqueCounter/nSomas)*100;
%     clear Discrim_Percent
end
end