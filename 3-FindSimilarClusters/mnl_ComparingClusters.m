function mnl_ComparingClusters(efPxTrace)

%% Extract the imput matrix and weights (magnitude)
dim=efPxTrace(1).dim;
nTrace=size(efPxTrace,2);
nXFP=size(efPxTrace(1).VecNormMean,2);
InputWeights=nan(nTrace,1);
InputMatrix=nan(nTrace,nXFP); % Done this way so if a value is missing we can see it
for i=1:nTrace
    if isnan(efPxTrace(i).VecNormMean)
        InputValues=zeros(1,nXFP);
    else
        InputValues=efPxTrace(i).VecNormMean;
    end
    InputMatrix(i,:)=InputValues;
    %Add the number of voxels as the weight
    InputWeights(i,1)=efPxTrace(i).NormMeanMagnitude;
    clear InputValues
end

%% Evaluate the dCrawler
EuThreshVals=0.05:0.01:1;
szTh=size(EuThreshVals,2);
n=1;
for i=1:szTh
    EuThresh=EuThreshVals(i);
    %Clustering Algorithm
    [ClusterIDs,Centroids]=mnl_WeightedEuclideanCrawler_v2(InputMatrix,InputWeights,EuThresh);
    %Create the Clusters Structure
    Clusters=struct('Traces',[],'Centroid',[]);
    maxClusterID=max(ClusterIDs);
    for i=1:maxClusterID
        %Has the cluster been merged?
        if ~isnan(Centroids(i,:))
            Rows=find(ClusterIDs==i);
            Clusters(i).Traces=Rows;
            Clusters(i).Centroid=Centroids(i,:);
        end
    end
    ClusterNum=maxClusterID;
    InitialClusterNum=ClusterNum;
    %Evaluate Per Cell
    [~,F1score_dCrawler(:,n),Recall_dCrawler(:,n),Precision_dCrawler(:,n),TruePositiveRate_dCrawler(:,n),FalsePositiveRate_dCrawler(:,n)]=mnl_EvaluateClassifierPerCell(efPxTrace,ClusterIDs,Clusters,dim,'n',EuThreshVals(i),1);
    Setting_dCrawler(n,1)=EuThresh;
end
%Now find the optimum thresh from the mean and median
F1scoreMean=mean(F1score_dCrawler,1,'omitnan');
F1scoreMedian=median(F1score_dCrawler,1,'omitnan');
F1scoreStd=std(F1score_dCrawler,1,'omitnan');
[M,I]=max(F1scoreMean);
[MedM,MedI]=max(F1scoreMedian);
maxI=MedI;
OptF1mean=M(1);
OptF1median=MedM(1);
OptThresh=EuThreshVals(I(1)); %Defaults to the lowest value
OptThreshMedian=EuThreshVals(MedI(1)); %Defaults to the lowest value

ClusteringMetrics(1).Name='dCrawler';
ClusteringMetrics(1).OptF1mean=M(1);
ClusteringMetrics(1).OptF1median=MedM(1);
ClusteringMetrics(1).OptimumSettings_Mean=OptThresh; %Defaults to the lowest value
ClusteringMetrics(1).OptimumSettings_Median=OptThreshMedian; %Defaults to the lowest value
%% Evaluate mean shift cluster
%Mean Shift Cluster tuning
MeanShiftSettings=[];


%Now find the optimum thresh from the mean and median
F1scoreMean=mean(F1score_MeanShiftCluster,1,'omitnan');
F1scoreMedian=median(F1score_MeanShiftCluster,1,'omitnan');
F1scoreStd=std(F1score_MeanShiftCluster,1,'omitnan');

[M,I]=max(F1scoreMean);
[MedM,MedI]=max(F1scoreMedian);
maxI=MedI;
OptF1mean=M(1);
OptF1median=MedM(1);
OptThresh=MeanShiftSettings(I(1),:); %Defaults to the lowest value
OptThreshMedian=MeanShiftSettings(MedI(1),:); %Defaults to the lowest value

ClusteringMetrics(2).Name='Mean Shift Clustering';
ClusteringMetrics(2).OptF1mean=M(1);
ClusteringMetrics(2).OptF1median=MedM(1);
ClusteringMetrics(2).OptimumSettings_Mean=OptThresh; %Defaults to the lowest value
ClusteringMetrics(2).OptimumSettings_Median=OptThreshMedian;
%% K means 
% K means Cluster tuning
KMeansSettings=[];


%Now find the optimum thresh from the mean and median
F1scoreMean=mean(F1score_KMeans,1,'omitnan');
F1scoreMedian=median(F1score_KMeans,1,'omitnan');
F1scoreStd=std(F1score_KMeans,1,'omitnan');

[M,I]=max(F1scoreMean);
[MedM,MedI]=max(F1scoreMedian);
maxI=MedI;
OptF1mean=M(1);
OptF1median=MedM(1);
OptThresh=KMeansSettings(I(1),:); %Defaults to the lowest value
OptThreshMedian=KMeansSettings(MedI(1),:); %Defaults to the lowest value

ClusteringMetrics(3).Name='K Means Clustering';
ClusteringMetrics(3).OptF1mean=M(1);
ClusteringMetrics(3).OptF1median=MedM(1);
ClusteringMetrics(3).OptimumSettings_Mean=KMeansSettings(I(1)); %Defaults to the lowest value
ClusteringMetrics(3).OptimumSettings_Median=KMeansSettings(MedI(1)); %Defaults to the lowest value
%% Evaluate k++
% K plus plus Cluster tuning
KPlusSettings=[];


%Now find the optimum thresh from the mean and median
F1scoreMean=mean(F1score_KPlus,1,'omitnan');
F1scoreMedian=median(F1score_KPlus,1,'omitnan');
F1scoreStd=std(F1score_KPlus,1,'omitnan');

[M,I]=max(F1scoreMean);
[MedM,MedI]=max(F1scoreMedian);
maxI=MedI;
OptF1mean=M(1);
OptF1median=MedM(1);
OptThresh=KMeansSettings(I(1),:); %Defaults to the lowest value
OptThreshMedian=KMeansSettings(MedI(1),:); %Defaults to the lowest value

ClusteringMetrics(4).Name='K Plus Plus Means Clustering';
ClusteringMetrics(4).OptF1mean=M(1);
ClusteringMetrics(4).OptF1median=MedM(1);
ClusteringMetrics(4).OptimumSettings_Mean=KPlusSettings(I(1)); %Defaults to the lowest value
ClusteringMetrics(4).OptimumSettings_Median=KPlusSettings(MedI(1)); %Defaults to the lowest value
%% Evaluate DBSCAN
% DBScan tuning
DBScanSettings=[];


%Now find the optimum thresh from the mean and median
F1scoreMean=mean(F1score_DBScan,1,'omitnan');
F1scoreMedian=median(F1score_DBScan,1,'omitnan');
F1scoreStd=std(F1score_DBScan,1,'omitnan');

[M,I]=max(F1scoreMean);
[MedM,MedI]=max(F1scoreMedian);
maxI=MedI;
OptF1mean=M(1);
OptF1median=MedM(1);
OptThresh=DBScanSettings(I(1),:); %Defaults to the lowest value
OptThreshMedian=DBScanSettings(MedI(1),:); %Defaults to the lowest value

ClusteringMetrics(4).Name='DBScan Clustering';
ClusteringMetrics(4).OptF1mean=M(1);
ClusteringMetrics(4).OptF1median=MedM(1);
ClusteringMetrics(4).OptimumSettings_Mean=DBScanSettings(I(1)); %Defaults to the lowest value
ClusteringMetrics(4).OptimumSettings_Median=DBScanSettings(MedI(1)); %Defaults to the lowest value

%% ROC curves
cmap=colormap(jet(5));
figure('Name','ROC curves')
%For dCrawler
cVal=1;
nRows=size(FalsePositiveRate_dCrawler,1);
mFPR=mean(FalsePositiveRate_dCrawler,1,'omitnan');
mTPR=mean(TruePositiveRate_dCrawler,1,'omitnan');
[smFPR,I]=sort(mFPR);%sort FPR
smTPR=mTPR(I);%sort TPR by FPR values
AUC=round(trapz(smFPR,smTPR),3);

Points_dCrawler=scatter(mFPR,mTPR,cmap(cVal,:),'.','MarkerFaceAlpha',0.25); %Plots the raw values
hold on
Line_dCrawler=plot(smFPR,smTPR,'Color',cmap(cVal,:),'LineWidth',2); %Plots a fitted line
str=sprintf('%s%s','AUC_dCrawler = ',num2str(AUC));
text(0.8,0.1*cVal,str)
OptPoint_Mean=scatter(mFPR(maxI(1)),mTPR(maxI(1)),20,cmap(cVal));
text(mFPR(maxI(1))+0.05,mTPR(maxI(1)),"<--- Optimum Threshold (via F1 score)")

ClusteringMetrics(cVal).AUC=AUC;

%For Mean Shift Clustering
cVal=2;
nRows=size(FalsePositiveRate_MeanShiftCluster,1);
mFPR=mean(FalsePositiveRate_MeanShiftCluster,1,'omitnan');
mTPR=mean(TruePositiveRate_MeanShiftCluster,1,'omitnan');
[smFPR,I]=sort(mFPR);%sort FPR
smTPR=mTPR(I);%sort TPR by FPR values
AUC=round(trapz(smFPR,smTPR),3);

Points_MeanShiftCluster=scatter(mFPR,mTPR,cmap(cVal,:),'.','MarkerFaceAlpha',0.25); %Plots the raw values
hold on
Line_MeanShiftCluster=plot(smFPR,smTPR,'Color',cmap(cVal,:),'LineWidth',2); %Plots a fitted line
str=sprintf('%s%s','AUC_MeanShiftCluster = ',num2str(AUC));
text(0.8,0.1*cVal,str)
OptPoint_Mean=scatter(mFPR(maxI(1)),mTPR(maxI(1)),20,cmap(cVal));
text(mFPR(maxI(1))+0.05,mTPR(maxI(1)),"<--- Optimum Threshold (via F1 score)")
ClusteringMetrics(cVal).AUC=AUC;

%For Kmeans Clustering
cVal=3;
nRows=size(FalsePositiveRate_KMeans,1);
mFPR=mean(FalsePositiveRate_KMeans,1,'omitnan');
mTPR=mean(TruePositiveRate_KMeans,1,'omitnan');
[smFPR,I]=sort(mFPR);%sort FPR
smTPR=mTPR(I);%sort TPR by FPR values
AUC=round(trapz(smFPR,smTPR),3);

Points_KMeans=scatter(mFPR,mTPR,cmap(cVal,:),'.','MarkerFaceAlpha',0.25); %Plots the raw values
hold on
Line_KMeans=plot(smFPR,smTPR,'Color',cmap(cVal,:),'LineWidth',2); %Plots a fitted line
str=sprintf('%s%s','AUC_KMeans = ',num2str(AUC));
text(0.8,0.1*cVal,str)
OptPoint_Mean=scatter(mFPR(maxI(1)),mTPR(maxI(1)),20,cmap(cVal));
text(mFPR(maxI(1))+0.05,mTPR(maxI(1)),"<--- Optimum Threshold (via F1 score)")
ClusteringMetrics(cVal).AUC=AUC;

%For K++means Clustering
cVal=4;
nRows=size(FalsePositiveRate_KPlus,1);
mFPR=mean(FalsePositiveRate_KPlus,1,'omitnan');
mTPR=mean(TruePositiveRate_KPlus,1,'omitnan');
[smFPR,I]=sort(mFPR);%sort FPR
smTPR=mTPR(I);%sort TPR by FPR values
AUC=round(trapz(smFPR,smTPR),3);

Points_KPlus=scatter(mFPR,mTPR,cmap(cVal,:),'.','MarkerFaceAlpha',0.25); %Plots the raw values
hold on
Line_KPlus=plot(smFPR,smTPR,'Color',cmap(cVal,:),'LineWidth',2); %Plots a fitted line
str=sprintf('%s%s','AUC_KMeans = ',num2str(AUC));
text(0.8,0.1*cVal,str)
OptPoint_Mean=scatter(mFPR(maxI(1)),mTPR(maxI(1)),20,cmap(cVal));
text(mFPR(maxI(1))+0.05,mTPR(maxI(1)),"<--- Optimum Threshold (via F1 score)")
ClusteringMetrics(cVal).AUC=AUC;

%For DBScan
cVal=5;
nRows=size(FalsePositiveRate_DBScan,1);
mFPR=mean(FalsePositiveRate_DBScan,1,'omitnan');
mTPR=mean(TruePositiveRate_DBScan,1,'omitnan');
[smFPR,I]=sort(mFPR);%sort FPR
smTPR=mTPR(I);%sort TPR by FPR values
AUC=round(trapz(smFPR,smTPR),3);

Points_DBScan=scatter(mFPR,mTPR,cmap(cVal,:),'.','MarkerFaceAlpha',0.25); %Plots the raw values
hold on
Line_DBScan=plot(smFPR,smTPR,'Color',cmap(cVal,:),'LineWidth',2); %Plots a fitted line
str=sprintf('%s%s','AUC_KMeans = ',num2str(AUC));
text(0.8,0.1*cVal,str)
OptPoint_Mean=scatter(mFPR(maxI(1)),mTPR(maxI(1)),20,cmap(cVal));
text(mFPR(maxI(1))+0.05,mTPR(maxI(1)),"<--- Optimum Threshold (via F1 score)")
ClusteringMetrics(cVal).AUC=AUC;

%Tidying up the graph
xlim([0 1])
xlabel('False Positive Rate')
ylim([0 1])
ylabel('True Positive Rate')
%Add legend
legend([Line_dCrawler,Line_MeanShiftCluster,Line_KMeans,Line_KPlus,Line_DBScan],{'dCrawler','Mean Shift Clustering','K Means', 'K++', 'DBScan'})
%% Precision Recall curve
figure('Name','Precision Recall Curve')
%dCrawler
cVal=1;
for i=1:nRows
    plot(Recall_dCrawler(i,:),Precision_dCrawler(i,:),"Color",cmap(cVal,:),'LineWidth',0.25)
    hold on
end
mRecall=mean(Recall_dCrawler,1,'omitnan');
mPrecision=mean(Precision_dCrawler,1,'omitnan');
PR_Line_dCrawler=plot(mRecall,mPrecision,'Color',cmap(cVal,:),'LineWidth',2);
scatter(mRecall(maxI(1)),mPrecision(maxI(1)),20,cmap(cVal));
text(mRecall(maxI(1))-0.7,mPrecision(maxI(1)),"Optimum Threshold for dCrawler (via F1 score) --->")

%Mean Shift Clustering
cVal=2;
for i=1:nRows
    plot(Recall_MeanShiftCluster(i,:),Precision_MeanShiftCluster(i,:),"Color",cmap(cVal,:),'LineWidth',0.25)
    hold on
end
mRecall=mean(Recall_MeanShiftCluster,1,'omitnan');
mPrecision=mean(Precision_MeanShiftCluster,1,'omitnan');
PR_Line_MeanShiftCluster=plot(mRecall,mPrecision,'Color',cmap(cVal,:),'LineWidth',2);
scatter(mRecall(maxI(1)),mPrecision(maxI(1)),20,cmap(cVal));
text(mRecall(maxI(1))-0.7,mPrecision(maxI(1)),"Optimum Threshold for Mean Shift Clustering (via F1 score) --->")

%K Means
cVal=3;
for i=1:nRows
    plot(Recall_KMeans(i,:),Precision_KMeans(i,:),"Color",cmap(cVal,:),'LineWidth',0.25)
    hold on
end
mRecall=mean(Recall_KMeans,1,'omitnan');
mPrecision=mean(Precision_KMeans,1,'omitnan');
PR_Line_KMeans=plot(mRecall,mPrecision,'Color',cmap(cVal,:),'LineWidth',2);
scatter(mRecall(maxI(1)),mPrecision(maxI(1)),20,cmap(cVal));
text(mRecall(maxI(1))-0.7,mPrecision(maxI(1)),"Optimum Threshold for K Means Clustering (via F1 score) --->")

%K Plus
cVal=4;
for i=1:nRows
    plot(Recall_KPlus(i,:),Precision_KPlus(i,:),"Color",cmap(cVal,:),'LineWidth',0.25)
    hold on
end
mRecall=mean(Recall_KPlus,1,'omitnan');
mPrecision=mean(Precision_KPlus,1,'omitnan');
PR_Line_KPlus=plot(mRecall,mPrecision,'Color',cmap(cVal,:),'LineWidth',2);
scatter(mRecall(maxI(1)),mPrecision(maxI(1)),20,cmap(cVal));
text(mRecall(maxI(1))-0.7,mPrecision(maxI(1)),"Optimum Threshold for K Plus Clustering (via F1 score) --->")

%DBScan
cVal=5;
for i=1:nRows
    plot(Recall_DBScan(i,:),Precision_DBScan(i,:),"Color",cmap(cVal,:),'LineWidth',0.25)
    hold on
end
mRecall=mean(Recall_DBScan,1,'omitnan');
mPrecision=mean(Precision_DBScan,1,'omitnan');
PR_Line_DBScan=plot(mRecall,mPrecision,'Color',cmap(cVal,:),'LineWidth',2);
scatter(mRecall(maxI(1)),mPrecision(maxI(1)),20,cmap(cVal));
text(mRecall(maxI(1))-0.7,mPrecision(maxI(1)),"Optimum Threshold forDBScan (via F1 score) --->")

%Details for graph
title('Precision Recall curves')
xlim([0 1])
xlabel('Recall')
ylim([0 1])
ylabel('Precision')
legend([PR_Line_dCrawler,PR_Line,MeanShiftCluster,PR_Line_MeanShiftCluster,PR_Line_Kmeans,PR_Line_KPlus,PR_Line_DBScan],{'dCrawler','Mean Shift Clustering','K means','K Plus','DBScan'})