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
    for j=1:maxClusterID
        %Has the cluster been merged?
        if ~isnan(Centroids(j,:))
            Rows=find(ClusterIDs==j);
            Clusters(j).Traces=Rows;
            Clusters(j).Centroid=Centroids(j,:);
        end
    end
    ClusterNum=maxClusterID;
    InitialClusterNum=ClusterNum;
    %Evaluate Per Cell
    [~,F1score_dCrawler(:,n),Recall_dCrawler(:,n),Precision_dCrawler(:,n),TruePositiveRate_dCrawler(:,n),FalsePositiveRate_dCrawler(:,n)]=mnl_EvaluateClassifierPerCell(efPxTrace,ClusterIDs,Clusters,dim,'n',EuThreshVals(i),1);
    Setting_dCrawler(n,1)=EuThresh;
    n=n+1;
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
ClusteringMetrics(1).OptimumSettings_MeanIndex = I;
ClusteringMetrics(1).OptimumSettings_MedianIndex = maxI;
%% 
% Evaluate mean shift cluster
% %Mean Shift Cluster tuning
vec = [0.14252453, 0.19953435, 0.25654416, 0.31355398, 0.37056379, 0.4275736 ];
Bandwidths = vec;
szTh = size(Bandwidths,2);
n=1;

for i=1:szTh
    Bandwidth = Bandwidths(i);
    %Meanshift
    [ClusterIDs,Centroids] = mnl_MeanShift(InputMatrix,Bandwidth);
    %Create Cluster struct
    Clusters = struct('Traces',[],'Centroid',[]);
    maxClusterID= max(ClusterIDs);
    for j=1:maxClusterID
        Rows = find(ClusterIDs==j);
        Clusters(j).Traces=Rows;
        Clusters(j).Centroid=Centroids(j,:);
    end
    ClusterNum=maxClusterID;
    InitialClusterNum=ClusterNum;
    %Evaluate per cell
    [~,F1score_MeanShift(:,n),Recall_MeanShift(:,n),Precision_MeanShift(:,n),TruePositiveRate_MeanShift(:,n),FalsePositiveRate_MeanShift(:,n)]=mnl_EvaluateClassifierPerCell(efPxTrace,ClusterIDs,Clusters,dim,'n',Bandwidths(i),2);
    Setting_MeanShift(n,1)=Bandwidth;
    n=n+1;
end

%Now find the optimum thresh from the mean and median
F1scoreMean=mean(F1score_MeanShift,1,'omitnan');
F1scoreMedian=median(F1score_MeanShift,1,'omitnan');
F1scoreStd=std(F1score_MeanShift,1,'omitnan');

[M,I]=max(F1scoreMean);
[MedM,MedI]=max(F1scoreMedian);
maxI=MedI;
OptF1mean=M(1);
OptF1median=MedM(1);
OptThresh=Bandwidths(I(1)); %Defaults to the lowest value
OptThreshMedian=Bandwidths(MedI(1)); %Defaults to the lowest value

ClusteringMetrics(2).Name='Mean Shift Clustering';
ClusteringMetrics(2).OptF1mean=M(1);
ClusteringMetrics(2).OptF1median=MedM(1);
ClusteringMetrics(2).OptimumSettings_Mean=OptThresh; %Defaults to the lowest value
ClusteringMetrics(2).OptimumSettings_Median=OptThreshMedian;
ClusteringMetrics(2).OptimumSettings_MeanIndex = I;
ClusteringMetrics(2).OptimumSettings_MedianIndex = maxI;
%% K means 
% K means Cluster tuning
KMeansKVals=1:200;
szK=size(KMeansKVals,2);
n=1;
for i=1:szK
    K=KMeansKVals(i);
    %Clustering Algorithm
    [ClusterIDs,Centroids] = kmeans(InputMatrix, K, 'Distance', 'sqeuclidean');
    %Create the Clusters Structure
    Clusters=struct('Traces',[],'Centroid',[]);
    maxClusterID=max(ClusterIDs);
    for j=1:maxClusterID
        Rows=find(ClusterIDs==j);
        Clusters(j).Traces=Rows;
        Clusters(j).Centroid=Centroids(j,:);
    end
    ClusterNum=maxClusterID;
    InitialClusterNum=ClusterNum;
    %Evaluate Per Cell
    [~,F1score_KMeans(:,n),Recall_KMeans(:,n),Precision_KMeans(:,n),TruePositiveRate_KMeans(:,n),FalsePositiveRate_KMeans(:,n)]=mnl_EvaluateClassifierPerCell(efPxTrace,ClusterIDs,Clusters,dim,'n',KMeansKVals(i),1);
    Setting(n,1)=K;
    n=n+1;
end

%Now find the optimum thresh from the mean and median
F1scoreMean=mean(F1score_KMeans,1,'omitnan');
F1scoreMedian=median(F1score_KMeans,1,'omitnan');
F1scoreStd=std(F1score_KMeans,1,'omitnan');

[M,I]=max(F1scoreMean);
[MedM,MedI]=max(F1scoreMedian);
maxI=MedI;
OptF1mean=M(1);
OptF1median=MedM(1);
OptThresh=KMeansKVals(I(1)); %Defaults to the lowest value
OptThreshMedian=KMeansKVals(MedI(1)); %Defaults to the lowest value

ClusteringMetrics(3).Name='K Means Clustering';
ClusteringMetrics(3).OptF1mean=M(1);
ClusteringMetrics(3).OptF1median=MedM(1);
ClusteringMetrics(3).OptimumSettings_Mean=OptThresh; %Defaults to the lowest value
ClusteringMetrics(3).OptimumSettings_Median=OptThreshMedian; %Defaults to the lowest value
ClusteringMetrics(3).OptimumSettings_MeanIndex = I;
ClusteringMetrics(3).OptimumSettings_MedianIndex = maxI;

%% Evaluate DBSCAN
% DBScan tuning
EpsVals=0.05:0.05:5;
MinPtsVals=1;
szEps=size(EpsVals,2);
n=1;
for i=1:szEps
    Eps=EpsVals(i);
        %DBSCAN
        D = pdist2(InputMatrix,InputMatrix);
        DB = dbscan(D,Eps,MinPtsVals);
         ClusterIDs = DB+1; % DBSCAN assigns cluster IDs starting from 0, so we add 1 to match the previous implementation
        %Create the Clusters Structure
        Clusters=struct('Traces',[],'Centroid',[]);
        maxClusterID=max(ClusterIDs);
        for k=1:maxClusterID
            %Has the cluster been merged?
            Rows=find(ClusterIDs==k);
            Clusters(k).Traces=Rows;
            Clusters(k).Centroid=mean(InputMatrix(Rows,:),1);
        end
        ClusterNum=maxClusterID;
        InitialClusterNum=ClusterNum;
        %Evaluate Per Cell
        [~,F1score_DBSCAN(:,n),Recall_DBSCAN(:,n),Precision_DBSCAN(:,n),TruePositiveRate_DBSCAN(:,n),FalsePositiveRate_DBSCAN(:,n)]=mnl_EvaluateClassifierPerCell(efPxTrace,ClusterIDs,Clusters,dim,'n',EpsVals(i),1);
        Setting(n,:)=Eps;
        n=n+1;
end


%Now find the optimum thresh from the mean and median
F1scoreMean=mean(F1score_DBSCAN,1,'omitnan');
F1scoreMedian=median(F1score_DBSCAN,1,'omitnan');
F1scoreStd=std(F1score_DBSCAN,1,'omitnan');

[M,I]=max(F1scoreMean);
[MedM,MedI]=max(F1scoreMedian);
maxI=MedI;
OptF1mean=M(1);
OptF1median=MedM(1);
OptThresh=EpsVals(I(1)); %Defaults to the lowest value
OptThreshMedian=EpsVals(MedI(1)); %Defaults to the lowest value

ClusteringMetrics(4).Name='DBScan Clustering';
ClusteringMetrics(4).OptF1mean=M(1);
ClusteringMetrics(4).OptF1median=MedM(1);
ClusteringMetrics(4).OptimumSettings_Mean=OptThresh; %Defaults to the lowest value
ClusteringMetrics(4).OptimumSettings_Median=OptThreshMedian; %Defaults to the lowest value
ClusteringMetrics(4).OptimumSettings_MeanIndex = I;
ClusteringMetrics(4).OptimumSettings_MedianIndex = maxI;

%% ROC curves
figure('Name','ROC curves')
cmap_vec = [0 0 1;0 1 1;0 1 0;1 0.5 0;1 0 0];
cmap=cmap_vec;
%For dCrawler
cVal=1;
nRows=size(FalsePositiveRate_dCrawler,1);
mFPR=mean(FalsePositiveRate_dCrawler,1,'omitnan');
mTPR=mean(TruePositiveRate_dCrawler,1,'omitnan');
[smFPR,I]=sort(mFPR);%sort FPR
smTPR=mTPR(I);%sort TPR by FPR values
AUC=round(trapz(smFPR,smTPR),3);

Points_dCrawler=scatter(mFPR, mTPR, 'b', '.','MarkerFaceAlpha',0.25); %Plots the raw values
hold on
Line_dCrawler=plot(smFPR,smTPR,'Color',cmap(cVal,:),'LineWidth',2); %Plots a fitted line
str=sprintf('%s%s','AUC-dCrawler = ',num2str(AUC));
text(0.8,0.1*cVal,str)
OptPoint_Mean=scatter(mFPR(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),mTPR(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),10,cmap(cVal));
text(mFPR(ClusteringMetrics(cVal).OptimumSettings_MedianIndex)+0.05,mTPR(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),"<--- Optimum Threshold (via F1 score)")

ClusteringMetrics(cVal).AUC=AUC;

%For Mean Shift Clustering
cVal=2;
nRows=size(FalsePositiveRate_MeanShift,1);
mFPR=mean(FalsePositiveRate_MeanShift,1,'omitnan');
mTPR=mean(TruePositiveRate_MeanShift,1,'omitnan');
[smFPR,I]=sort(mFPR);%sort FPR
smTPR=mTPR(I);%sort TPR by FPR values
AUC=round(trapz(smFPR,smTPR),3);

Points_MeanShiftCluster=scatter(mFPR,mTPR,'r','.','MarkerFaceAlpha',0.25); %Plots the raw values
hold on
Line_MeanShiftCluster=plot(smFPR,smTPR,'Color',cmap(cVal,:),'LineWidth',2); %Plots a fitted line
str=sprintf('%s%s','AUC-MeanShiftCluster = ',num2str(AUC));
text(0.8,0.1*cVal,str)
OptPoint_Mean=scatter(mFPR(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),mTPR(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),10,cmap(cVal,:));
text(mFPR(ClusteringMetrics(cVal).OptimumSettings_MedianIndex)+0.05,mTPR(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),"<--- Optimum Threshold (via F1 score)")
ClusteringMetrics(cVal).AUC=AUC;

%For Kmeans Clustering
cVal=3;
nRows=size(FalsePositiveRate_KMeans,1);
mFPR=mean(FalsePositiveRate_KMeans,1,'omitnan');
mTPR=mean(TruePositiveRate_KMeans,1,'omitnan');
[smFPR,I]=sort(mFPR);%sort FPR
smTPR=mTPR(I);%sort TPR by FPR values
AUC=round(trapz(smFPR,smTPR),3);

Points_KMeans=scatter(mFPR,mTPR,'r','.','MarkerFaceAlpha',0.25); %Plots the raw values
hold on
Line_KMeans=plot(smFPR,smTPR,'Color',cmap(cVal,:),'LineWidth',2); %Plots a fitted line
str=sprintf('%s%s','AUC-KMeans = ',num2str(AUC));
text(0.8,0.1*cVal,str)
OptPoint_Mean=scatter(mFPR(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),mTPR(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),20,cmap(cVal,:));
text(mFPR(ClusteringMetrics(cVal).OptimumSettings_MedianIndex)+0.05,mTPR(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),"<--- Optimum Threshold (via F1 score)")
ClusteringMetrics(cVal).AUC=AUC;

% % %For K++means Clustering
% % cVal=4;
% % nRows=size(FalsePositiveRate_KPlus,1);
% % mFPR=mean(FalsePositiveRate_KPlus,1,'omitnan');
% % mTPR=mean(TruePositiveRate_KPlus,1,'omitnan');
% % [smFPR,I]=sort(mFPR);%sort FPR
% % smTPR=mTPR(I);%sort TPR by FPR values
% % AUC=round(trapz(smFPR,smTPR),3);
% 
% Points_KPlus=scatter(mFPR,mTPR,cmap(cVal,:),'.','MarkerFaceAlpha',0.25); %Plots the raw values
% hold on
% Line_KPlus=plot(smFPR,smTPR,'Color',cmap(cVal,:),'LineWidth',2); %Plots a fitted line
% str=sprintf('%s%s','AUC_KMeans = ',num2str(AUC));
% text(0.8,0.1*cVal,str)
% OptPoint_Mean=scatter(mFPR(maxI(1)),mTPR(maxI(1)),20,cmap(cVal));
% text(mFPR(maxI(1))+0.05,mTPR(maxI(1)),"<--- Optimum Threshold (via F1 score)")
% ClusteringMetrics(cVal).AUC=AUC;

%For DBScan
cVal=4;
nRows=size(FalsePositiveRate_DBSCAN,1);
mFPR=mean(FalsePositiveRate_DBSCAN,1,'omitnan');
mTPR=mean(TruePositiveRate_DBSCAN,1,'omitnan');
[smFPR,I]=sort(mFPR);%sort FPR
smTPR=mTPR(I);%sort TPR by FPR values
AUC=round(trapz(smFPR,smTPR),3);

Points_DBScan=scatter(mFPR,mTPR,'cyan','.','MarkerFaceAlpha',0.25); %Plots the raw values
hold on
Line_DBScan=plot(smFPR,smTPR,'Color',cmap(cVal,:),'LineWidth',2); %Plots a fitted line
str=sprintf('%s%s','AUC-DBScan = ',num2str(AUC));
text(0.8,0.1*cVal,str)
OptPoint_Mean=scatter(mFPR(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),mTPR(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),20,cmap(cVal));
text(mFPR(ClusteringMetrics(cVal).OptimumSettings_MedianIndex)+0.07,mTPR(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),"<--- Optimum Threshold (via F1 score)")
ClusteringMetrics(cVal).AUC=AUC;
% 
%Tidying up the graph
title('ROC curves')
xlim([0 1])
xlabel('False Positive Rate')
ylim([0 1])
ylabel('True Positive Rate')
%Add legend
legend([Line_dCrawler,Line_KMeans,Line_MeanShiftCluster,Line_DBScan],{'dCrawler','K Means','Mean-Shift','DBScan'})
%% Precision Recall curve
figure('Name','Precision Recall Curve')
%dCrawler
cVal=1;
% for i=1:nRows
%     plot(Recall_dCrawler(i,:),Precision_dCrawler(i,:),"Color",cmap(cVal,:),'LineWidth',0.01)
%     hold on
% end
mRecall=mean(Recall_dCrawler,1,'omitnan');
mPrecision=mean(Precision_dCrawler,1,'omitnan');
PR_Line_dCrawler=plot(mRecall,mPrecision,'Color',cmap(cVal,:),'LineWidth',2);
hold on
scatter(mRecall(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),mPrecision(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),20,cmap(cVal));
text(mRecall(ClusteringMetrics(cVal).OptimumSettings_MedianIndex)-0.7,mPrecision(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),"Optimum Threshold for dCrawler (via F1 score) --->")

%Mean Shift Clustering
cVal=2;
% for i=1:nRows
%     plot(Recall_MeanShift(i,:),Precision_MeanShift(i,:),"Color",cmap(cVal,:),'LineWidth',0.01)
%     hold on
% end
mRecall=mean(Recall_MeanShift,1,'omitnan');
mPrecision=mean(Precision_MeanShift,1,'omitnan');
PR_Line_MeanShiftCluster=plot(mRecall,mPrecision,'Color',cmap(cVal,:),'LineWidth',2);
scatter(mRecall(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),mPrecision(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),20,cmap(cVal));
text(mRecall(ClusteringMetrics(cVal).OptimumSettings_MedianIndex)-0.7,mPrecision(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),"Optimum Threshold for Mean Shift Clustering (via F1 score) --->")

%K Means
cVal=3;
% for i=1:nRows
%     plot(Recall_KMeans(i,:),Precision_KMeans(i,:),"Color",cmap(cVal,:),'LineWidth',0.01)
%     hold on
% end
mRecall=mean(Recall_KMeans,1,'omitnan');
mPrecision=mean(Precision_KMeans,1,'omitnan');
PR_Line_KMeans=plot(mRecall,mPrecision,'Color',cmap(cVal,:),'LineWidth',2);
scatter(mRecall(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),mPrecision(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),20,cmap(cVal));
text(mRecall(ClusteringMetrics(cVal).OptimumSettings_MedianIndex)-0.7,mPrecision(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),"Optimum Threshold for K Means Clustering (via F1 score) --->")

% %K Plus
% %cVal=4;
% %for i=1:nRows
% %    plot(Recall_KPlus(i,:),Precision_KPlus(i,:),"Color",cmap(cVal,:),'LineWidth',0.25)
%     hold on
% end
% mRecall=mean(Recall_KPlus,1,'omitnan');
% mPrecision=mean(Precision_KPlus,1,'omitnan');
% PR_Line_KPlus=plot(mRecall,mPrecision,'Color',cmap(cVal,:),'LineWidth',2);
% scatter(mRecall(maxI(1)),mPrecision(maxI(1)),20,cmap(cVal));
% text(mRecall(maxI(1))-0.7,mPrecision(maxI(1)),"Optimum Threshold for K Plus Clustering (via F1 score) --->")

%DBScan
cVal=4;
% for i=1:nRows
%     plot(Recall_DBSCAN(i,:),Precision_DBSCAN(i,:),"Color",cmap(cVal,:),'LineWidth',0.01)
%     hold on
% end
mRecall=mean(Recall_DBSCAN,1,'omitnan');
mPrecision=mean(Precision_DBSCAN,1,'omitnan');
PR_Line_DBScan=plot(mRecall,mPrecision,'Color',cmap(cVal,:),'LineWidth',2);
scatter(mRecall(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),mPrecision(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),20,cmap(cVal));
text(mRecall(ClusteringMetrics(cVal).OptimumSettings_MedianIndex)-0.7,mPrecision(ClusteringMetrics(cVal).OptimumSettings_MedianIndex),"Optimum Threshold forDBScan (via F1 score) --->")

%Details for graph
title('Precision Recall curves')
xlim([0 1])
xlabel('Recall')
ylim([0 1])
ylabel('Precision')
legend([PR_Line_dCrawler,PR_Line_KMeans,PR_Line_MeanShiftCluster,PR_Line_DBScan],{'dCrawler','K means','Mean-Shift','DBScan'})
