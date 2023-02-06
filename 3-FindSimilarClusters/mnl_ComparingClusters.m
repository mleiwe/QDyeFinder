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
        InputValues=Trace(i).VecNormMean;
    end
    InputMatrix(i,:)=InputValues;
    %Add the number of voxels as the weight
    InputWeights(i,1)=Trace(i).NormMeanMagnitude;
    clear InputValues
end

%% Evaluate the dCrawler
EuThreshVals=0.05:0.01:1;
szTh=size(EuThreshVals,2);
n=1;
for i=1:szTh
    EuThresh=EuThreshVals(i);
    [ClusterIDs,Centroids]=mnl_WeightedEuclideanCrawler_v2(InputMatrix,InputWeights,EuThresh);
    %Evaluate Per Cell
    [~,F1score_dCrawler(:,n),Recall_dCrawler(:,n),Precision_dCrawler(:,n),TruePositiveRate_dCrawler(:,n),FalsePositiveRate_dCrawler(:,n)]=mnl_EvaluateClassifierPerCell(efPxTrace,ClusterIDs,Clusters,dim,'n',EuThreshVals(i),1);
    Setting_dCrawler(n,1)=EuThresh;
end

%% Evaluate mean shift cluster

%% Evaluate k-means

%% Evaluate k++

%% Evaluate DBSCAN


%% Precision Recall Curves