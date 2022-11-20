function [ClusterIDs,Centroids]=mnl_WeightedEuclideanCrawler_v2(InputMatrix,InputWeights,EuThresh)
%% Step 1 - Initial Crawl
[~,nCol]=size(InputMatrix);
%disp('Performing the initial crawl...')
[Centroids,ClusterIDs]=mnl_EuclideanCrawl_Weighted(InputMatrix,InputWeights,EuThresh);
ClusterNum=size(Centroids,1);
%fprintf('%d%s',ClusterNum,' clusters in the initial crawl')
%% Step 2 - Adjusting clusters
%disp('Adjusting Clusters...')
[Centroids,ClusterIDs]=mnl_EuclideanAdjust(InputMatrix,InputWeights,ClusterIDs,Centroids);
ClusterNum=size(Centroids,1);
%% Step 3 - Merge Clusters if they are too close
ClusterThresh=EuThresh;
%disp('Merging Nearby Clusters')
mxClusterNum=ClusterNum;
%Measure the distances between clusters
mnLoop=1;
while 1
    nMerged=0;
    i=1;
    while i<=mxClusterNum
        Pos1=Centroids(i,:);
        tCent=Centroids;
        tCent(i,:)=nan(1,nCol);
        [DistMatrix]=mnl_GroupEuD_v2(Pos1,tCent);
        [loc]=find(DistMatrix<ClusterThresh);
        if isempty(loc)~=1
            nMerged=nMerged+1;
            %Now merge the closest one
            [~,MinLoc]=min(DistMatrix(loc));
            PotentialMerge=[i;loc(MinLoc)];
            NewClusterNumber=min(PotentialMerge);
            %Calculate the new centroid position
            tCentroid=nan(2,nCol);
            tWeight=nan(2,1);
            for j=1:2
                tCentroid(j,:)=Centroids(PotentialMerge(j),:);
                %Get the appropriate weights
                Weight=0;
                PointsIdx=ClusterIDs==PotentialMerge(j);
                nP=sum(PointsIdx);
                for k=1:nP
                    Weight=Weight+InputWeights(PotentialMerge(j));
                end
                tWeight(j,1)=Weight;
            end
            TempCentroidPos=mnl_WeightedMean(tCentroid,tWeight);
            clear tCentroid tWeight
            %Now recode the clusters
            for j=1:2
                if NewClusterNumber==PotentialMerge(j)
                    %Is the point close enough to the centroid?
                    [RowList]=find(ClusterIDs==PotentialMerge(j));
                    tMatrix=InputMatrix(RowList,:);
                    [DistMatrix]=mnl_GroupEuD_v2(TempCentroidPos,tMatrix);
                    %Allocate those within with the new number
                    [locInc]=find(DistMatrix<=EuThresh);
                    [locExc]=find(DistMatrix>EuThresh);
                    if isempty(locInc)==0
                        ClusterIDs(RowList(locInc))=NewClusterNumber;
                    end
                    %Score those outside as zero
                    if isempty(locExc)==0
                        ClusterIDs(RowList(locExc))=0;
                    end
                    clear tMatrix
                elseif NewClusterNumber~=PotentialMerge(j)
                    [RowList]=find(ClusterIDs==PotentialMerge(j));
                    tMatrix=InputMatrix(RowList,:);
                    [DistMatrix]=mnl_GroupEuD_v2(TempCentroidPos,tMatrix);
                    %Allocate those within with the new number
                    [locInc]=find(DistMatrix<=EuThresh);
                    [locExc]=find(DistMatrix>EuThresh);
                    if isempty(locInc)==0
                        ClusterIDs(RowList(locInc))=NewClusterNumber;
                    end
                    %Score those outside as zero
                    if isempty(locExc)==0
                        ClusterIDs(RowList(locExc))=0;
                    end
                    clear tMatrix
                    %Score the discarded centroids as nan
                    Centroids(PotentialMerge(j),:)=nan(1,nCol);
                end
            end
            clear TempCentroidPos
            %Update the centroids
            %Re-calculate the new centroid position
            [RowList]=find(ClusterIDs==NewClusterNumber);
            tPoints=InputMatrix(RowList,:);
            tWeight=InputWeights(RowList);
            Centroids(NewClusterNumber,:)=mnl_WeightedMean(tPoints,tWeight);
            clear tPoints tWeight          
        end  
        i=i+1;
    end
    %Are There any Points that are unlabelled?
    [RowList]=find(ClusterIDs==0);
    %If yes then we need to crawl through again
    if isempty(RowList)==0
        [Centroids,ClusterIDs]=mnl_CrawlEmptyOnly(InputMatrix,InputWeights,RowList,ClusterIDs,Centroids,EuThresh);
        %Adjust the new clusters
        %fprintf('\n%s\n','Readjusting Clusters After Merging...')
        [Centroids,ClusterIDs]=mnl_EuclideanAdjust(InputMatrix,InputWeights,ClusterIDs,Centroids);
    end
    %Now check if another loop is needed
    if nMerged>=1
        %txt=sprintf('%s%d%s%d%s','Merge Loop ',mnLoop,' Merged ',nMerged,' Cluster Sets');
        %fprintf('%s\n',txt)
        mnLoop=mnLoop+1;
        mxClusterNum=max(ClusterIDs);
        i=1; %Now Start from the begining
        %if it has looped more than 20 times reduce the merging threshold by 1%
        if mnLoop>=21
            ClusterThresh=ClusterThresh*0.95;
        end
    elseif nMerged==0
        break
    end
end
end
%% Sub functions
%Weighted Mean Function
function [wMean]=mnl_WeightedMean(Points,Weights)
%Remove nans
PointsIdx=~isnan(Points(:,1));
WeightsIdx=~isnan(Weights(:,1));
cPoints=Points(PointsIdx,:);
cWeights=Weights(WeightsIdx,:);

nPoints=size(cPoints,1);
nWeights=size(cWeights,1);
if nPoints~=nWeights
    error('The number of points and the number of weights do not match')
end
tPoints=cPoints.*cWeights;
sWeights=sum(cWeights);
sPoints=sum(tPoints,1);
wMean=sPoints./sWeights;
end
%Measure the Euclidean Distance
function [DistMatrix]=mnl_GroupEuD_v2(Pos,Matrix)
% New Method
%Calculate the distances
DistMatrix=sqrt(sum((Matrix-Pos).^2,2));
end
%Initial Crawl Function
function [Centroids,ClusterIDs]=mnl_EuclideanCrawl_Weighted(InputMatrix,InputWeights,EuThresh)
i=1;
ClusterNum=1;
nCol=size(InputMatrix,2);
nRow=size(InputMatrix,1);
IntermediateMatrix=InputMatrix;
IntermediateWeights=InputMatrix;
Points=nan(nRow,nCol);
Weights=nan(nRow,1);
n_CP=1;%Cluster Point Counter
ClusterIDs=zeros(nRow,1);
% Asign the first values
NewPoints=InputMatrix(i,:);
NewWeights=InputWeights(i);
IntermediateMatrix(i,:)=nan(1,nCol);
ClusterIDs(i)=ClusterNum;
Points(n_CP,:)=NewPoints;
Weights(n_CP)=NewWeights;
Centroids(n_CP,:)=mnl_WeightedMean(Points,Weights);
%Centroids(i,:)=Points(1:n_CP,:);
n_CP=n_CP+1;
TracesPerCluster=1;
TracesOfCluster(TracesPerCluster)=i;
%% Step 1 - The first Crawl through
while i<=nRow
    %Calculate the centroid
    %Centroid=mean(Points,1,'omitnan');
    Centroid=mnl_WeightedMean(Points,Weights);
    Centroids(ClusterNum,:)=Centroid;
    % Check if the points are still close enough to the centroid
    [tDistMatrix]=mnl_GroupEuD_v2(Centroid,Points); %This is faster
    [Pos]=find(tDistMatrix>EuThresh);
    if isempty(Pos)==0
        NumberOutside=size(Pos,2);
        %Return Values back for assessment
        for j=1:NumberOutside
            %Get the Row Id
            RowID=TracesOfCluster(Pos(j));
            %Return the colour to the comparing matrix
            IntermediateMatrix(RowID,:)=InputMatrix(RowID,:);
            IntermediateWeights(RowID)=InputWeights(RowID);
            %Return the cluster ids
            ClusterIDs(RowID)=0;
            %Now Mark the points for deletion
            Points(Pos(j),:)=nan(1,nCol); %initially mark them NaN
            Weights(Pos(j))=NaN;
            TracesOfCluster(Pos(j))=NaN;
        end
        %Remove NaNs and reset the counters
        idx_nan=~isnan(Points(:,1));
        n=sum(idx_nan);%Number of Points left
        tPoints=nan(nRow,nCol);
        tPoints(1:n,:)=Points(idx_nan,:);
        tWeights=nan(nRow,1);
        tWeights(1:n)=Weights(idx_nan,:);
        n_CP=n+1;

        idx_nanClust=~isnan(TracesOfCluster);
        TracesOfCluster=TracesOfCluster(idx_nanClust);
        Points=tPoints;
        Weights=tWeights;
        TracesPerCluster=TracesPerCluster-NumberOutside;
        %Now Recalculate the centroid
        Centroid=mnl_WeightedMean(Points,Weights);
        clear Pos tPoints tWeights
    end    
    %Calculate the distance to the Center
    [DistMatrix]=mnl_GroupEuD_v2(Centroid,IntermediateMatrix); %Faster
    %Find the minumum distance
    [val,~]=min(DistMatrix);
    %Is it below the Euclidean Distance
    if val<=EuThresh
        %Speed up the calculations by also adding in the points at the same distance
        minPoints=DistMatrix==val;
        n_minPoints=sum(minPoints);
        %Update the locations
        ClusterIDs(minPoints)=ClusterNum;
        NewPoints=IntermediateMatrix(minPoints,:);
        NewWeights=IntermediateWeights(minPoints);
        Points(n_CP:n_CP+n_minPoints-1,:)=NewPoints;
        Weights(n_CP:n_CP+n_minPoints-1)=NewWeights;
        n_CP=n_CP+n_minPoints;
        IntermediateMatrix(minPoints,:)=nan(n_minPoints,nCol);       
        IntermediateWeights(minPoints)=nan(n_minPoints,1);
        loc=find(minPoints==1);
        TracesOfCluster(TracesPerCluster+1:TracesPerCluster+n_minPoints)=loc;
        TracesPerCluster=TracesPerCluster+n_minPoints;
    else
        %Save as a cluster
        Centroids(ClusterNum,:)=Centroid;
        clear TracesOfCluster NewPoints NewWeights
        %Make A new Cluster
        ClusterNum=ClusterNum+1;
        %Move to the next value
        i=1;
        while 1            
            if i<=nRow
                if ~isnan(IntermediateMatrix(i,:)) %keep cycling through until we find a value that hasn't been allocated to a cluster
                    [a,~]=find(ClusterIDs~=0);
                    numDefined=size(a,1);
                    %mnl_InsertProgressTrackerInLoops(numDefined,nRow)                   
                    break
                else
                    i=i+1;
                end
            else
                break
            end
        end
        %Redefine if i<=nRow
        if i<=nRow
            n_CP=1;
            Points=nan(nRow,nCol);
            Weights=nan(nRow,1);
            NewPoints=IntermediateMatrix(i,:);
            NewWeights=IntermediateWeights(i);
            IntermediateMatrix(i,:)=nan(1,nCol);
            IntermediateWeights(i)=NaN;
            ClusterIDs(i)=ClusterNum;
            Points(n_CP,:)=NewPoints;
            Weights(n_CP)=NewWeights;
            n_CP=n_CP+1;
            TracesPerCluster=1;
            TracesOfCluster(1)=i;
        else
            break
        end
    end    
end
end
%Adjust Crawl Function
function [Centroids,ClusterIDs]=mnl_EuclideanAdjust(InputMatrix,InputWeights,ClusterIDs,Centroids)
i=1;
nLoop=0;
nRow=size(InputMatrix,1);
while i<=nRow %NB Cannot just loop it need to keep checking at each stage until all are associated with their nearest centroid
    [DistMatrix]=pdist2(InputMatrix,Centroids); %Measure the distances to each centroid
    [~,loc]=min(DistMatrix,[],2);%Find the closest cluster
    idx=ClusterIDs~=loc;%find the ones that belong to a different cluster than their closest one
    nAdjust=sum(idx);%How many need to be adjusted
    if nAdjust>0
        %Re-assign the Cluster IDs
        ClusterIDs(idx)=loc(idx);
        %Re-calculate the new centroid positions
        tn=size(Centroids,1);
        for j=1:tn
            idx_loc=ClusterIDs==j;%All the points belonging to that cluster
            MeanVal=mnl_WeightedMean(InputMatrix(idx_loc,:),InputWeights(idx_loc));
            Centroids(j,:)=MeanVal;
        end
        %Display loop information
        nLoop=nLoop+1;
        %str=sprintf('%s%d%s%d%s','Done ',nLoop,' Loops with ',nAdjust,' Adjustments');
        %fprintf('%s\n',str)
        i=1;
    elseif nAdjust==0
        i=nRow+1;
    end
end
end
%Crawl through the empty ones only
function [Centroids,ClusterIDs]=mnl_CrawlEmptyOnly(InputMatrix,InputWeights,RowList,ClusterIDs,Centroids,EuThresh)
%RowList is the Unassigned points
%% Key Numbers
nRow=size(InputMatrix,1);
nCol=size(InputMatrix,2);
%% Find the new cluster number
ClustNumList=unique(ClusterIDs,'rows');
BiggestNumber=max(ClusterIDs);
c=1;
for j=1:BiggestNumber
    [loc]=find(ClustNumList==j);
    if isempty(loc)==1
        MissingClustNum(c)=j;
        c=c+1;
    end
end
nMissingClusters=c-1;
%% Create a New Matrix for unassigned points
TempPoints=InputMatrix(RowList,:);
TempWeights=InputWeights(RowList);
%Make new clusters
if nMissingClusters>0
    ClusterNum=MissingClustNum(1);
else
    ClusterNum=max(ClustNumList)+1;
end
%% Pre-allocate and first values
i=1;
IntermediateMatrix=nan(nRow,nCol);
IntermediateWeights=nan(nRow,1);
Points=nan(nRow,nCol);
NewPoints=TempPoints(i,:);
Weights=nan(nRow,1);
NewWeights=TempWeights(i);
n_CP=1;%Cluster Point Counter
% Asign the first values
IntermediateMatrix(RowList,:)=TempPoints;
IntermediateWeights(RowList,1)=TempWeights;
ClusterIDs(RowList(1))=ClusterNum;
Points(n_CP,:)=NewPoints;
Weights(n_CP)=NewWeights;
Centroids(ClusterNum,:)=mnl_WeightedMean(Points(1:n_CP,:),Weights(1:n_CP));
%Centroids(ClusterNum,:)=Points(1:n_CP,:);
n_CP=n_CP+1;
TracesPerCluster=1;
TracesOfCluster(TracesPerCluster)=RowList(i);
%% Now crawl
while i<=nRow
    Centroid=mnl_WeightedMean(Points,Weights);
    %Centroid=mean(Points,1,'omitnan');
    Centroids(ClusterNum,:)=Centroid;
    %% Check if the points are still close enough to the centroid
    [tDistMatrix]=mnl_GroupEuD_v2(Centroid,Points); %This is faster
    [Pos]=find(tDistMatrix>EuThresh);
    if isempty(Pos)==0
        NumberOutside=size(Pos,2);
        %Return Values back for assessment
        for j=1:NumberOutside
            %Get the Row Id
            RowID=TracesOfCluster(Pos(j));
            %Return the colour to the comparing matrix
            IntermediateMatrix(RowID,:)=InputMatrix(RowID,:);
            IntermediateWeights(RowID,:)=InputWeights(RowID,:);
            %Return the cluster ids
            ClusterIDs(RowID)=0;
            %Now Mark the points for deletion
            Points(Pos(j),:)=nan(1,nCol); %initially mark them NaN
            TracesOfCluster(Pos(j))=NaN;
        end
        %Remove NaNs and reset the counters
        idx_nan=~isnan(Points(:,1));
        n=sum(idx_nan);%Number of Points left
        tPoints=nan(nRow,nCol);
        tWeights=nan(nRow,1);
        tPoints(1:n,:)=Points(idx_nan,:);
        tWeights(1:n,:)=Weights(idx_nan);
        n_CP=n+1;
        idx_nanClust=~isnan(TracesOfCluster);
        TracesOfCluster=TracesOfCluster(idx_nanClust);
        Points=tPoints;
        Weights=tWeights;
        TracesPerCluster=TracesPerCluster-NumberOutside;
        %Now Recalculate the centroid
        Centroid=mnl_WeightedMean(Points,Weights);
        %Centroid=mean(Points,1,'omitnan');
        clear Pos tPoints tWeights
    end    
    %% Find the next point
    %Calculate the distance to the Center
    [DistMatrix]=mnl_GroupEuD_v2(Centroid,IntermediateMatrix); %Faster
    %Find the minumum distance
    [val,~]=min(DistMatrix);
    if val<=EuThresh
        %Speed up the calculations by also adding in the points at the same distance
        minPoints=DistMatrix==val;
        n_minPoints=sum(minPoints);
        %Update the locations
        ClusterIDs(minPoints)=ClusterNum;
        NewPoints=IntermediateMatrix(minPoints,:);
        NewWeights=IntermediateWeights(minPoints,1);
        Points(n_CP:n_CP+n_minPoints-1,:)=NewPoints;
        Weights(n_CP:n_CP+n_minPoints-1,1)=NewWeights;
        n_CP=n_CP+n_minPoints;
        IntermediateMatrix(minPoints,:)=nan(n_minPoints,nCol);       
        IntermediateWeights(minPoints)=nan(n_minPoints,1);
        loc=find(minPoints==1);
        TracesOfCluster(TracesPerCluster+1:TracesPerCluster+n_minPoints)=loc;
        TracesPerCluster=TracesPerCluster+n_minPoints;
    else
        %% Cluster Fully identified
        %Save as a cluster
        Centroids(ClusterNum,:)=Centroid;
        clear TracesOfCluster NewPoints
        %Find the new cluster number
        ClustNumList=unique(ClusterIDs,'rows');
        BiggestNumber=max(ClusterIDs);
        c=1;
        MissingClustNum=[];
        for k=1:BiggestNumber
            [loc]=find(ClustNumList==k);
            if isempty(loc)==1
                MissingClustNum(c)=k;
                c=c+1;
            end
        end
        nMissingClusters=c-1;
        if c>1
            ClusterNum=min(MissingClustNum);
        else
            ClusterNum=BiggestNumber+1;
        end
        TracesPerCluster=1;
        %Move to the next value
        i=1;
        while 1
            if i<=nRow
                if isnan(IntermediateMatrix(i,:))==0 %keep cycling through until we find a value that hasn't been allocated to a cluster
                    [a,~]=find(ClusterIDs~=0);
                    numDefined=size(a,1);
                    %mnl_InsertProgressTrackerInLoops(numDefined,nRow)
                    break
                else
                    i=i+1;
                end
            else
                break
            end
        end
        %Redefine if i<=nRow
        if i<=nRow
            n_CP=1;
            Points=nan(nRow,nCol);
            Weights=nan(nRow,1);
            NewPoints=IntermediateMatrix(i,:);
            NewWeights=IntermediateWeights(i,1);
            IntermediateMatrix(i,:)=nan(1,nCol);
            IntermediateWeights(i)=NaN;
            ClusterIDs(i)=ClusterNum;
            Points(n_CP,:)=NewPoints;
            Weights(n_CP)=NewWeights;
            n_CP=n_CP+1;
            TracesPerCluster=1;
            TracesOfCluster(1)=i;
        else
            break
        end
    end
end
end