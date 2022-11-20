function [FinalClusters,FinalClusterIDs,NumTraceLim,MinDendriteLength,InitialClusterNum,Y]=mnl_EuclideanCrawler_Mag_Weighted_Axons(Trace,EuThresh,dim,FigureDisplay,Y,MinPoints,MinLength)
%Function to create clusters based on the Euclidean Distances weighted by
%the norm mean magnitude of each trace. "_Axons" has ammended the figure plotting
%a little by making the diameter of the cluster fragments a little wider
%and also using red instead of green. "_Axons" is recommended for thin axons, if
%using dendrites v3 is recommended
%final figure
%
% Inputs
%  Trace - the trace information
%  EuThresh - the euclidean threshold you want to use
%  dim - the dimesnions of the image (used to create the images)
%  FigureDisplay - Do you want to display the results of each cluster?(y/n)
%  Y - the tSNE display of the points
%  MinPoints - the minimum number of points required to be a cluster
%  MinLength - the mimumum cumulative length of the points to be a cluster
%
% Outputs
%  FinalClusters - Structure with the trace information sorted per cluster
%  FinalClusterIDs - The id of each trace in a single matrix
%  NumTraceLim - The minimum number of traces required to be considered a true cluster
%  MinDendriteLength - The minumum length of dendrite to be considered a true cluster
%  InitialClusterNum - The number of clusters created before filtering

%% Create the Input Matrix
nTrace=size(Trace,2);
nXFP=size(Trace(1).VecNormMean,2);
InputWeights=nan(nTrace,1);
InputMatrix=nan(nTrace,nXFP); % Done this way so if a value is missing we can see it
for i=1:nTrace
    if isnan(Trace(i).VecNormMean)
        InputValues=zeros(1,nXFP);
    else
        InputValues=Trace(i).VecNormMean;
    end
    InputMatrix(i,:)=InputValues;
    %Add the number of voxels as the weight
    InputWeights(i,1)=Trace(i).NormMeanMagnitude;
    clear InputValues
end
%% Run the Weighted Crawler
[ClusterIDs,Centroids]=mnl_WeightedEuclideanCrawler_v2(InputMatrix,InputWeights,EuThresh);
disp('Creating the Cluster Structure')
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
%% Filter 1 - At least n length and nPoints
disp('Filtering the Clusters...')
%NumTraceLim=10;
NumTraceLim=MinPoints;
%MinDendriteLength=1000;
MinDendriteLength=MinLength;
fC=0;
FinalClusterIDs=zeros(nTrace,1);
FinalClusters=struct('Traces',[],'Centroid',[]);
for i=1:ClusterNum
    NumTraces=size(Clusters(i).Traces,1);
    if NumTraces>=NumTraceLim
        Length=0;
        for j=1:NumTraces
            t=Clusters(i).Traces(j);
            %len=Trace(t).LengthVx;
            len=Trace(t).LengthReal;
            Length=Length+len;
        end
        if Length>=MinDendriteLength
            fC=fC+1;
            FinalClusters(fC).Traces=Clusters(i).Traces;
            FinalClusters(fC).Centroid=Clusters(i).Centroid;
            for j=1:NumTraces
                t=Clusters(i).Traces(j);
                FinalClusterIDs(t)=fC;
            end
        end
    end
end
ClusterNum=fC;
%% Images
if strcmp(FigureDisplay,'y')==1
    %% Plot a tSNE for the groups
    figure('Name','tSNE of filtered Traces')
    cmap=mnl_GenerateShuffledColourmap(ClusterNum);
    %If Y isn't specified
    if exist('Y','var')==0
        Iter=1000;
        Perp=10;
        Options=statset('MaxIter',Iter);
        Y=tsne(InputMatrix,'Options',Options,'Perplexity',Perp);
    elseif isempty(Y)==1
        Iter=1000;
        Perp=10;
        Options=statset('MaxIter',Iter);
        Y=tsne(InputMatrix,'Options',Options,'Perplexity',Perp);
    end
    PointMap=[0 0 0;cmap];
    gscatter(Y(:,1),Y(:,2),FinalClusterIDs,PointMap,'.',10)
    savefig('tSNE.fig')
    close all
    %% Plot a joint figure with all the traces together
    figure('Name','Plot of Traces - Together')
    Scale=Trace(1).Scale;
    %First draw all the traces in black
    NumTraces=size(Trace,2);
    for i=1:NumTraces
        Points=Trace(i).Points;
        Diameters=Trace(i).Diameter;
        mnl_Nested_PlotTraces2D(Points,Diameters,Scale,[0 0 0])
        hold on
    end    
    for i=1:ClusterNum
        NumTraces=size(FinalClusters(i).Traces,1);
        if NumTraces>=NumTraceLim
            for j=1:NumTraces
                t=FinalClusters(i).Traces(j);
                Points=Trace(t).Points;
                Diameters=Trace(t).Diameter;
                mnl_Nested_PlotTraces2D(Points,Diameters,Scale,cmap(i,:))
                hold on
            end
        end
    end
    xlim([0 dim(1)])
    ylim([0 dim(2)])
    axis ij
    savefig('AllTogether.fig')
    close all
    %% Now the map plots
    disp('Mapping Each Cluster')
    ChosenClusters=1:1:max(FinalClusterIDs);
    mnl_ImagesOfClusters_Axons(Trace,dim,FinalClusters,ChosenClusters,0.75,2);
end
end
%% Subfunctions
function mnl_Nested_PlotTraces2D(Points,Diameters,Scale,colour)
PointList=Points(:,1:2);
szD=size(Diameters);
[nD,p]=max(szD);
Diameter=zeros(nD,1);
if p==1
    for j=1:nD
        Diameter(j,1)=Diameters(j,1)./Scale(1);
    end
elseif p==2
    for j=1:nD
        Diameter(j,1)=Diameters(1,j)./Scale(1);
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
patch(X,Y,colour,'EdgeColor','none')
clear X Y PointList Diameter
end
