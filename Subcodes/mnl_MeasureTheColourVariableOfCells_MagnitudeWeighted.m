function [efPxTrace]=mnl_MeasureTheColourVariableOfCells_MagnitudeWeighted(efPxTrace,CellFigs)
%Functions to measure the accuracy of traces, the true grouping of the
%traces/fragments is mentioned in the SetId of the structure. And the
%colour distance is calculated from the weighted colour mean of the Set (NB weighted to the magnitude of each fragment).
%
% Inputs
%  efPxTrace - Structure with the colours etc...
%  CellFigs - Whether you want to see the distributions of each cell
% Outputs
%  efPxTrace - updated to now include the colour distance of each trace to
%              the weighted mean colour of their set
NumTraces=size(efPxTrace,2);
dim=efPxTrace(1).dim;
Scale=efPxTrace(1).Scale;
%% Create a Cell Id List and other basic info
BkgRemovedMatrix=nan(NumTraces,dim(3));
RawMeanMatrix=nan(NumTraces,dim(3));
n=0;
for i=1:NumTraces
    setname=efPxTrace(i).SetId;
    if isempty(setname)==0
        n=n+1;
        CellID_List_NoEmptys{n}=setname;
    end
    CellID_List{i}=efPxTrace(i).SetId;
    %Basic Matrices
    BkgRemovedMatrix(i,:)=efPxTrace(i).BkgRemovedMean;
    RawMeanMatrix(i,:)=efPxTrace(i).RawMean;
end
VecNormScale=[0 1];
NormMeanScale=[0 1];
BkgRemovedScale=[0 max(BkgRemovedMatrix(:))];
RawMeanScale=[0 max(RawMeanMatrix(:))];
%% Associate the Traces to each cell
CellIDs=unique(CellID_List_NoEmptys);
NumCells=size(CellIDs,2);
for i=1:NumCells
    test=CellIDs{i};
    tCellIDs(i)=str2double(test(6:end));
end
[~,B]=sort(tCellIDs);
for i=1:NumCells
    fCellIDs{i}=CellIDs{B(i)};
end
for i=1:NumCells
    String1=fCellIDs{i};
    idx=find(strcmp(String1,CellID_List)==1);
    Cells(i).Id=String1;
    Cells(i).Traces=idx;
    CellIDs{i}=String1;
end
%% Now Create the centroid and measure the distance of each trace to the centre
Max_nT=1;
RawCentroids=nan(dim(3),NumCells);
BkgRemovedCentroids=nan(dim(3),NumCells);
NormMeanCentroids=nan(dim(3),NumCells);
Centroids=nan(dim(3),NumCells);
wCentroids=nan(dim(3),NumCells);
for i=1:NumCells
    if strcmp(CellFigs,'y')==1
        fign=sprintf('%s%s','Basic Info - ',Cells(i).Id);
        fign2=sprintf('%s%s',Cells(i).Id,'-BasicInfo');
        CellFig=figure('Name',fign,'WindowState','fullscreen');
    end
    Traces=Cells(i).Traces;
    nT=size(Traces,2);
    Lengths=[];
    nVx=[];
    
    MagnitudeList=[];
    VecNormMatrix=nan(nT,dim(3));
    NormMeanMatrix=nan(nT,dim(3));
    BkgRemovedMatrix=nan(nT,dim(3));
    RawMeanMatrix=nan(nT,dim(3));
    for j=1:nT
        %Get Colours
        VecNormMatrix(j,:)=efPxTrace(Traces(j)).VecNormMean;
        NormMeanMatrix(j,:)=efPxTrace(Traces(j)).NormMean;
        BkgRemovedMatrix(j,:)=efPxTrace(Traces(j)).BkgRemovedMean;
        RawMeanMatrix(j,:)=efPxTrace(Traces(j)).RawMean;
        %Get the number of voxels
        nVx(j,1)=size(efPxTrace(Traces(j)).AllVoxels,1);
        %Get Magnitude
        MagnitudeList(j,1)=efPxTrace(Traces(j)).NormMeanMagnitude;
    end
    %Calculate the weighted mean and each colour distance per trace
    %Now Calculations
    TotalMagnitude=sum(MagnitudeList);
    %Get centroids (arithmetic mean) and key values
    Cells(i).Centroid=mean(VecNormMatrix,1,'omitnan');
    Cells(i).RawCentroid=mean(RawMeanMatrix,1,'omitnan');
    Cells(i).BkgRemovedCentroid=mean(BkgRemovedMatrix,1,'omitnan');
    Cells(i).NormMeanCentroid=mean(NormMeanMatrix,1,'omitnan');
    Cells(i).NormMeanMagnitude=MagnitudeList;
    Cells(i).NumberVoxels=nVx;

    %Get Weighted Centroid 
    Weights=MagnitudeList/TotalMagnitude;
    szVN=size(VecNormMatrix);
    wVecNormMatrix=nan(szVN);
    for j=1:szVN(1)
        wVecNormMatrix(j,:)=VecNormMatrix(j,:)*Weights(j);
    end
    swVecNormMatrix=sum(wVecNormMatrix,1);
    Cells(i).WeightedCentroid=swVecNormMatrix;
    RawCentroids(:,i)=Cells(i).RawCentroid';
    BkgRemovedCentroids(:,i)=Cells(i).BkgRemovedCentroid';
    NormMeanCentroids(:,i)=Cells(i).NormMeanCentroid';
    Centroids(:,i)=Cells(i).Centroid';
    wCentroids(:,i)=Cells(i).WeightedCentroid';        
    %Measure Distance to centroid
    [Dist]=mnl_GroupEuD_v2(Cells(i).Centroid,VecNormMatrix);
    Cells(i).Distances=Dist;
    %Measure Distance to Weighted Centroid
    [Dist]=mnl_GroupEuD_v2(Cells(i).WeightedCentroid,VecNormMatrix);
    Cells(i).WeightedDistances=Dist;
    for j=1:nT %Add the distances to efPxTrace
        tNum=Traces(j);
        efPxTrace(tNum).DistToWeightedCentroid_Mag=Dist(j);
    end
    %Update the Max_nT
    if Max_nT<nT
        Max_nT=nT;
    end  
    %Now if individual figures are specified
    if strcmp(CellFigs,'y')==1
        %create a color map with 51 points for 0 to 0.5 euclidean distances
        cmap=flipud(magma(50));
        cmap(51,:)=[0 0 0];
        DistanceRange=linspace(0,0.5,51);
        %Plot the individual fragments with colour traces
        subplot(1,3,1)
        for j=1:nT
            %Determine the colour to use
            Cdist=Dist(j);
            [~,cNum]=min(abs(DistanceRange-Cdist));
            mnl_PlotSingleTrace(efPxTrace(Traces(j)),dim,Scale,cmap(cNum,:))
            hold on
            %Now write the trace number near the trace
            str=sprintf('%s%d','#',j);
            mnl_AddTextToFig(efPxTrace(Traces(j)).Points,str,dim)
        end
        hTrace=gca;
        colormap(hTrace,cmap)
        c=colorbar('Ticks',0:0.2:1,'TickLabels',{'0','0.1','0.2','0.3','0.4','0.5'});
        c.Label.String = 'Colour Distance in EuD';
        subplot(2,6,3)
        imagesc(RawMeanMatrix)
        hRawMean=gca;
        colormap(hRawMean,magma)
        colorbar
        title('Raw Mean')
        subplot(2,6,4)
        imagesc(BkgRemovedMatrix)
        hBkgMean=gca;
        colormap(hBkgMean,magma)
        colorbar
        title('Background Substracted')
        subplot(2,6,9)
        imagesc(NormMeanMatrix,NormMeanScale)
        hNormMean=gca;
        colormap(hNormMean,magma)
        colorbar
        title('Normalised Mean')
        subplot(2,6,10)
        imagesc(VecNormMatrix,VecNormScale)
        hVecNormMean=gca;
        colormap(hVecNormMean,magma)
        colorbar
        title('Vector Normalised')
        %Plot distances as a boxplot
        subplot(1,3,3)
        mnl_boxplot2(Dist,'','Distance to centroid','y','y')

        savefig(CellFig,fign);
        mnl_ExportEPSdense(CellFig,fign2);
        close(CellFig)
    end
end

h=figure('Name','CentroidColours');
colormap(magma)
subplot(5,1,1)
imagesc(RawCentroids);
colorbar
title('Raw Mean Cell Centroid')
subplot(5,1,2)
imagesc(BkgRemovedCentroids);
colorbar
title('Bkg Removed Mean Cell Centroid')
subplot(5,1,3)
imagesc(NormMeanCentroids,[0 1]);
colorbar
title('Norm Mean Cell Centroid')
subplot(5,1,4)
imagesc(Centroids,[0 1]);
colorbar
title('Vec Norm Mean Cell Centroid')
subplot(5,1,5)
imagesc(wCentroids,[0 1]);
colorbar
title('Weighted VecNorm Mean Centroid')
savefig(h,'MeanForEachCell')
close all
%% Create a matrix for boxplots - Using weighted mean
BoxPlotMatrix=nan(Max_nT,NumCells);
MagMatrix=nan(Max_nT,NumCells);
nVxMatrix=nan(Max_nT,NumCells);
AllDist=[];
wAllDist=[];
wBoxPlotMatrix=nan(Max_nT,NumCells);
for i=1:NumCells
    %Get the valuess
    wDist=Cells(i).WeightedDistances;
    Mags=Cells(i).NormMeanMagnitude;
    nVoxels=Cells(i).NumberVoxels;
    szD=size(wDist,1);
    %Assign to the box plot matrices
    wBoxPlotMatrix(1:szD,i)=wDist;
    MagMatrix(1:szD,i)=Mags;
    nVxMatrix(1:szD,i)=nVoxels;
    wAllDist=[wAllDist;wDist];
end
%Now Weighted
wMeanPerCell=mean(wBoxPlotMatrix,1,'omitnan');
wOverallMeanDist=mean(wMeanPerCell,'omitnan');
wOverallStdList=std(wMeanPerCell,'omitnan');
wAverageOfEachTrace=mean(wAllDist,'omitnan');
wStdOfEachTrace=std(wAllDist,'omitnan');
%Distance to weighted mean figure
h=figure('Name','Per Cell');
fn=sprintf('%s%s%s','DistancePerCell');
mnl_boxplot2(wBoxPlotMatrix,CellIDs,'Distance to  Weighted Cell Centroid','y','y');
grid on
ylim([0 1.4])
savefig(h,fn,'compact');
%Magnitude of fragments sorted by cell
h=figure('Name','Magnitude Per Cell');
fn=sprintf('%s','MagnitudeOfEachTrace');
subplot(2,1,1)
mnl_boxplot2(MagMatrix,CellIDs,'Magnitude of Trace','y','y');
grid on
title('Magnitude of Traces')
subplot(2,1,2) %Scatter to see the relationship between the magnitude and colour distance
idx=~isnan(MagMatrix); %the non-nan values
AllMagVals(:,1)=MagMatrix(idx);
AllColourDistances(:,1)=wBoxPlotMatrix(idx);
scatter(AllMagVals,AllColourDistances,'.')
ylabel('Colour distance to weighted cell mean')
xlabel('Magnitude of trace')
savefig(h,fn,'compact');
%Number of voxels per trace to colour distance (scatter plot)
h=figure('Name','Number of Voxels Per Cell');
fn=sprintf('%s','NumberOfVoxelsofEachTrace');
subplot(2,1,1)
mnl_boxplot2(nVxMatrix,CellIDs,'nVx Per Trace','y','y');
grid on
title('Number of voxels per Trace')
subplot(2,1,2)
AllnVx(:,1)=nVxMatrix(idx);
scatter(AllnVx,AllColourDistances,'.')
ylabel('Colour distance to weighted cell mean')
xlabel('nVx per trace')
savefig(h,fn,'compact');

h=figure('Name','All Traces');
mnl_boxplot2(wAllDist,'Each Trace','Average Distance to Cell Centroid','y','y');
grid on
ylim([0 1.4])
title('Weighted Mean')
savefig(h,'DistancePerTrace','compact');
%% Now plot the EuDistance to the Magnitude of the trace
%Each cell as a subplot
nCol=ceil(sqrt(NumCells));
nRow=round(sqrt(NumCells));
%Mag vs EuD
AllX=[];AllY=[];
h=figure('Name','Magnitude vs Euclidean Distance');
for i=1:NumCells
    subplot(nRow,nCol,i)
    Y=Cells(i).WeightedDistances;
    X=Cells(i).NormMeanMagnitude;
    scatter(X,Y,'.')
    xlabel('Magnitude')
    ylabel('Eu Distance')
    xlim([0 2])
    ylim([0 1.4])
    title(CellIDs{i})
    AllX=[AllX;X];AllY=[AllY;Y];
    clear X Y
end
savefig(h,'Magnitude_VsEuDist_PerCell')
%nVx vs EuD
AllX=[];AllY=[];
h=figure('Name','Number of Voxels vs Euclidean Distance');
for i=1:NumCells
    subplot(nRow,nCol,i)
    Y=Cells(i).WeightedDistances;
    X=Cells(i).NumberVoxels;
    scatter(X,Y,'.')
    xlabel('Number of Voxels')
    ylabel('Eu Distance')
    ylim([0 1.4])
    title(CellIDs{i})
    AllX=[AllX;X];AllY=[AllY;Y];
    clear X Y
end
savefig(h,'NumVoxels_VsEuDist_PerCell')
close all
end

%%Subfunctions
function [DistMatrix]=mnl_GroupEuD_v2(Pos,Matrix)
%% New Method
%Calculate the distances
DistMatrix=sqrt(sum((Matrix-Pos).^2,2));
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
function mnl_PlotSingleTrace(efPxTrace,dim,Scale,ChosenColour)
nTraces=size(efPxTrace,2);
%% Now make the volume for each soma
for i=1:nTraces
    PointList=efPxTrace(i).Points(:,1:2);
    Diameters=efPxTrace(i).Diameter;
    szD=size(Diameters);
    [nD,p]=max(szD);
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
    patch(X,Y,ChosenColour,'EdgeColor','none')
    clear X Y PointList Diameter
end
axis equal
axis ij
xlim([0 dim(1)])
ylim([0 dim(2)])
end
function mnl_AddTextToFig(Points,str,dim)
MeanXY=mean(Points,1,'omitnan');
MeanXY(1)=MeanXY(1)+20;
if MeanXY(1)>=dim(2)
    MeanXY(1)=MeanXY(1)-20;
end
if MeanXY(2)>=dim(1)
    MeanXY(2)=MeanXY(2)-20;
end
text(MeanXY(1),MeanXY(2),str,'Color','k');
end