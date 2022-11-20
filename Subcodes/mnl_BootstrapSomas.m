function [Dimension]=mnl_BootstrapSomas(Somas)
nTrials=200;
nCells=[10 25 50 100 150 200 250 300];
sz1=size(nCells,2);
EuThresh=0:0.025:0.5;
szTh=size(EuThresh,2);
nSomas=size(Somas,2);
maxDim=size(Somas(1).VecNormMean,2); %The maximum number of dimensions
Chans=1:1:maxDim;
for i=1:maxDim
    [PermCombo]=mnl_PermsNoColumnBias(Chans,i);
    Perms(i).PermCombo=PermCombo;
end
%% Create the matrix to compare
SomaVals=nan(nSomas,maxDim);
for i=1:nSomas
    SomaVals(i,:)=Somas(i).VecNormMean;
end
%% Evaluate each permutation
j=1;
for i=2:maxDim
    if i==3
        ck=1;
    end
    nComb=size(Perms(i).PermCombo,1);
    CombNum=randperm(nComb,1); %Randomly pick one combination
    Comb=Perms(i).PermCombo(CombNum,:);
    Dimension(i).Comb(j).Comb=Comb;
    % Sub Select Depending on the Permutations
    tSomaVals=SomaVals(:,Comb);
    for k=1:sz1
        %How many cells to choose
        nCell=nCells(k);
        UniqueStructure(k).NumberOfCellsChosen=nCell;
        TempMatrix=nan(nTrials,szTh);
        for m=1:nTrials
            %Shuffle the Cells
            ShuffledMatrix=[];
            Chosen_Cells=randperm(nSomas,nCell); %Which ones are chosen
            ShuffledMatrix=tSomaVals(Chosen_Cells,:);
            %For each threshold
            [PcUnique]=mnl_CalculateThePercentageUnique(ShuffledMatrix,EuThresh);
            for n=1:szTh
                UniqueStructure(k).Trial(m).Thresh(n).Threshold=EuThresh(n);
                TempMatrix(m,n)=PcUnique(n);
                UniqueStructure(k).Trial(m).Thresh(n).PercentUnique=PcUnique(n);
            end
        end
        %Now Work out the mean and standard deviation
        MeanList=nanmean(TempMatrix,1);
        StdList=nanstd(TempMatrix,1);
        for m=1:szTh
            UniqueStructure(k).PerThresh(m).Mean=MeanList(m);
            UniqueStructure(k).PerThresh(m).Std=StdList(m);
        end
        Dimension(i).Comb(j).NumberOfCells(k)=UniqueStructure(k);
    end
end
%% Plot a figure for each number of dimensions
for i=2:maxDim
    fign=sprintf('%d%s',i,' dimensions');
    h=figure('Name',fign);
    %Make Sure We have the x axis labels
    if exist('EuThresh','var')==0
        szTh=size(Dimension(i).Comb.NumberOfCells(1).Trial(1).Thresh,2);
        for j=1:szTh
            EuThresh(j)=Dimension(i).Comb.NumberOfCells(1).Trial(1).Thresh(j).Threshold;
        end
    end
    %How many cell groups are there?
    if exist('nCells','var')==0
        sz1=size(Dimension(i).Comb.NumberOfCells,2);
        for j=1:sz1
            nCells(j)=Dimension(i).Comb.NumberOfCells(j).NumberOfCellsChosen;
        end
    else
        sz1=size(nCells,2);
    end
    cmap=colormap(jet(sz1));
    for j=1:sz1
        pleg=sprintf('%d%s',nCells(j),' Cells');
        for k=1:szTh
            yMean(k)=Dimension(i).Comb.NumberOfCells(j).PerThresh(k).Mean;
            yStd(k)=Dimension(i).Comb.NumberOfCells(j).PerThresh(k).Std;
        end
        %Plot the Std errors
        P2=patch([EuThresh fliplr(EuThresh)], [yMean+yStd fliplr(yMean-yStd)], cmap(j,:),'EdgeColor','none');
        hold on
        P2.FaceAlpha=0.2;
        %Plot Mean
        plot(EuThresh,yMean,'Color',cmap(j,:),'LineWidth',2,'DisplayName',pleg)
    end
    legend
    xlim([0 max(EuThresh)])
    ylim([0 100])
    title('Percent Unique Per Trial')
    xlabel('Euclidean Distance Threshold')
    ylabel('Percent Unique Per Trial')
    savefig(h,fign,'compact');
    mnl_ExportEPSdense(h,fign);
end
end
%% Sub functions
function [PcUnique]=mnl_CalculateThePercentageUnique(SomaVals,EuThresh)
%Now measure the distances
[~,Eu_Matrix]=mnl_GroupColourEuclidean_ListAndMatrix(SomaVals);
%% For Each EuThresh
nThresh=size(EuThresh,2);
nSomas=size(SomaVals,1);
PcUnique=nan(1,nThresh);
for i=1:nThresh
    Thresh=EuThresh(i);
    %Calculate the percent unique
    UniqueCounter=0;
    for j=1:nSomas
        EuVals=Eu_Matrix(j,:);
        [~,loc]=find(EuVals<=Thresh);
        nPairs=find(loc~=j);
        if isempty(nPairs)==1
            UniqueCounter=UniqueCounter+1;
        end
    end
    PcUnique(i)=(UniqueCounter/nSomas)*100;
end
end