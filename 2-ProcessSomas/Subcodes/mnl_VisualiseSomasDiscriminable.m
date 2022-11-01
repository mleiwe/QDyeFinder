function mnl_VisualiseSomasDiscriminable(Somas)
EuThresh=0:0.01:0.5;
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

%% For Each Permutation
for i=2:maxDim
    nComb=size(Perms(i).PermCombo,1);
    fprintf('%s%d%s\n','For ',i,' dimensions...')
    for j=1:nComb
        Comb=Perms(i).PermCombo(j,:);
        % Sub Select Depending on the Permutations
        tSomaVals=SomaVals(:,Comb);
        %Eval
        [Discrim_Pc,Discrim_PcStd,PcUnique]=mnl_EvaluateThePermutation(tSomaVals,EuThresh);
        Dimension(i).Comb(j).Discrim_Pc=Discrim_Pc;
        Dimension(i).Comb(j).Discrim_PcStd=Discrim_PcStd;
        Dimension(i).Comb(j).Matrix=tSomaVals;
        Dimension(i).Comb(j).PcUnique=PcUnique;
        mnl_InsertProgressTrackerInLoops(j,nComb)
    end
end
%% Plot Unique
figure
cmap=colormap(lines(maxDim));
for i=2:maxDim
    sz=size(Dimension(i).Comb,2);
    for j=1:sz
        Vals(j,:)=Dimension(i).Comb(j).PcUnique;
    end
    mVals=mean(Vals,1);
    stVals=std(Vals,0,1);
    %Plot plus and minus StDev
    P2=patch([EuThresh fliplr(EuThresh)], [mVals+stVals fliplr(mVals-stVals)], cmap(i,:));
    hold on
    P2.FaceAlpha=0.2;
    %Plot Mean
    plot(EuThresh,mVals,'Color',cmap(i,:),'LineWidth',2)
    clear Vals mVals stVals
end
ylim([0 100])
xlim([0 max(EuThresh)])
title('Percent Unique')
xlabel('Euclidean Distance Threshold')
ylabel('Percent Unique')
%% Plot Discrim
figure
cmap=colormap(lines(maxDim));
for i=2:maxDim
    sz=size(Dimension(i).Comb,2);
    for j=1:sz
        Vals(j,:)=Dimension(i).Comb(j).Discrim_Pc;
        Vals_Std(j,:)=Dimension(i).Comb(j).Discrim_PcStd;
    end
    mVals=mean(Vals,1);
    stVals=mean(Vals_Std,1);
    %Plot plus and minus StDev
    P2=patch([EuThresh fliplr(EuThresh)], [mVals+stVals fliplr(mVals-stVals)], cmap(i,:));
    hold on
    P2.FaceAlpha=0.2;
    %Plot Mean
    pleg=sprintf('%d%s',i,' XFPs');
    plot(EuThresh,mVals,'Color',cmap(i,:),'LineWidth',2,'DisplayName',pleg);
    clear Vals mVals stVals
end
xlim([0 max(EuThresh)])
ylim([0 100])
title('Percent Discriminable')
xlabel('Euclidean Distance Threshold')
ylabel('Percent Discriminable')
legend
%The Original Code
%[Dimension2]=mnl_OriginalCode(Dimension,EuThresh,Perms);i
end
%% Sub-codes
function [Discrim_Pc,Discrim_PcStd,PcUnique]=mnl_EvaluateThePermutation(SomaVals,EuThresh)
%Now measure the distances
[EuD_All,Eu_Matrix]=mnl_GroupColourEuclidean_ListAndMatrix(SomaVals);
%% For Each EuThresh
nThresh=size(EuThresh,2);
nSomas=size(SomaVals,1);
PcUnique=nan(1,nThresh);
for i=1:nThresh
    Thresh=EuThresh(i);
    nSame=0;
    nDiff=0;
    %Calculate the percent discriminable
    szAll=size(EuD_All,1);
    [r,~]=find(EuD_All<=Thresh);
    nSame=size(r,1);
    nDiff=szAll-nSame;
    %Discrim_Pc(i)=(nDiff/szAll)*100;
    %Calculate the percent unique
    UniqueCounter=0;
    for j=1:nSomas
        EuVals=Eu_Matrix(j,:);
        [~,loc]=find(EuVals<=Thresh);
        nPairs=find(loc~=j);
        szNP=size(nPairs,2);
        if isempty(nPairs)==1
            UniqueCounter=UniqueCounter+1;
            Discrim_Percent(j)=100;
        else
            Discrim_Percent(j)=((nSomas-1-szNP)/(nSomas-1))*100;
        end
    end
    Discrim_Pc(i)=mean(Discrim_Percent);
    Discrim_PcStd(i)=std(Discrim_Percent);
    PcUnique(i)=(UniqueCounter/nSomas)*100;
    clear Discrim_Percent
end
end