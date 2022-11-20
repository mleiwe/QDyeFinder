function [OptF1mean,OptThresh,OptF1median,OptThreshMedian]=mnl_OptimumEuThreshViaF1Star_MagWeighted(efPxTrace)
EuThreshVals=0.05:0.01:1;
szTh=size(EuThreshVals,2);
dim=efPxTrace(1).dim;
n=1;
for i=1:szTh
    EuThresh=EuThreshVals(i);
    [FinalClusters,FinalClusterIDs,~,~,~]=mnl_EuclideanCrawler_Mag_Weighted_v3(efPxTrace,EuThreshVals(i),dim,'n',[],1,0);
    %Per Trace
    %[Kappa(:,i),F1score(:,i),Recall(:,i),Precision(:,i),TruePositiveRate(:,i),FalsePositiveRate(:,i)]=mnl_EvaluateClassifierPerTrace(efPxTrace,fClusterIDs,fClusters,dim,'n');
    %Per Cell
    [Kappa(:,n),F1score(:,n),Recall(:,n),Precision(:,n),TruePositiveRate(:,n),FalsePositiveRate(:,n)]=mnl_EvaluateClassifierPerCell(efPxTrace,FinalClusterIDs,FinalClusters,dim,'n',EuThreshVals(i),1);
    Setting(n,1)=EuThreshVals(i);
    n=n+1;
end
%% Find the optimum value from F1 scores
F1scoreMean=mean(F1score,1,'omitnan');
F1scoreMedian=median(F1score,1,'omitnan');
F1scoreStd=std(F1score,1,'omitnan');
%Now find the optimum thresh from the mean and median
[M,I]=max(F1scoreMean);
[MedM,MedI]=max(F1scoreMedian);
maxI=MedI;
OptF1mean=M(1);
OptF1median=MedM(1);
OptThresh=EuThreshVals(I(1)); %Defaults to the lowest value
OptThreshMedian=EuThreshVals(MedI(1)); %Defaults to the lowest value
%% ROC curves
figure('Name','ROC curves')
nRows=size(FalsePositiveRate,1);
mFPR=mean(FalsePositiveRate,1,'omitnan');
mTPR=mean(TruePositiveRate,1,'omitnan');
scatter(mFPR,mTPR,'.k')
hold on
xlim([0 1])
xlabel('False Positive Rate')
ylim([0 1])
ylabel('True Positive Rate')
%Plot the fitted line
[smFPR,I]=sort(mFPR);%sort FPR
smTPR=mTPR(I);%sort TPR by FPR values
plot(smFPR,smTPR,'r');
AUC=round(trapz(smFPR,smTPR),3);
str=sprintf('%s%s','AUC = ',num2str(AUC));
text(0.8,0.1,str)
OptPoint_Mean=scatter(mFPR(maxI(1)),mTPR(maxI(1)),20,"g");
text(mFPR(maxI(1))+0.05,mTPR(maxI(1)),"<--- Optimum Threshold (via F1 score)")
hold off
%% Precision Recall curve
%Now plot the precision recall curve
figure('Name','Precision Recall Curve')
for i=1:nRows
    plot(Recall(i,:),Precision(i,:),"Color",[0.5 0.5 0.5])
    hold on
end
mRecall=mean(Recall,1,'omitnan');
mPrecision=mean(Precision,1,'omitnan');
plot(mRecall,mPrecision,'r')
xlim([0 1])
xlabel('Recall')
ylim([0 1])
ylabel('Precision')
scatter(mRecall(maxI(1)),mPrecision(maxI(1)),20,"g");
text(mRecall(maxI(1))-0.7,mPrecision(maxI(1)),"Optimum Threshold (via F1 score) --->")
title('Precision Recall curve')
%% For F1 Score plot
nCells=size(F1score,1);
figure('Name','F1 Score Means')
errorbar(EuThreshVals,F1scoreMean,F1scoreStd,'k')
hold on
plot(EuThreshVals,F1scoreMean,'r')

figure('Name','F1 Score Medians')
%Scatter the individual cells
c=[0.5 0.5 0.5];
for i=1:szTh
    scatter(EuThreshVals(i),F1score(:,i),2,c)
    hold on
end
%Plot the median and inter quartile range
F1iqr=iqr(F1score);
errorbar(EuThreshVals,F1scoreMedian,F1iqr,'k')
plot(EuThreshVals,F1scoreMedian,'r')

figure('Name','F1 Score Boxplots')
for i=1:szTh
    ThreshLegend{i}=num2str(EuThreshVals(i));
end
mnl_boxplot(F1score,ThreshLegend,'F1 Score');
end