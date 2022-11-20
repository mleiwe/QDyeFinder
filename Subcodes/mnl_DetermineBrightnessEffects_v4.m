function [BrightnessThresh]=mnl_DetermineBrightnessEffects_v4(FragmentValues,FragLength,EuclideanThreshold,ErrorThresh)
%This is the version 2 of determining what is an appropriate brightness
%threshold to use based on the data provided. It was initially designed to
%be compatible with mnl_ProcessTracesFromNeurolucida_DendritesSomas_v13. And
%as such relies on the Traces already being subdivided into fragments of
%the minimum fragment length with those Euclidean distances being measured.
%(i.e. use codes mnl_SplitTracesIntoFragments or mnl_CalculateMaxValues.
%This is an improvement of v3 by streamlining the decay functions etc...

%% Pre-Allocation Period
nTraces=size(FragmentValues,2);
%% Step One - Scroll Through the fragment values section to find the dividing
%point where there are the most fragments at the specified FragLength. Then
%record the vector magnitude and the percentage of traces that are within
%the EuThresh
Splits=struct('PointNum',[],'Length',[],'NormMeanMagnitude',[],'FullTraceMagnitude',[],'EuDistance',[]);
for i=1:nTraces
    FullLength=FragmentValues(i).All.Length;
    if FullLength>=FragLength
        nDiffLengths=size(FragmentValues(i).FragSize,2);
        nFrag=nan(nDiffLengths,1);
        for j=1:nDiffLengths
            Lengths=FragmentValues(i).FragSize(j).Length;
            idx=find(Lengths>=FragLength);
            nFrag(j)=size(idx,2);
        end
        [~,I]=max(nFrag);
        Splits(i).PointNum=FragmentValues(i).FragSize(I).PointNum;
        Splits(i).Length=FragmentValues(i).FragSize(I).Length;
        Splits(i).NormMeanMagnitude=FragmentValues(i).FragSize(I).NormMeanMagnitude;
        Splits(i).FullTraceMagnitude=FragmentValues(i).All.NormMeanMagnitude;
        Splits(i).EuDistance=FragmentValues(i).FragSize(I).EuDistance;
    else
        Splits(i).PointNum=NaN;
        Splits(i).Length=NaN;
        Splits(i).NormMeanMagnitude=NaN;
        Splits(i).FullTraceMagnitude=FragmentValues(i).All.NormMeanMagnitude;
        Splits(i).EuDistance=NaN;
    end
end
%% Step Two - Plot percentage of fragments within the threshold (y) and vector magnitude (x)
figure('Name','Plotting Individual Curves')
[cmap]=mnl_GenerateShuffledColourmap(nTraces);
AllX=[];
AllY=[];
Xvals=nan(nTraces,201);
Yvals=nan(nTraces,201);
for i=1:nTraces
    if ~isnan(Splits(i).EuDistance)
        X(:,1)=Splits(i).NormMeanMagnitude;
        Y(:,1)=Splits(i).EuDistance;
        scatter(X,Y,2,cmap(i,:),'filled')
        hold on
        %How many are within the Threshold?
        idx=find(Y<=EuclideanThreshold);
        nY=size(Y,1);
        n=size(idx,1);
        %Add extra info
        Splits(i).PcWithinThresh=(n/nY)*100;    
        AllX=[AllX;X];
        AllY=[AllY;Y];
        clear X Y n nY
    else
        Splits(i).PcWithinThresh=NaN;
    end
end
a=xlim;
%Fit mean
beta0=[1,0.01];
[b,modelfunc]=mnl_FitExponentialDecay_v2(AllX,AllY,beta0);%y=CoEff(1)*exp(-CoEff(2)*X)
mXvals=linspace(a(1),a(2),201)';
mYvals=modelfunc(b,mXvals);
plot(mXvals,mYvals,'r','LineWidth',2)
plot([a(1) a(2)],[EuclideanThreshold EuclideanThreshold],'--k')
xlabel('Normalised Mean Magnitude')
ylabel('Euclidean Distance to Whole Fragment')
%% Now express this as a percent error
nPoints=501;
Xints=linspace(0,max(AllX(:)),nPoints)';
%Now find the error rate at each threshold
ErrX=nan(nPoints,1);
ErrY=nan(nPoints,1);
for i=1:nPoints
    %Error rate
    idx=AllX>Xints(i);
    tY=AllY(idx);
    idx2=tY>EuclideanThreshold;
    ErrY(i)=sum(idx2)/sum(idx)*100;
end
%Now plot it
figure('Name','Error rate at each magnitude')
scatter(Xints,ErrY,'.')
%Fit a decay curve
%[coefficients]=mnl_FitExponentialDecay([Xints ErrY],beta0);
[b2,modelfunc]=mnl_FitExponentialDecay_v2(Xints,ErrY,beta0);%y=CoEff(1)*exp(-CoEff(2)*X)
mXvals=linspace(a(1),a(2),201)';
mYvals=modelfunc(b2,mXvals);
hold on
plot(mXvals,mYvals,'r','LineWidth',1)
plot([0 max(AllX(:))],[ErrorThresh ErrorThresh],'--k')
xlabel('Normalised Mean Magnitude')
ylabel('Error Rate(%)')
% Then choose which magnitude is sufficient
BrightnessThresh=(log(ErrorThresh/b2(1)))/(b2(2)*-1);
text(BrightnessThresh+0.2,ErrorThresh+1,sprintf('%s%d','Minimum Magnitude Threshold = ',BrightnessThresh))
savefig('MinimumMagnitudeOfTracePlot');
end