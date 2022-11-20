function [LengthThresh]=mnl_EvaluateFragmentsColourConsistency5(FragmentValues,ColourConsistencyThreshold)
% Attempt to explore how consistent the colour is for a trace.
%
%Inputs
% FragmentValues - Structure created by mnl_SplitTracesIntoFragments
%
%Outputs
% LengthThresh - The length at which the colour is sufficiently consistent
%
%Marcus Leiwe, Kyushu University - 24th December 2020
%% Step One - Calculate the Decay of each trace
%Organise each calculation by it's trace
NumTraces=size(FragmentValues,2);
DecayPoints=struct('Points',[],'Length',[],'EuDistance',[]);
for i=1:NumTraces
    nJ=size(FragmentValues(i).FragSize,2);
    n=1;
    for j=1:nJ
        nK=size(FragmentValues(i).FragSize(j).EuDistance,2);
        for k=1:nK
            x_Points(n)=j;
            x_Length(n)=FragmentValues(i).FragSize(j).Length(k);
            y(n)=FragmentValues(i).FragSize(j).EuDistance(k);
            n=n+1;
        end        
    end
    DecayPoints(i).Points=x_Points;
    DecayPoints(i).Length=x_Length;
    DecayPoints(i).EuDistance=y;
    clear x_Points x_Length y 
end
% Fit Decay Curves
MaxLength=0;
CoEffMatrix_Length=nan(NumTraces,2);
for i=1:NumTraces
    %For Points
    X=DecayPoints(i).Points;
    Y=DecayPoints(i).EuDistance;
    if size(X,2)>10 %Need at least 10 points
        %For Length
        X=DecayPoints(i).Length;
        beta0=[1,0.01];
        [coefficients,modelfunc]=mnl_FitExponentialDecay_v2(X',Y',beta0);%y=CoEff(1)*exp(-CoEff(2)*X)
        %[coefficients]=mnl_FitExponentialDecay([X' Y'],beta0); %y=CoEff(1)*exp(-CoEff(2)*X)+CoEff(3)
        CoEffMatrix_Length(i,:)=coefficients;
        mxLV=max(DecayPoints(i).Length);
        if mxLV>MaxLength
            MaxLength=mxLV;
        end
    end
end
%Now plot them
figure('Name','Fitted Decay Curves')
X=linspace(0,MaxLength,201)';
MedianCoEffs=median(CoEffMatrix_Length,'omitnan');
for i=1:NumTraces
    %tY=CoEffMatrix_Length(i,1)*exp(-CoEffMatrix_Length(i,2)*X)+CoEffMatrix_Length(i,3);
    if ~isnan(CoEffMatrix_Length(i,:))
        tY=modelfunc(CoEffMatrix_Length(i,:),X);
        plot(X,tY,'k','LineWidth',0.25)
        hold on
    end
    clear tY
end
%tY=MedianCoEffs(1)*exp(-MedianCoEffs(2)*X)+MedianCoEffs(3);
tY=modelfunc(MedianCoEffs,X);
plot(X,tY,'r','LineWidth',2)
ylim([0 0.5])
title('Per Length of Fragment in Microns')
% Then choose which magnitude is the lowest
LengthThresh=(log(ColourConsistencyThreshold/MedianCoEffs(1)))/(MedianCoEffs(2)*-1);
text(LengthThresh+0.1,ColourConsistencyThreshold,sprintf('%s%d','Minimum Consistent Length = ',LengthThresh))
savefig('MinimumLengthOfTracePlot');
end
