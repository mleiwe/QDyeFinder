function [Sig2NoisePeak]=mnl_EvaluateCumulativeSignal2Noise(Noise,Signal,ChanNum)
%mnl_EvaluateCumulativeSignal2Noise is designed to evaluate where the
%difference is between two different distributions of different sample
%sizes. It then plots several graphs to visualise the differences

%Inputs
% Noise - Single column of data as baseline
% Signal - Single column of data to compare it
%Outputs
% As of yet undecided, will decide later
%
% Written by Marcus Leiwe, Kyushu University 2020
%% Basic Info
fn=sprintf('%s%d','Value Distributions for Channel ',ChanNum);
szBkg=size(Noise,1);
szSig=size(Signal,1);
maxSig=max(Signal);
maxBkg=max(Noise);
figure('Name',fn,'Units','normalized','Position',[0.05 0.33 0.9 0.33])

%% Step One - Make the Cumulative Density Plot
subplot(1,2,1)
mnl_CumulativePlot3(Noise,Signal);
legend('Background','Signal')
xlabel('Intensity')
ylim([50 100])
title('Cumulative Density Plot')

%% Step 4 - Difference Between Intensity Values
[Snoise,Ynoise]=mnl_SortNoNaNs(Noise);%Sort Background and make cumulative percentages
[Ssignal,Ysignal]=mnl_SortNoNaNs(Signal);%Sort Signal and make cumulative percentages
% Now find the values at the various percentiles
Xcumulative=0:1:100;
sz=size(Xcumulative,2);
for i=1:sz
    PctVal=Xcumulative(i);
    %For Noise
    [~,idx]=min(abs(Ynoise-PctVal)); %The value closest to the percentile
    tNoise=Snoise(idx);
    %For Signal
    [~,idx]=min(abs(Ysignal-PctVal));
    tSignal=Ssignal(idx);
    % Signal to Noise Ratio at this Point
    S2N(i)=tSignal/tNoise;
end
subplot(1,2,2)
plot(S2N,Xcumulative)
xlabel('Signal / Noise')
ylabel('Cumulative Percent')
title('Signal to Noise Ratio')
%% Calculate the signal to noise peak (the second peak)
Sig2NoisePeak=max(S2N(81:101));
end

function [Sdata,Ydata]=mnl_SortNoNaNs(Data)
% pData=The cumulative percentage of the data
s1=length(Data(~isnan(Data))); 
if s1~=0
    %Removing NaNs
    Sdata=sort(Data);
    Sdata=Sdata(1:s1,1); %Removes the NaNs at the end
    %Get Y values (cumulative % values)
    index=[];
    index=isnan(Sdata)~=1; %Index of real values
    TotVals=sum(index);
    n=1;
    for i2=1:s1(1)
        if index(i2)==1
            Ydata(n)=n*(1/TotVals)*100;
            n=n+1;
        end
    end
end
end