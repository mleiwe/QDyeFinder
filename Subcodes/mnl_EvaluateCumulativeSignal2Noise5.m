function [Sig2NoisePeak,NoisePercentiles,SignalPercentiles]=mnl_EvaluateCumulativeSignal2Noise5(Noise,Signal,ChanNum)
%mnl_EvaluateCumulativeSignal2Noise is designed to evaluate where the
%difference is between two different distributions of different sample
%sizes. It then plots three graphs to visualise the differences. 
% The first plot is a boxplot, the second is the cumulative distribution, the third is the signal to noise ratio at each percentil
% These are then saved as a MATLAB figure and a png with the format "Sig2Noise_ChanX"

%Inputs
% Noise - Single column of data as baseline
% Signal - Single column of data to compare it
%Outputs
% NoisePercentiles - Cumulative distribution percentiles for the background
% SignalPercentiles - Cumulative distribution percentiles for the signal
% Sig2NoisePeak - the peak signal to noise ratio between the 81st and 100th percentile
%
% Written by Marcus Leiwe, Kyushu University, 2021
%% Basic Info
fn=sprintf('%s%d','Value Distributions for Channel ',ChanNum);
szBkg=size(Noise,1);
szSig=size(Signal,1);
h=figure('Name',fn,'Units','normalized','Position',[0.05 0.33 0.9 0.33]);
%% Step One - Make a box plot for the signal vs the noise
bigN=max([szBkg szSig]);
DataSet=nan(bigN,2);
DataSet(1:szBkg,1)=Noise;
DataSet(1:szSig,2)=Signal;
subplot(1,3,1)
mnl_boxplot2(DataSet,{'Noise','Signal'},'Intensity Value','y','n')
%% Step Two - Make the Cumulative Density Plot
subplot(1,3,2)
mnl_CumulativePlot3(Noise,Signal);
legend('Background','Signal')
xlabel('Intensity')
ylim([50 100])
title('Cumulative Density Plot')
%% Step Three - Difference Between Intensity Values
[Snoise,Ynoise]=mnl_SortNoNaNs(Noise);%Sort Background and make cumulative percentages
[Ssignal,Ysignal]=mnl_SortNoNaNs(Signal);%Sort Signal and make cumulative percentages
% Now find the values at the various percentiles
Xcumulative=0:1:100;
sz=size(Xcumulative,2);
%Pre-allocate the Noise and Signal Percentiles
NoisePercentiles=nan(sz,2);
SignalPercentiles=nan(sz,2);
S2N=nan(sz,1);
for i=1:sz
    PctVal=Xcumulative(i);
    %For Noise
    [~,idx]=min(abs(Ynoise-PctVal)); %The value closest to the percentile
    tNoise=Snoise(idx);
    NoisePercentiles(i,:)=[PctVal tNoise];
    %For Signal
    [~,idx]=min(abs(Ysignal-PctVal));
    tSignal=Ssignal(idx);
    SignalPercentiles(i,:)=[PctVal tSignal];
    % Signal to Noise Ratio at this Point
    S2N(i)=tSignal/tNoise;
end
subplot(1,3,3)
plot(S2N,Xcumulative)
xlabel('Signal / Noise')
ylabel('Cumulative Percent')
title('Signal to Noise Ratio')
%% Calculate the signal to noise peak (the second peak)
List80to100=S2N(81:101);%Just look at the results from the 80th percentile onwards
idx=List80to100~=Inf;
FiltList=List80to100(idx);
if sum(idx)>1
    Sig2NoisePeak=max(FiltList);
else
    Sig2NoisePeak=5;
    disp("Warning, there is no background value so the signal to noise values are all infinite, replacing with 5")
end
%% Save the figure
sn=sprintf('%s%d','Sig2Noise_Chan',ChanNum);
savefig(h,sn);
saveas(h,sn,'png');
end
function [Sdata,Ydata]=mnl_SortNoNaNs(Data)
% pData=The cumulative percentage of the data
s1=length(Data(~isnan(Data))); 
if s1~=0
    %Removing NaNs
    Sdata=sort(Data);
    Sdata=Sdata(1:s1,1); %Removes the NaNs at the end
    %Get Y values (cumulative % values)
    index=~isnan(Sdata); %Index of real values
    TotVals=sum(index);
    n=1;
    for i2=1:s1(1)
        if index(i2)==1
            Ydata(n)=n*(1/TotVals)*100;
            n=n+1;
        end
    end
    clear index
end
end