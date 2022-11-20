function [Cells]=mnl_GeneratePossionsNchannels(NumPoints,Spreads,Ndim)
% Code for creating  N dimensional poisson distributions at
% user-specified spreads.
%
%Inputs
% NumPoints=number of cells
% Spread=[pdf variances] e.g. [0 0.1 0.2 0.5 1 2 4 8]
% Ndim - Number of dimensions
%
%Outputs
% Cells=Structure containing cells and the spread information
% Structure Details - Copy Number = The specified poisson spread
%                   - XFPvals = n*Ndim 
%                   - NormXFPvals = vector normalised values (n*dim)
%% Generate Poisson Matrices
maxSp=max(Spreads);
if maxSp*2>50
    MaxXFP=maxSp*2;
else
    MaxXFP=50;
end
x=linspace(0,MaxXFP,MaxXFP+1); %Theoretical maximum number of "fluorescent proteins" in the distribution
sz=size(Spreads);
for i2=1:sz(2)
    for i=1:MaxXFP
        MatlabPoisson(i,i2)=poisspdf(x(i),Spreads(i2));
    end
    legnames{i2}=sprintf('%s%d','Copy Number ',Spreads(i2));
end
% figure
% plot(x,MatlabPoisson)
% legend(legnames);
%% Allocate Values for each cell
for i=1:sz(2) % For each Vector Conc
    a(:,i)=MatlabPoisson(:,i).*NumPoints;
    for j=1:Ndim %for each colour
        RealNum=0;
        Counter=1;
        for k=1:MaxXFP %Number of Copies of XFP in the cell - Theoretical Max needs to be very unlikely
            if round(MatlabPoisson(k,i)*NumPoints)>=1 %And if the value will be sufficient
                if RealNum<NumPoints %If there are more cells to label
                    NumCopies=round(MatlabPoisson(k,i)*NumPoints); %The number of times this value will appear
                    t2=ones(NumCopies,1)*k-1; %Say there are k-1 XFPs being produced, first one is zero
                    Temp(Counter:(Counter+NumCopies-1))=t2;
                    %                     for m=1:round(MatlabPoisson(k,i)*NumPoints)
                    %                         Temp(Counter,j)=k-1;
                    %                         Counter=Counter+1;
                    %                         RealNum=RealNum+1;
                    %                     end
                    Counter=Counter+NumCopies;
                    RealNum=RealNum+NumCopies;
                end
            else
                chk=1;
            end
        end
        if Counter<=NumPoints %If at the end there is still some cells missing....
            Temp(Counter:NumPoints)=0;%then give them a zero
        elseif RealNum>NumPoints %Or if the number of cells has gone over the limit
            Temp=Temp(1:NumPoints);
        end
        Cells(i).CopyNumber=Spreads(i);
        Cells(i).XFPvals(:,j)=Temp(randperm(length(Temp))); %Shuffle the expression        
    end
    clear Temp
    data=Cells(i).XFPvals;
    Cells(i).NormXFPvals=mnl_NormaliseVectors(data);
    %Switch NaNs to Zeros
    data=Cells(i).NormXFPvals;
    data(isnan(data))=0;
    Cells(i).NormXFPvals=data;
    clear data    
end
%b=round(a);
% %% Statistical Evaluation
% n=1;
% for i=1:sz(2)
%     data=Cells(i).XFPvals;
%     [data]=mnl_NomaliseVectors(Cells(i).XFPvals);% Normalise Vectors
%     %data=mnl_RatioVectors(Cells(i).RGB);
%     %data=mnl_RatioVectors2(Cells(i).RGB);
%     data=Convert2Ratios(Cells(i).XFPvals);
%     [CC,CClist]=mnl_GroupColourCrossCorr(data);% Cross Correlation
%     [EuD_all,~,EuD_mean,EuD_allMean]=mnl_GroupColourEuclidean(data); %Euclidean Distances
%     Cells(i).CC=CC;
%     Cells(i).CClist=CClist;
%     Cells(i).EuD_all=EuD_all;
%     Cells(i).EuD_mean=EuD_mean;
%     Cells(i).EuD_All_Mean=EuD_allMean;
%     sz3=size(CClist);
%     if n<sz3(1)
%         n=sz3(1);
%     end
% end
% CCcombined=NaN(n,sz(2));
% EuD_all_combined=NaN(n,sz(2));
% sz2=size(Cells(1).EuD_mean);
% EuD_mean_combined=NaN(sz2(1),sz(2));
% for i=1:sz(2)
%     sz2=size(Cells(i).CClist,1);
%     sz3=size(Cells(i).EuD_all,1);
%     sz4=size(Cells(i).EuD_mean,1);
%     CCcombined(1:sz2,i)=Cells(i).CClist;
%     EuD_all_combined(1:sz3,i)=Cells(i).EuD_all;
%     EuD_mean_combined(1:sz4,i)=Cells(i).EuD_mean;
% end
% mnl_boxplot(CCcombined,{'0' '0.1' '0.2' '0.5' '1' '2' '4' '6' '8' '10' '20'},'Cross Correlation Value');% Stats Graph
% mnl_boxplot(EuD_all_combined,{'0' '0.1' '0.2' '0.5' '1' '2' '4' '6' '8' '10' '20'},'Euclidean Distance');% Stats Graph
% mnl_boxplot(EuD_mean_combined,{'0' '0.1' '0.2' '0.5' '1' '2' '4' '6' '8' '10' '20'},'Euclidean Distance');% Stats Graph
% figure('Name','CC measurements')
% mnl_CumulativePlot3(Cells(1).CClist,Cells(2).CClist,Cells(3).CClist,Cells(4).CClist,Cells(5).CClist,Cells(6).CClist,Cells(7).CClist,Cells(8).CClist,Cells(9).CClist,Cells(10).CClist,Cells(11).CClist)
% figure('Name','Euclidean Distance All')
% mnl_CumulativePlot3(Cells(1).EuD_all,Cells(2).EuD_all,Cells(3).EuD_all,Cells(4).EuD_all,Cells(5).EuD_all,Cells(6).EuD_all,Cells(7).EuD_all,Cells(8).EuD_all,Cells(9).EuD_all,Cells(10).EuD_all,Cells(11).EuD_all)
% figure('Name','Euclidean Distance to mean')
% mnl_CumulativePlot3(Cells(1).EuD_mean,Cells(2).EuD_mean,Cells(3).EuD_mean,Cells(4).EuD_mean,Cells(5).EuD_mean,Cells(6).EuD_mean,Cells(7).EuD_mean,Cells(8).EuD_mean,Cells(9).EuD_mean,Cells(10).EuD_mean,Cells(11).EuD_mean)
end
function [data]=mnl_RatioVectors(Group)
maxVals=max(Group);
data=Group./maxVals;
end
function [data]=mnl_RatioVectors2(Group)
sz=size(Group);
data=zeros(sz);
for i=1:sz(1)
    MaxVal=max(Group(i,:));
    data(i,:)=Group(i,:)/MaxVal;
    if isnan(data(i,:))==1
        data(i,:)=0;
    end
end
end
function [nGroup]=Convert2Ratios(Group)
%Function to convert the copy numbers to a bespoke ratio for each copy
%number set
sz=size(Group);
nGroup=zeros(sz);
for i=1:sz(1)
    base=sum(Group(i,:));
    if base~=0
        nGroup(i,:)=Group(i,:)/base;
    else
        nGroup(i,:)=zeros(1,sz(2));
    end
end
end
