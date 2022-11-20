function [BkgMedian,MaxValues,Trace,Sig2NoisePeak]=mnl_ExtractVoxelValuesForTraces_v7(Trace,Data,BkgMedian,MagThresh,NumLengthThresh)
%New extraction method to calculate the mean cell values following this
%pathway. But now it also filters according to the brightness magnitude too
% 1) Visualise the Colour Distributions to see if the channels have
% a sufficient signal to noise ratio
% 2) Find the mean of the background voxels
% 3) Calculate the average colour of each segment (Cells.RawMean)
% 4) Subtract the background mean from the averaged segment (Cells.BkgRemovedMean)
% 5) Calculate the Maximum Value for each channel of the averaged traces
% 6) Normalise to this Maximum Value (Cells.NormMean)
% 7) Calculates the vector magnitude
% 8) Vector Normalise (Cells.VecNormMean)
% 9) Filter out the dim Traces (<MagThresh)
% 10) Barcode of each Process
%
%Inputs
% Trace - Structured Matrix with the Somatic Data containing
%           Points
%           Diameter
%           Original Trace
%           LengthVx
%           LengthReal
%           AllVoxels
%           Scale
%           Dim
% Data - The Image in its 4 dimensions (x,y,c,z)
% BkgMedian - The Background median values calculated beforehand. If not it
% will be calculated elsewhere
% MagThresh - The magnitude threshold
% NumLengthThresh - The minimum length threshold
%
%Outputs
% BkgMean - The mean background value of each colour
% MaxValues - The maximum value of each channel from the average colour of each soma
% Trace - Updated structure with the raw, background subtracted, normalised, and vector normalised mean colours of each cell
%
% Figures are also produced to evaluate the quality of each channel and to
% visualise the spread across all the dimensions

%% Basic Information
nTrace=size(Trace,2); %The number of somas
dims=size(Data); %Assumes data is [y x c z]
RawMatrix=nan(nTrace,dims(3));
BkgSubMatrix=nan(nTrace,dims(3));
NormMatrix=nan(nTrace,dims(3));
if isempty(BkgMedian)==1
    %% Step 1 - Visualise the Colour Distributions
    %First get all the Voxels for the Somas
    n=1; %Counter for the Voxel List
    for i=1:nTrace %For each Soma
        %Get the Voxel List of that Trace
        temp=Trace(i).AllVoxels;
        szT=size(temp,1);
        AllVx(n:n+szT-1,:)=temp;
        n=n+szT;
        clear temp szT
    end
    AllVx=unique(AllVx(:,1:3),'rows'); %Remove Duplicates
    nVx=size(AllVx,1); %Final Number of voxels
    sAllVx=sortrows(AllVx,3); %sorted along z
    %Make a BW Mask
    BW=zeros(dims(1),dims(2),dims(4));
    for i=1:dims(4)
        index=find(sAllVx(:,3)==i);
        BW(sAllVx(index,2),sAllVx(index,1),i)=1;
    end
    
    %Now the background locations
    n=1;
    for i=1:dims(4)
        tBW=BW(:,:,i);
        [r,c]=find(tBW==0);
        szRC=size(r,1);
        n2=n+szRC-1;
        z=ones(szRC,1)*i;
        AllBkg(n:n2,:)=[r,c,z];
        n=n+szRC;
        clear tBW r c szRC n2 z
    end
    nBkg=n-1;
    %Now sub-sample to speed it up
    [SampleN_Vx]=mnl_DetermineSampleSize(99,nVx,1); %Signal Subsampling
    [SampleN_Bkg]=mnl_DetermineSampleSize(99,nBkg,1); %Backgroud Subsampling
    Chosen_Vx=randperm(nVx,round(SampleN_Vx)); %Which ones are chosen
    Chosen_Bkg=randperm(nBkg,round(SampleN_Bkg)); %Which ones are chosen
    % For the Signal - NB for this and the background it is quicker to assign it directly but the loop is because the memory won't direct extract  more than a million or so values
    for i=1:round(SampleN_Vx)
        Pos=AllVx(Chosen_Vx(i),:);
        if Pos(1)>dims(2)
            Pos(1)=dims(2);
        elseif Pos(1)<1
            Pos(1)=1;
        end
        
        if Pos(2)>dims(1)
            Pos(2)=dims(1);
        elseif Pos(2)<1
            Pos(2)=1;
        end
        
        if Pos(3)>dims(4)
            Pos(3)=dims(4);
        elseif Pos(3)<1
            Pos(3)=1;
        end
        temp=Data(Pos(2),Pos(1),:,Pos(3)); %Remember it is flipped in the data
        Signal(i,:)=temp;
    end
    % For the Background
    for i=1:round(SampleN_Bkg)
        Pos=AllBkg(Chosen_Bkg(i),:);
        temp=Data(Pos(2),Pos(1),:,Pos(3));%Remember it is flipped in the data
        Background(i,:)=temp;
    end
    % So now we can evaluate the quality of the signal
    for i=1:dims(3) %for each channel
        [Sig2NoisePeak(i)]=mnl_EvaluateCumulativeSignal2Noise5(Background(:,i),Signal(:,i),i);
    end
    %% Step2 - Find the mean background value
    mBkg=median(Background);
    BkgMedian=mBkg;
else
    mBkg=BkgMedian;
end
%% Step 3 - Calculate the mean value of each Cell (Trace.RawMean)
for i=1:nTrace
    temp=round(Trace(i).AllVoxels); %Get the Voxels for that Trace
    %Filter the NaNs if they exist
    a=~isnan(temp(:,1));
    temp=temp(a,:);
    %Now get the colours
    szT=size(temp,1);
    Values=zeros(szT,dims(3));
    for j=1:szT
        if temp(j,3)>dims(4)
            temp(j,3)=dims(4);
        elseif temp(j,3)<1
            temp(j,3)=1;
        end
        if temp(j,2)>dims(1)
            temp(j,2)=dims(1);
        elseif temp(j,2)<1
            temp(j,2)=1;
        end
        if temp(j,1)>dims(2)
            temp(j,1)=dims(2);
        elseif temp(j,1)<1
            temp(j,1)=1;
        end
        Values(j,1:dims(3))=Data(temp(j,2),temp(j,1),:,temp(j,3));
    end
    mV=mean(Values,'omitnan');
    Trace(i).RawMean=mV;
    smV=mV-mBkg;
    %If it is a negative value convert it to 0
    for j=1:dims(3)
        if smV(j)<0
            smV(j)=0;
        end
    end
    %% Step 4 - Subtract the Background Mean Value (Trace.BkgRemovedMean)
    Trace(i).BkgRemovedMean=smV;
    %Assign to the Required Matrices
    RawMatrix(i,:)=mV;
    BkgSubMatrix(i,:)=smV;
    clear temp Values    
end
%% Step 5 - Find the Maximum Value of each channel
MxVals=max(BkgSubMatrix);
MaxValues=MxVals;
%% Steps 6 and 7 - Normalise each Channel to the nex maximum value then calculate the vector magnitude
for i=1:nTrace
    smV=Trace(i).BkgRemovedMean;
    nsmV=smV./MxVals;
    mnsmV=vectorNorm(nsmV);
    Trace(i).NormMean=nsmV;
    Trace(i).NormMeanMagnitude=mnsmV;
    NormMatrix(i,:)=nsmV;
    MagMatrix(i,1)=mnsmV;
end
%% Step 8 - Vector Normalise each Channel
[VecNormMatrix]=mnl_NormaliseVectors(NormMatrix);
for i=1:nTrace
    Vals=VecNormMatrix(i,:);
    if isnan(Vals)==1
        VecNormMatrix(i,:)=zeros(1,dims(3)); %If all the values are zero it switches from NaN to 0
    end
    Trace(i).VecNormMean=VecNormMatrix(i,:);
end
%% Step 9 - Now filter out the short and dim traces
%Add the length to the MagMatrix
for i=1:nTrace
    MagMatrix(i,2)=Trace(i).LengthReal;
end
idx=MagMatrix(:,1)>=MagThresh & MagMatrix(:,2)>=NumLengthThresh;
fTrace=Trace(idx);
fRawMatrix=RawMatrix(idx,:);
fBkgSubMatrix=BkgSubMatrix(idx,:);
fNormMatrix=NormMatrix(idx,:);
fVecNormMatrix=VecNormMatrix(idx,:);

%Clear old values
clear Trace RawMatrix BkgSubMatrix NormMatrix VecNormMatrix
%Provide Update the values
Trace=fTrace;
RawMatrix=fRawMatrix;BkgSubMatrix=fBkgSubMatrix;NormMatrix=fNormMatrix;VecNormMatrix=fVecNormMatrix;
%% Step 10 - Now Visualise these Colours as a Barcode
cmap=magma(256);
colormap(cmap);
figure('Name','Barcodes of the Mean Colours','Units','normalized','Position',[0.05 0.25 0.9 0.5])
%Raw Colours
subplot(1,4,1)
imagesc(RawMatrix)
colormap(cmap);
colorbar
ylabel('Trace #')
xlabel('Channel #')
title('Raw Values')
%Background Subtracted
subplot(1,4,2)
imagesc(BkgSubMatrix)
colormap(cmap);
colorbar
ylabel('Trace #')
xlabel('Channel #')
title('Bkg Subtracted Values')
%Normalised
subplot(1,4,3)
imagesc(NormMatrix)
colorbar
ylabel('Trace #')
xlabel('Channel #')
title('Normalised to Max Averaged Soma Colour')
%Vector Normalised
subplot(1,4,4)
imagesc(VecNormMatrix,[0 1])
colormap(cmap);
colorbar
ylabel('Trace #')
xlabel('Channel #')
title('Vector Normalised Values')
end