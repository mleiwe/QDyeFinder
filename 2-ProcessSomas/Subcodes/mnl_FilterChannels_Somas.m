function [BadChannels,SampleBkg,SampleVx]=mnl_FilterChannels_Somas(Somas,Data,Scale)
% Function to automatically filter out the bad channels
%
% Inputs
%  Trace - the Structured Trace extracted from "mnl_ReadXML_BranchesSeperately_WhileLoop"
%  Scale - the Scale of the image
%  Data - the chromatically corrected image data
%
% Outputs
%  Bad Channels - a list of the bad channels
%  SampleBkg - A chosen sample of background voxels
%  SampleVx - A chose sample of trace voxels

%% Extract the Raw Colour Values Of the Traces
szT=size(Somas,2);
dim=size(Data);
n=1; %counter for collating all the voxels
AllVx=[];
for i=1:szT %Do this for each soma
    %% Extract the Mask for the Whole Fragment
    NumPoints=size(Somas(i).dxyz,1);
    for j=1:NumPoints
        %Convert the surface points to voxels
        dxyz(j,1)=round(Somas(i).dxyz(j,1));
        dxyz(j,2)=(round(Somas(i).dxyz(j,2)/Scale(1)))+1; %Add one to compensate
        dxyz(j,3)=(round(Somas(i).dxyz(j,3)/Scale(2))*-1)+1;%Flip the Y axis, add one to compensate
        dxyz(j,4)=(round(Somas(i).dxyz(j,4)/Scale(3))*-1)+1;%Flip the Z axis, add one to compensate
        if dxyz(j,1)<1
            dxyz(j,1)=1;
        end
        if dxyz(j,2)>dim(1)
            dxyz(j,2)=dim(1);
        elseif dxyz(j,2)<1
            dxyz(j,2)=1;
        end
        if dxyz(j,3)>dim(2)
            dxyz(j,3)=dim(2);
        elseif dxyz(j,3)<1
            dxyz(j,3)=1;
        end
        if dxyz(j,4)>dim(4)
            dxyz(j,4)=dim(4);
        elseif dxyz(j,4)<1
            dxyz(j,4)=1;
        end
    end
    VxList=dxyz(:,2:4);
    AllVx=[AllVx;VxList];
    clear VxList
end
%% Calculate the Background Values
%Find the Voxels of Interest
AllVx=unique(AllVx,'rows'); %Remove Duplicates
nVx=size(AllVx,1); %Final Number of voxels

% Now the background locations
prompt='Do you want to manually select the background? (y/n)';
UserInput=input(prompt,'s');
ManBkg=strcmp(UserInput,'y');
if ManBkg==1
    %Make a MIP from the Data
    for i=1:dim(1)
        for j=1:dim(2)
            a(1:dim(3),1:dim(4))=Data(i,j,:,:);
            MIP(i,j)=max(a(:));
        end
    end    
    %Draw an ROI around the background of the MIP
    figure('Name','Please Draw an ROI of the background')
    imagesc(MIP)
    title('Grayscale MIP')
    [x,y,Bkg_BW,xi,yi]=roipoly; %returns the XData and YData in x and y, the mask image in BW, and the polygon coordinates in xi and yi.
    %Now assign the locations to an AllBkg list
    [r,c]=find(Bkg_BW==1);
    nVx=size(r,1);
    nBkg=1;
    for i=1:nVx
        AllBkg(nBkg:nBkg+dim(4)-1,1)=r(i);
        AllBkg(nBkg:nBkg+dim(4)-1,2)=c(i);
        AllBkg(nBkg:nBkg+dim(4)-1,3)=1:1:dim(4);
        nBkg=nBkg+dim(4);
    end
    nBkg=nBkg-1;
else
    disp('Finding the Background...')
    TotVx=dim(1)*dim(2)*dim(4);%Total Number of Voxels of Interest
    nBkgVx=TotVx-nVx; %Number of background voxels
    [SampleN_Bkg]=mnl_DetermineSampleSize(95,nBkgVx,1); %Backgroud Subsampling
    xyzdim=[dim(2),dim(1),dim(4)];
    [Chosen_Bkg]=mnl_SubSampleBkgFromMask(round(SampleN_Bkg),xyzdim,AllVx);
end
%Now sub-sample all voxels to speed it up too
[SampleN_Vx]=mnl_DetermineSampleSize(99,nVx,1); %Signal Subsampling
Chosen_Vx=randperm(nVx,round(SampleN_Vx)); %Which ones are chosen
%% Now Compare the channels
% For the Signal - NB for this and the background it is quicker to assign it directly but the loop is because the memory won't direct extract  more than a million or so values
for i=1:round(SampleN_Vx)
    Pos=AllVx(Chosen_Vx(i),:);
    SampleVx(i,:)=Pos;
    idx=find(Pos==0);
    if isempty(idx)==0
        Pos(idx)=1;
    end
    temp=Data(Pos(2),Pos(1),:,Pos(3)); %Remember it is flipped in the data
    Signal(i,:)=temp;
end
% For the Background
for i=1:round(SampleN_Bkg)
    Pos=Chosen_Bkg(i,:);
    SampleBkg(i,:)=Chosen_Bkg(i,:);
    temp=Data(Pos(2),Pos(1),:,Pos(3));%Remember it is flipped in the data
    Background(i,:)=temp;
end
% Plot the differences
figure('Name','Scatter Plots to look for correlation')
c=1;
for i=1:dim(3)
    for j=1:dim(3)
        subplot(dim(3),dim(3),c)
        scatter(Background(:,j),Background(:,i),'.k')
        hold on
        scatter(Signal(:,j),Signal(:,i),'.r')
        c=c+1;
        
    end
end
% So now we can evaluate the quality of the signal
prompt='Do you want to filter the channels manually (y/n)';
UserInput=input(prompt,'s');
Manual=strcmp(UserInput,'y');
cGood=1;
cBad=1;
Sig2NoiseThresh=3;
BadChannels=[];
GoodChannels=[];
for i=1:dim(3) %for each channel
    [Sig2NoisePeak(i)]=mnl_EvaluateCumulativeSignal2Noise(Background(:,i),Signal(:,i),i); %The signal to noise peak anywhere between the 81st to 100th percentile
    if Manual==1
        txt=sprintf('%s%d%s','Do you want to remove channel ',i,'?, y/n');
        prompt=txt;
        a=input(prompt,'s');
        Bd=strcmp(a,'y');
        Gd=strcmp(a,'n');
        if Bd==1
            BadChannels(cBad)=i;
            cBad=cBad+1;
        elseif Gd==1
            GoodChannels(cGood)=i;
            cGood=cGood+1;
        else
            disp('Wrong Input Detected - one chance to correct it')
            txt=sprintf('%s%d%s','Do you want to remove channel ',i,'?, y/n');
            prompt=txt;
            a=input(prompt,'s');
            Bd=strcmp(a,'y');
            Gd=strcmp(a,'n');
            if Bd==1
                BadChannels(cBad)=i;
                cBad=cBad+1;
            elseif Gd==1
                GoodChannels(cGood)=i;
                cGood=cGood+1;
            end
        end
    else
        %Automated Detection
        if Sig2NoisePeak(i)>=Sig2NoiseThresh
            GoodChannels(cGood)=i;
            cGood=cGood+1;
        else
            BadChannels(cBad)=i;
            cBad=cBad+1;
        end
    end 
end
figure('Name','All the Signal to Noise Peaks')
bar(Sig2NoisePeak)
xlabel('Channel Number')
ylabel('Signal to Noise')
end