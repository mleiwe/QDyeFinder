function []=mnl_MakeMaxNormalisedMIPs(Data,n)
sz=size(Data);
NumChan=sz(3);
NumStacks=sz(4);
NumVx=sz(1)*sz(2)*NumStacks;

for i=1:NumChan
    %Make a single z stack for the colour
    tData=nan(sz(1),sz(2),sz(4));
    for j=1:sz(4)
        tData(:,:,j)=Data(:,:,i,j);
    end
    %Find the Nth percentile
    Data2=reshape(tData,[NumVx,1]);
    NormValue=prctile(Data2,n);
    %Make a MIP
    for j=1:sz(1)
        for k=1:sz(2)
            tMIP(j,k)=max(tData(j,k,:));
        end
    end    
    %Normalise to Nth Percentile
    tNormMIP=tMIP/NormValue;
    %If image goes above 1 make it 1
    idx=tNormMIP>1;
    tNormMIP(idx)=1;
    %Assign to Colour Stack
    NormMIPs(:,:,i)=tNormMIP;
    %Now Save the Image
    fn=sprintf('%s%d%s%d%s','NormMIPto',n,'percentile_Channel',i,'.tiff');
    imwrite(tNormMIP,fn);
end
%% Make a figure
figure('Name','Normalised MIPs')
cmap=magma(256);
MxVal=max(NormMIPs(:));
TotalMax=0;
NumR=round(sqrt(NumChan+1));
NumC=ceil(sqrt(NumChan+1));
colormap(cmap)
for i=1:NumChan+1
    if i<=NumChan
        st=sprintf('%s%d','Channel ',i);
        subplot(NumR,NumC,i)
        %imagesc(NormMIPs(:,:,i),[0 MxVal])
        imagesc(NormMIPs(:,:,i))
        colorbar
        axis equal
        axis off
        title(st)
    else
        subplot(NumR,NumC,i)
        imagesc(NormMIPs(:,:,i-1),[0 MxVal])
        colorbar
    end
end