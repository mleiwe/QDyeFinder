function []=mnl_TracesChannelMIPs_TracesAndSomas(dim,efPxTrace,Somas,BkgMean,MaxValues)
%function that creates normalised MIPs for each channel based on the values
%from the trace extraction
% Inputs
%  dim - dimensions of the original image [x y c z]
%  data - the image matrix
%  BkgMean - the background values to subtract for each channel
%  MaxValues - the maximum values to normalise to

nTrace=size(efPxTrace,2);
nSoma=size(Somas,2);
%Per Channel
for i=1:dim(3)
    BW=zeros(dim(2),dim(1)); %Create a BW mask
    AllVx=[]; %The list of the key xy positions
    %Per Trace - get the important voxels
    for j=1:nTrace
        tAllVx=efPxTrace(j).AllVoxels(:,1:2);
        fAllVx=unique(tAllVx,'rows');
        tfAllVx=[AllVx;fAllVx];
        AllVx=unique(tfAllVx,'rows');
        clear tAllVx fAllVx tfAllVx
    end
    %Per Soma - get the Soma voxels
    for j=1:nSoma
        tAllVx=Somas(j).VxList(:,1:2);
        fAllVx=unique(tAllVx,'rows');
        tfAllVx=[AllVx;fAllVx];
        AllVx=unique(tfAllVx,'rows');
        clear tAllVx fAllVx tfAllVx
    end
    numVx=size(AllVx,1);
    %Now get the max value at each of the points
    for j=1:numVx
        BW(AllVx(j,2),AllVx(j,1))=max(data(AllVx(j,2),AllVx(j,1),i,:));
    end
    %Make a Tiff
    fn=sprintf('%s%d%s','TracesMIP_Raw_Channel',i,'.tiff');
    imwrite(uint16(BW),fn);
end
end