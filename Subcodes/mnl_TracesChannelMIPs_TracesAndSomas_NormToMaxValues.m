function []=mnl_TracesChannelMIPs_TracesAndSomas_NormToMaxValues(dim,data,BkgMean,MaxValues)
%function that creates normalised MIPs for each channel based on the values
%from the trace extraction
% Inputs
%  dim - dimensions of the original image [x y c z]
%  data - the image matrix
%  BkgMean - the background values to subtract for each channel
%  MaxValues - the maximum values to normalise to

%% Code
data=double(data);
for i=1:dim(3)
    %% Step One - reduce to a per channel MIP
    temp=zeros(dim(2),dim(1),dim(4));
    for j=1:dim(4)
        temp(:,:,j)=data(:,:,i,j);
    end
    tMIP=max(temp,[],3);
    %% Step Two - Subtrace the background value
    bMIP=tMIP-BkgMean(i);
    %% Step Three - Divided by the highest trace value
    nbMIP=bMIP./MaxValues(i);
    %scale with some wriggle room
    sf=(2^16)/4;
    snbMIP=nbMIP*sf;
    %% Step Fours - Make a Tiff
    fn=sprintf('%s%d%s','TraceNorm_MIP_Raw_Channel',i,'.tiff');
    imwrite(uint16(snbMIP),fn); 
end
end