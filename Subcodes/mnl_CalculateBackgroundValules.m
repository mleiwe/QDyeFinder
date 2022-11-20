function [mBkgVals,sBkgVals,Thresh]=mnl_CalculateBackgroundValules(SampleBkgVxList,cData,Th)
%Function to calculate the background values for an image
% Inputs
%   SampleBkgVxList - The list of x,y, and z points of the background
%   cData - The image (x*y*c*z)
%   Th - The number of Sd above the mean to use as the threshold value
%
% Outputs
%   mBkgVals - The mean value
%   Thresh - Threshold values to be used later
dim=size(cData);
szVx=size(SampleBkgVxList,1);
BkgVals=nan(szVx,dim(3));
for i=1:szVx
    tVals=cData(SampleBkgVxList(i,2),SampleBkgVxList(i,1),:,SampleBkgVxList(i,3)); %Remember it is flipped
    BkgVals(i,:)=tVals;
end
mBkgVals=mean(BkgVals);
sBkgVals=std(BkgVals);
Thresh=mBkgVals+(Th*sBkgVals);
end