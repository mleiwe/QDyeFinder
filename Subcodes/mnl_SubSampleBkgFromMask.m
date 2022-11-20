function [BkgVx]=mnl_SubSampleBkgFromMask(nBkg,dim,AllVx)
%Function to extract a set number(nBkg) of values that are not in the mask
%Inputs
%nBkg - The number of background voxels that you want
%dim - The dimensions of the image
%marked as 0
%AllVx - The list of voxels that are 1 in the BW
%Outputs
% BkgVx - The [x,y,z] components of the background selected
C=1; %C is the counter
szBW=dim;
BkgVx=nan(nBkg,3);
tic
while C<=nBkg
    tempXYZ=ceil(rand(1,3).*szBW);
    %Are the tempXYZ co-ordinates present in the mask
    [tf,~]=ismember(tempXYZ,AllVx,'rows');
    if tf==0       
        %Is the position already present in the BkgVx list
        if isempty(BkgVx)==1
            BkgVx=tempXYZ;
            mnl_InsertProgressTrackerInLoops(C,nBkg)
            C=C+1;
        else
            [tf,~]=ismember(tempXYZ,BkgVx,'rows');
            if tf==0
                BkgVx(C,:)=tempXYZ;
                mnl_InsertProgressTrackerInLoops(C,nBkg)
                C=C+1;
            end
        end    
    end
end
toc
end