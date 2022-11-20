function [Somas,BkgMean,BkgStd,Thresh]=mnl_ProcessTracesFromNeurolucida_Somas(fname,cData,Scale,dim)
% Process the traces that have been extracted from Neurolucida and then get
% the corresponding raw voxel values for each trace (NB Trace=the
% neurolucida identified fragment/branch/axon. This version also splits the
% image into 2 along the y axis. But this can be removed quite easily.
%
% Inputs
% fname - The name of the xml document (NB the .xml portion is not
% necessary)
% cData - The chromatically corrected and linearly unmixed 4D image that
% was used to trace the data
% Scale - the scale of the x*y*z dimensions of the image
% dim - The dimensions of the image
%
% Outputs
% Somas - Structured Array Containing...
%       Points - List of the Points (in um)
%       Diameter - Diameter of the trace at that point
%       OriginalTrace - The number of the Original Trace
%       LengthVx - The Length of the Process in Voxels
%       LengthReal - The length is real distances
%       RawColourValues - List of the raw colour values of the assigned voxels
%       NormColourValues - The colour values that have been normalised to a bespoke value depending on the noise
%       NormColourThreshold - The Threshold value of the background used
%       BkgMean - The mean background colour. NB if values are below zero they get re-coded to zero 
%
% What to do next? Process all the different images and their traces. Then
% go to mnl_RecordTheColourSpreadOfTraces2 which will combine the different
% Traces into a single structured matrix
%Marcus Leiwe, Kyushu University - 29th June 2020
%% Step 1 - Read from .XML file
%[Trace,Somas]=mnl_ReadXML_DendritesAndSomas(fname);
[Somas]=mnl_ReadXML_Somas_v2(fname);
%% Step 2 -  Filter to Remove channels
[BadChannels,SampleBkgVxList,SampleVx]=mnl_FilterChannels_Somas(Somas,Data,Scale);
%Remove Channels
if isempty(BadChannels)==0
    AllChans=1:dim(3);
    GoodChans=find(AllChans~=BadChannels);
    fcData=cData(:,:,GoodChans,:);
    dim=size(fcData);
    cData=fcData;
    clear fcData
end
%Calculate Mean Background
Th=0;
[BkgMean,BkgStd,Thresh]=mnl_CalculateBackgroundValules(SampleBkgVxList,cData,Th);


%% Step 6 - Get the values for the somas
[Somas]=mnl_CalculateSomaValues(Somas,Scale,cData,BkgMean);
%% Save key variables
save('ProcessedSomas.mat','Somas','Scale','BkgMean','dim','BkgStd','Thresh')
%Evaluate Tetbow
mnl_VisualiseSomasDiscriminable(Somas); %Percent Discriminable
[Dimension]=mnl_BootstrapSomas(Somas); %Percent Unique

end

