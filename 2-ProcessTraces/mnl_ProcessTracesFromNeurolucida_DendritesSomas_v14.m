function [efPxTrace,Somas,BkgSigCumulativeDistribution,NumLengthThresh,MagThresh,BkgMedian,MaxValues]=mnl_ProcessTracesFromNeurolucida_DendritesSomas_v14(fname,cData,Scale,dim)
% Process the traces that have been extracted from Neurolucida and then get
% the corresponding raw voxel values for each trace (NB Trace=the
% neurolucida identified fragment/branch/axon.
%
% Version 14 follows version 13, but the difference is now that any nans in the image are now ignored in the calculation of the means.
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
% efPxTrace - Structured matrix containing...
%       Points - List of the Points
%       Diameter - Diameter of the trace at that point
%       OriginalTrace - The number of the Original Trace
%       LengthVx - The Length of the Process in Voxels
%       LengthReal - The length is real distances
%       AllVoxels - All the voxels for this trace
%       SetId - A group identifier that is specified in Neurolucida. e.g. Cell 1
%       TypeId - What type of trace it is, i.e. Apical Dendtrite, Dendrite, Axon, etc... 
%       Scale - The size of each voxel
%       dim - The dimensions of the image
%       RawMean - List of the raw colour values of the assigned voxels
%       BkgRemovedMean - The mean colour after background subtraction
% Somas - Structured matrix containing...
%       dxyz - the diameter, x,y, and z positions in um NB not in Vx
%       SetId - A group identifier that is specified in Neurolucida. e.g. Cell 1
%       VxList - All the voxels for this trace
%       Scale - The size of each voxel
%       dim - The dimensions of the image
%       RawMean - List of the raw colour values of the assigned voxels
%       BkgSubMean - The mean colour after background subtraction
% BkgSigCumulativeDistribution - Structured matrix containing ...
%       BackgroundPercentiles - n*2 matrix containing the percentile value
%       (1st column), and the intensity value (2nd column) for the
%       background voxels
%       SignalPercentiles - n*2 matrix containing the percentile value
%       (1st column), and the intensity value (2nd column) for the
%       mean intensity of each fragment
%       NB this is sorted per channel
% NumLengthThresh - The length threshold at which the each trace must be longer
% MagThresh - The magnitude threshold
% BkgMedian- The median background value in all channels
% MaxValues - The maximum value across all traces
%
% What to do next? Process all the different images and their traces. Then
% go to mnl_RecordTheColourSpreadOfTraces2 which will combine the different
% Traces into a single structured matrix
%Marcus Leiwe, Kyushu University - 3rd February 2022
%% Pre-decisions
FilteredPlot=0;
VoxelPlot=0;
EvaluateTracesPerCell='y';
SingleChannelMIPs='n';
%% Step 1 - Read from .XML file
[Trace,Somas]=mnl_ReadXML_DendritesAndSomas4(fname);
%Correct Neurolucida tracing
prompt='Do you want to correct the automatic Neurolucida tracing? (y/n)';
UserInput=input(prompt,'s');
TraceEval=strcmp(UserInput,'y');
if TraceEval==1
    [Trace]=mnl_AdjustTraceLocations(Trace,cData,Scale);
end
%% Step 2 -  Filter to Remove channels
[BadChannels,SampleBkgVxList,~,BkgSigCumulativeDistribution]=mnl_FilterChannels_FragmentBased(Trace,Scale,cData);
%Remove Channels
if isempty(BadChannels)==0
    AllChans=1:dim(3);
    GoodChans=find(AllChans~=BadChannels);
    fcData=cData(:,:,GoodChans,:);
    dim=size(fcData);
    cData=fcData;
    clear fcData
end
close all
%%  Step 3 - Filter to remove bad traces
[Trace]=mnl_CalculateLengths(Trace,Scale);% Measure Length of Traces
NumTraceLimit=size(Trace,2);
%Step 3a - The minimum fragment length for colour consistency
prompt='Do you want to calculate the minimum fragment length for colour consistency? (y/n)';
UserInput=input(prompt,'s');
ColourEval=strcmp(UserInput,'y');
if ColourEval==1
    [FragmentValues,BkgMedian,MaxValues]=mnl_SplitTracesIntoFragments_v2(Trace,Scale,cData,SampleBkgVxList); %Make all the subfragments and measure their colour consistency and also their brightness magnitude
    ColourConsistencyThreshold=0.075;
    [NumLengthThresh]=mnl_EvaluateFragmentsColourConsistency5(FragmentValues,ColourConsistencyThreshold);%to calculate the length when the colours becomes consistent
    lim=round(NumLengthThresh); %minimum distance (in um) required to be a true trace
    fprintf('%s%d\n','The minimum length is ',NumLengthThresh)
else
    prompt='Then please insert the minimum length of a trace';
    lim=input(prompt,'s');
    lim=str2double(lim);
    NumLengthThresh=lim;
end
%Step 3b - Filter for Brightness
prompt='Do you want to determine the minimum brightness of the trace? (y/n)';
UserInput=input(prompt,'s');
if strcmp(UserInput,'y')==1
    %Magnitude threshold settings
    BrightnessEuclideanThreshold=0.2;
    ErrorThresh=3;
    FragLength=5;
    if FragLength<lim
        FragLength=lim;
    end
    %Split into fragments
    SubFrag='y';
    [FragmentValues,BkgMedian,MaxValues]=mnl_CalculateMaxValues_v3(Trace,Scale,cData,SampleBkgVxList,FragLength,SubFrag);
    % Now determing the magnitude threshold
    [MagThresh]=mnl_DetermineBrightnessEffects_v4(FragmentValues,NumLengthThresh,BrightnessEuclideanThreshold,ErrorThresh);
    MagnitudeThresh=MagThresh;
    fprintf('%s%d\n','The minimum brightness magnitude is',MagThresh)
else
    prompt='Then please input the desired magnitude threshold';
    UserInput=input(prompt,'s');
    MagThresh=str2double(UserInput);
    MagnitudeThresh=MagThresh;
    %Split the traces for future calculations
    FragLength=5;
    if FragLength<lim
        FragLength=lim; %Changes the fragment length to the length threshold if the length is greater than 5
    end
    %Split into fragments
    SubFrag='y';
    [FragmentValues,BkgMedian,MaxValues]=mnl_CalculateMaxValues_v3(Trace,Scale,cData,SampleBkgVxList,FragLength,SubFrag);
end
%Step 3c - Now remove the colour inconsistent long traces
LengthThresh=1000; %The length at which to begin checking
% LengthThresh=30; %The length at which to begin checking
InconsistentColourThresh=0.5;

[fTrace]=mnl_DoubleCheckLongTraceConsistency_v3(Trace,lim,LengthThresh,cData,BkgMedian,MaxValues,Scale,InconsistentColourThresh);

% Now visualise the Filtered Traces
szT=size(fTrace,2);
if FilteredPlot==1
    cmap=colormap(lines(szT));
    figure('Name','fTrace Locations')
    for i=1:szT
        scatter3(fTrace(i).RealPoints(:,1),fTrace(i).RealPoints(:,2),fTrace(i).RealPoints(:,3),'.','MarkerFaceColor',cmap(i,:))
        hold on
    end
    axis equal
end
%%  Step 4 - Create Mask
[efPxTrace,fPxTrace]=mnl_MakeSeperateMasksForEachImage(fTrace,Scale,dim);
% Now visualise the traces in voxel dimensions
if VoxelPlot==1
    szT=size(fPxTrace,2);
    cmap=colormap(lines(szT));
    figure('Name','fPxTrace Locations')
    for i=1:szT
        scatter3(fPxTrace(i).Points(:,1),fPxTrace(i).Points(:,2),fPxTrace(i).Points(:,3),'.','MarkerFaceColor',cmap(i,:))
        hold on
    end
    axis equal
    % Now visualise the expanded areas going to be used to get colour values
    szT=size(efPxTrace,2);
    cmap=colormap(lines(szT));
    figure('Name','efPxTrace -Voxels')
    for i=1:szT
        scatter3(efPxTrace(i).AllVoxels(:,1),efPxTrace(i).AllVoxels(:,2),efPxTrace(i).AllVoxels(:,3),'.','MarkerFaceColor',cmap(i,:))
        hold on
    end
    axis equal
    axis ij
end
%Store Scale and dim information
for i=1:szT
    efPxTrace(i).Scale=Scale; %Preallocate the scale of the whole image to the first array
    efPxTrace(i).dim=dim; %Preallocate the dimensions of the whole image to the first array
end
%%  Step 5 - Now Get the Colour Values
[BkgMedian,MaxValues,efPxTrace]=mnl_ExtractVoxelValuesForTraces_v7(efPxTrace,cData,BkgMedian,MagThresh,NumLengthThresh); %Calculate for Full Image
%% Step 6 - Get the values for the somas too
if isempty(Somas)==0
    [Somas]=mnl_CalculateSomaValues_v2(Somas,Scale,cData,BkgMedian);
end
%% Plot a final figure of processed traces
prompt='Do you want a figure of the final processed traces? (y/n)';
UserInput=input(prompt,'s');
FigYN=strcmp(UserInput,'y');
if FigYN==1
    mnl_PlotFinalTracingFigure(efPxTrace,Somas,dim,Scale)
end
%% Now Create Single Channel MIPs of just the trace voxels
if strcmp(SingleChannelMIPs,'y')==1
    %Max Trace Value Normalised (full image)
    mnl_TracesChannelMIPs_TracesAndSomas_NormToMaxValues(dim,cData,BkgMedian,MaxValues)
    %Max Voxel Normalised
    mnl_TracesChannelMIPs_TracesAndSomas(dim,efPxTrace,Somas,cData);
    mnl_TracesChannelMIPs_Processed_TracesAndSomas(dim,efPxTrace,Somas)
    
end
%% Step 7 - Evaluate the Distribution
if strcmp(EvaluateTracesPerCell,'y')==1
    %Check if the SetId is specified
    szT=size(efPxTrace,2);
    Sets=zeros(szT,1);
    for i=1:szT
        Sets(i)=~isempty(efPxTrace(i).SetId);
    end
    if sum(Sets)>=1
        mnl_MeasureTheColourVariableOfCells_MagnitudeWeighted(efPxTrace,'n');
    else
        disp("No SetIds specified skipping the step")
    end
end
% Now insert a recommended EuThresh
if ~isempty(Somas)
    [ColourDiversity_95]=mnl_EvaluateColourDistribution(dim(3),Somas);
end

%% Step 8 - Save the key Information
save('ProcessedTraces.mat','efPxTrace','Somas','NumLengthThresh','BkgMedian','MagnitudeThresh','dim','Scale','MaxValues')
end
%% Subfunctions
function mnl_PlotFinalTracingFigure(efPxTrace,Somas,dim,Scale)
h=figure('Name','The Final Processed Traces');
nSoma=size(Somas,2);
nTraces=size(efPxTrace,2);
TraceSetId=0;
if isempty(Somas)==0
    [cmap]=mnl_GenerateShuffledColourmap(nSoma);
    for i=1:nSoma
        if isempty(Somas(i).dxyz)==0
            SetId=Somas(i).SetId;
            %Plot the soma
            Points=Somas(i).dxyz(:,2:3)./Scale(1);
            Points(:,2)=Points(:,2)*-1; %Y axis needs to be inervted
            scatter(Points(:,1),Points(:,2),'.','MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:))
            hold on
            %Find and draw the matching traces
            for j=1:nTraces
                chk=strcmp(SetId,efPxTrace(j).SetId);
                if chk==1
                    TraceSetId=1;
                    TracePoints=efPxTrace(j).Points(:,1:2);
                    Diameters=efPxTrace(j).Diameter./Scale(1);
                    mnl_PlotSingleTrace2D(TracePoints,Diameters,cmap(i,:))
                end
            end
        end
    end
    %Now draw the unattached Traces
    if TraceSetId==0
        [cmap]=mnl_GenerateShuffledColourmap(nTraces);
    else
        [cmap]=zeros(nTraces,3);
    end
    for i=1:nTraces
        if isempty(efPxTrace(i).SetId)==1
            TracePoints=efPxTrace(i).Points(:,1:2);
            Diameters=efPxTrace(i).Diameter./Scale(1);
            mnl_PlotSingleTrace2D(TracePoints,Diameters,cmap(i,:))
            hold on
        end
    end
else
    %Now draw the unattached Traces
    for i=1:nTraces
        if isempty(efPxTrace(i).SetId)==1
            TracePoints=efPxTrace(i).Points(:,1:2);
            Diameters=efPxTrace(i).Diameter./Scale(1);
            mnl_PlotSingleTrace2D(TracePoints,Diameters,[0 0 0])
            hold on
        end
    end
end
axis equal
axis ij
xlim([0 dim(1)])
ylim([0 dim(2)])
grid on
%mnl_ExportEPSdense(h,'ProcessedTracesAndSomas')
savefig(h,'ProcessedTracesAndSomas')
end
function mnl_PlotSingleTrace2D(TracePoints,Diameters,Colour)
%Inputs
% TracePoints - The position of the points
% Diameters - The diameter of the trace
% Colour - the RGB values that you want to plot the trace in
%% Now make the volume for each trace
PointList=TracePoints(:,1:2);
szD=size(Diameters);
[nD,p]=max(szD);
%Work if the diameters are encoded in rows (p=1) or columns (p=2)
if p==1
    for j=1:nD
        Diameter(j,1)=Diameters(j,1);
    end
elseif p==2
    for j=1:nD
        Diameter(j,1)=Diameters(1,j);
    end
end
%Left Side
Xleft=PointList(:,1)-(Diameter/2);
%Right Side
Xright=PointList(:,1)+(Diameter/2);
X=[Xleft;flipud(Xright)];
%Sort Out Y
Y=[PointList(:,2);flipud(PointList(:,2))];
%Merge
patch(X,Y,Colour,'EdgeColor','none')
end