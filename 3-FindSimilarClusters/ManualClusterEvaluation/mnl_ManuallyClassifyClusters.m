function [GoodSingle,GoodMultiple,GoodAnatomicallyClose,Few,Scattered]=mnl_ManuallyClassifyClusters
%This function manually classifies the clusters. 
% The input is manually chosen tiff images of the clusters
% The output are lists of the clusters numbers and saved as both a
%csv, and as separate output variables

%% List of tiffs to loads
[files]=uipickfiles;
nFiles=size(files,2);
%Pre-allocate groups
GoodSingle=[];
GoodMultiple=[];
GoodAnatomicallyClose=[];
Few=[];
Scattered=[];
%Initialise figure
figure('Name','Select class... 1=Good Single, 2= Good Multiple, 3=Anatomically close, 4= Few, and 5=Scattered','Units','normalized','Position',[0.005 0.4 0.5 0.5])
%% Loop through each image
for i=1:nFiles
    %Get the file name
    file_full=files{i};
    [~,fn,~]=fileparts(file_full);
    ClustNum=regexp(fn,'\d*','Match');
    %Load the image
    Im=imread(file_full);
    image(Im)
    axis equal
    axis off
    tn=sprintf('%s%s','Cluster ',ClustNum{1});
    title(tn)
    PcComp=round(((i-1)/nFiles)*100,3);
    stn=sprintf('%d%s',PcComp,'% complete');
    subtitle(stn)
    %Classify the neuron
    while 1
        prompt='Select class... 1=Good Single, 2= Good Multiple, 3=Anatomically close, 4= Few, and 5=Scatttered';
        UserInput=input(prompt,'s');
        if strcmp(UserInput,'1')==1
            GoodSingle=[GoodSingle;ClustNum];
            break
        elseif strcmp(UserInput,'2')==1
            GoodMultiple=[GoodMultiple;ClustNum];
            break
        elseif strcmp(UserInput,'4')==1
            Few=[Few;ClustNum];
            break
        elseif strcmp(UserInput,'3')==1
            GoodAnatomicallyClose=[GoodAnatomicallyClose;ClustNum];
            break
        elseif strcmp(UserInput,'5')==1
            Scattered=[Scattered;ClustNum];
            break
        else
            disp('Input not recognised please try again')
        end
    end
end
%% Output as a CSV
%Convert to double
GoodSingle=str2double(GoodSingle);
GoodMultiple=str2double(GoodMultiple);
GoodAnatomicallyClose=str2double(GoodAnatomicallyClose);
Few=str2double(Few);
Scattered=str2double(Scattered);

%Calculate the longest string
nSingle=length(GoodSingle);
nMultiple=length(GoodMultiple);
nFew=length(Few);
nGoodAnatomicallyClose=length(GoodAnatomicallyClose);
nScattered=length(Scattered);
MaxN=max([nSingle,nMultiple,nFew,nGoodAnatomicallyClose,nScattered]);
%Output Matrix
OutputMatrix=nan(MaxN,5);
OutputMatrix(1:nSingle,1)=GoodSingle;
OutputMatrix(1:nMultiple,2)=GoodMultiple;
OutputMatrix(1:nGoodAnatomicallyClose,3)=GoodAnatomicallyClose;
OutputMatrix(1:nFew,4)=Few;
OutputMatrix(1:nScattered,5)=Scattered;
%Make a table
OutputTable=table(OutputMatrix(:,1),OutputMatrix(:,2),OutputMatrix(:,3),OutputMatrix(:,4),OutputMatrix(:,5),'VariableNames',{'Good Single Cluster','Good Multiple Cluster','Good but anatomical close','Few','Scattered'});
%Write to csv
writetable(OutputTable,'ClusterClassification.csv')




