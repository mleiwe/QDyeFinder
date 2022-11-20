function mnl_CombineManualClassificationGroups
%% Load in CSVs
disp('Please load in all the csv files')
[files]=uipickfiles;
%Pre-allocate groups for bars
Classes_Struct=struct('Thd',[],'DataMatrix',[], 'SingleCluster',[],'PartialCluster',[],'MultipleCluster',[], ...
    'Single_n',[],'Single_pc',[],'Partial_n',[],'Partial_pc',[],'Multiple_n',[],'Multiple_pc',[]);

%% Regroup each set
nfiles=length(files);
for i=1:nfiles
    %% Work out the threshold from the file name
    [path,name,extension]=fileparts(files{i});
    Thd_name=erase(name,'ClusterClassification_Thd_');
    %Now load the csv and convert to matrix
    filename=files{i};
    T=readtable(filename,'VariableNamingRule','preserve');
    [DataMatrix]=mnl_NestedExtractValuesFromCSV(T);
    %% Re-classify
    %Convert to individual variables
    %Get good single, good multiple, good anatomically close, few and scattered
    nMax=size(T,1);
    GoodSingle=DataMatrix(~isnan(DataMatrix(:,1)),1);
    GoodMultiple=DataMatrix(~isnan(DataMatrix(:,2)),2);
    GoodMultipleClose=DataMatrix(~isnan(DataMatrix(:,3)),3);
    Partial=DataMatrix(~isnan(DataMatrix(:,4)),4);
    Scattered=DataMatrix(~isnan(DataMatrix(:,5)),5);
    %Now recode
    Multiple=[GoodMultiple;GoodMultipleClose;Scattered];
    if length(Multiple)>nMax
        nMax=length(Multiple);
    end
    %% Export as CSV
    %Calculate the longest string
    nSingle=length(GoodSingle);
    nMultiple=length(Multiple);
    nPartial=length(Partial);
    nClust=sum([nSingle,nMultiple,nPartial]);
    %Output Matrix
    OutputMatrix=nan(nMax,3);
    OutputMatrix(1:nSingle,1)=GoodSingle;
    OutputMatrix(1:nPartial,2)=Partial;
    OutputMatrix(1:nMultiple,3)=Multiple;
    %Make a table
    OutputTable=table(OutputMatrix(:,1),OutputMatrix(:,2),OutputMatrix(:,3),'VariableNames',{'Good single cluster','Partial cluster','Multiple cluster'});
    %Write to csv
    csv_fn=sprintf('%s%s%s','ClusterClassification_ThreeGroups_Thd_',Thd_name,'.csv');
    writetable(OutputTable,csv_fn);
    %% Assign to structure
    % 'Single_n',[],'Single_pc',[],'Partial_n',[],'Partial_pc',[],'Multiple_n',[],'Multiple_pc',[]);
    Classes_Struct(i).Thd=Thd_name;
    Classes_Struct(i).DataMatrix=OutputMatrix;
    %Single
    Classes_Struct(i).SingleCluster=GoodSingle;
    Classes_Struct(i).Single_n=nSingle;
    Classes_Struct(i).Single_pc=(nSingle/nClust)*100;
    %Partial
    Classes_Struct(i).PartialCluster=Partial;
    Classes_Struct(i).Partial_n=nPartial;
    Classes_Struct(i).Partial_pc=(nPartial/nClust)*100;
    %Multiple
    Classes_Struct(i).MultipleCluster=Multiple;
    Classes_Struct(i).Multiple_n=nMultiple;
    Classes_Struct(i).Multiple_pc=(nMultiple/nClust)*100;
end
%% Plot figures
for i=1:nfiles
    Singles(i)=Classes_Struct(i).Single_pc;
    Partials(i)=Classes_Struct(i).Partial_pc;
    Multiples(i)=Classes_Struct(i).Multiple_pc;
end
%Subplots of pie charts
figure('Name','Pie Charts')
cmap=[1 1 1;0.5 0.5 0.5;0 0 0];
colormap(cmap)
labs={'Single','Partial','Multiple'};
t=tiledlayout('flow','TileSpacing','compact');
for i=1:nfiles
    PieVals=[Singles(i) Partials(i) Multiples(i)];    
    nexttile
    pie(PieVals)
    tn=Classes_Struct(i).Thd;
    title(tn)
end
lgd=legend(labs);
lgd.Layout.Tile = nfiles+1;

%Stacked bars
StackVals=nan(nfiles,3);
for i=1:nfiles
    StackVals(i,:)=[Singles(i) Partials(i) Multiples(i)];   
    xtcklab{i}=Classes_Struct(i).Thd;
end
figure('Name','Stacked Bars')
colormap(cmap)
ba=bar(StackVals,'stacked','FaceColor','flat');
for i=1:3
    ba(i).CData=cmap(i,:);
end
xticks([1:nfiles])
xticklabels(xtcklab)
legend(labs,'Location','eastoutside')
xlim([1-0.75 nfiles+0.75])
ylim([0 100])
xlabel('Th(d)')
ylabel('% of clusters')
end
function [DataMatrix]=mnl_NestedExtractValuesFromCSV(T)
%Get good single, good multiple, good anatomically close, few and scattered
%and put into columns
nMax=size(T,1);
GoodSingle=T(:,'Good Single Cluster');
GoodMultiple=T(:,'Good Multiple Cluster');
GoodMultipleClose=T(:,'Good but anatomical close');
Few=T(:,'Few');
Scattered=T(:,'Scattered');

DataMatrix=nan(nMax,5); %Pre-allocatie DataMatrux
%% Good single
ClassClusters=GoodSingle{:,1};
idx=~isnan(ClassClusters);
ClassClusters=ClassClusters(idx);
%Find corresponding traces and nTraces
nClust=length(ClassClusters);
%Save into data structure
DataMatrix(1:nClust,1)=ClassClusters;

%% Good multiple
ClassClusters=GoodMultiple{:,1};
idx=~isnan(ClassClusters);
ClassClusters=ClassClusters(idx);
%Find corresponding traces and nTraces
nClust=length(ClassClusters);
%Save into data structure
DataMatrix(1:nClust,2)=ClassClusters;

%% Good multiple close
ClassClusters=GoodMultipleClose{:,1};
idx=~isnan(ClassClusters);
ClassClusters=ClassClusters(idx);
%Find corresponding traces and nTraces
nClust=length(ClassClusters);
%Save into data structure
DataMatrix(1:nClust,3)=ClassClusters;

%% Few
ClassClusters=Few{:,1};
idx=~isnan(ClassClusters);
ClassClusters=ClassClusters(idx);
%Find corresponding traces and nTraces
nClust=length(ClassClusters);
%Save into data structure
DataMatrix(1:nClust,4)=ClassClusters;

%% Scattered
ClassClusters=Scattered{:,1};
idx=~isnan(ClassClusters);
ClassClusters=ClassClusters(idx);
%Find corresponding traces and nTraces
nClust=length(ClassClusters);
%Save into data structure
DataMatrix(1:nClust,5)=ClassClusters;
end