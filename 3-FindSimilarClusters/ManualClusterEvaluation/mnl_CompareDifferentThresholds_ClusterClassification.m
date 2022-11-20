function mnl_CompareDifferentThresholds_ClusterClassification
%% Load in each csv
c=1;
next=2;
Quant=struct('Threshold',[],'DataMatrix',[],'Single_n',[],'Single_pc',[],'Multiple_n',[],'Multiple_pc',[],...
    'MultipleClose_n',[],'MultipleClose_pc',[],'Few_n',[],'Few_pc',[],'Scattered_n',[],'Scattered_pc',[]);
%Load in the evaluations one by one
while c<next
    fprintf('%s\n','Select file of classified clusters')
    filename=uigetfile;
    T=readtable(filename,'VariableNamingRule','preserve');
    %Determine the threshold
    prompt='What is the threshold used by the dCrawler?';
    UserInput=input(prompt,'s');
    Thd=str2num(UserInput);
    %Get the data from the table
    [DataMatrix]=mnl_NestedExtractValuesFromCSV(T);
    
    %Count the number of clusters in each class and do percentages too
    nGoodSingle=sum(~isnan(DataMatrix(:,1)));
    nGoodMultiple=sum(~isnan(DataMatrix(:,2)));
    nGoodMultipleClose=sum(~isnan(DataMatrix(:,3)));
    nFew=sum(~isnan(DataMatrix(:,4)));
    nScattered=sum(~isnan(DataMatrix(:,5)));
    
    TotalClust=nGoodSingle+nGoodMultiple+nGoodMultipleClose+nFew+nScattered;
    pc_GoodSingle=(nGoodSingle/TotalClust)*100;
    pc_GoodMultiple=(nGoodMultiple/TotalClust)*100;
    pc_GoodMultipleClose=(nGoodMultipleClose/TotalClust)*100;
    pc_Few=(nFew/TotalClust)*100;
    pc_Scattered=(nScattered/TotalClust)*100;

    %Save into structure
    Quant(c).Threshold=Thd;
    Quant(c).DataMatrix=DataMatrix;
    Quant(c).Single_n=nGoodSingle;
    Quant(c).Single_pc=pc_GoodSingle;
    Quant(c).Multiple_n=nGoodMultiple;
    Quant(c).Multiple_pc=pc_GoodMultiple;
    Quant(c).MultipleClose_n=nGoodMultipleClose;
    Quant(c).MultipleClose_pc=pc_GoodMultipleClose;
    Quant(c).Few_n=nFew;
    Quant(c).Few_pc=pc_Few;
    Quant(c).Scattered_n=nScattered;
    Quant(c).Scattered_pc=pc_Scattered;

    %Finished?
    while 1
        prompt='Are there more Thresholds to load in? y/n';
        UserInput=input(prompt,'s');
        if strcmp(UserInput,'y')==1
            fprintf('%s\n','Okay next data set...')
            next=next+1;
            c=c+1;
            break
        elseif strcmp(UserInput,'n')==1
            c=c+1;
            break
        else
            fprintf('%s\n','Input not recognised please try again')
        end
    end
end

%% Make figures
nThresh=c-1;
for i=1:nThresh
    Singles(i)=Quant(i).Single_pc;
    Multiple(i)=Quant(i).Multiple_pc;
    MultipleClose(i)=Quant(i).MultipleClose_pc;
    Fews(i)=Quant(i).Few_pc;
    Scattereds(i)=Quant(i).Scattered_pc;
end
Multiples=MultipleClose+Multiple;

%Subplots of pie charts
figure('Name','Pie Charts')
cmap=flipud(magma(4));
colormap(cmap)
labs={'Good Single','Good Multiple','Few','Scattered'};
%labs={'Scattered','Few','Good Multiple','Good Single'};
t=tiledlayout('flow','TileSpacing','compact');
for i=1:nThresh
    PieVals=[Singles(i) Multiples(i) Fews(i) Scattereds(i)];
    %PieVals=[Scattereds(i) Fews(i) Multiples(i) Singles(i)];
    nexttile
    pie(PieVals)
    tn=num2str(round(Quant(i).Threshold,2));
    tn2=sprintf('%s%s','Th(d) = ',tn);
    title(tn2)
end
lgd=legend(labs);
lgd.Layout.Tile = nThresh+1;

%Stacked bars
for i=1:nThresh
    StackVals(i,:)=[Singles(i) Multiples(i) Fews(i) Scattereds(i)];   
    %StackVals(i,:)=[Scattereds(i) Fews(i) Singles(i) Multiples(i)];   
    xtcklab{i}=num2str(round(Quant(i).Threshold,3));
end
figure('Name','Stacked Bars')
colormap(cmap)
ba=bar(StackVals,'stacked','FaceColor','flat');
for i=1:4
    ba(i).CData=cmap(i,:);
end
xticks([1:nThresh])
xticklabels(xtcklab)
legend(labs,'Location','eastoutside')
xlim([1-0.75 nThresh+0.75])
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